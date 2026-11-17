"""按读段号去重:逐对流式处理,读段号索引放在磁盘 SQLite 里。

规则:
- 按输入顺序处理读段对,首次出现的读段号原样保留。
- 同一读段号再次出现,且 R1/R2 的序列与质量都完全相同 -> 剔除。
- 同一读段号再次出现但内容不同 -> 抛 DuplicateContentConflict,整次去重失败。
- R1/R2 必须逐对配对(读段号一致、记录数一致),总对数必须等于原入库对数。

去重键写入调用方指定的 SQLite 文件,页缓存有上限,不把全部读段号放在内存里。
"""
from __future__ import annotations

import hashlib
import os
import sqlite3
from itertools import zip_longest

from .demux import read_id

_MISSING = object()
_COMMIT_EVERY = 100_000


class DedupError(Exception):
    """去重失败。上层据此把任务标为 failed,不发布任何文件或计数。"""


class FastqFormatError(DedupError):
    """源 FASTQ 结构损坏。"""


class PairingError(DedupError):
    """R1/R2 记录数或读段号不能逐对配对。"""


class PairCountMismatch(DedupError):
    """读段对数与原入库对数不一致。"""


class DuplicateContentConflict(DedupError):
    """同一读段号出现多次,且内容并不相同。"""


class SourceChangedError(DedupError):
    """源文件在处理期间发生变化。"""


class _HashedFastq:
    """按 4 行一条读取未压缩 FASTQ,同时计算所读原始字节的 SHA-256。

    读盘是二进制的,因此哈希覆盖换行符等原始字节;用于处理结束后复核源文件
    是否仍是同一份内容。结构校验规则与拆样流程保持一致。
    """

    def __init__(self, path, label):
        self.path = path
        self.label = label
        self.hasher = hashlib.sha256()
        self.records = 0

    def __iter__(self):
        with open(self.path, "rb") as fh:
            while True:
                block = [fh.readline() for _ in range(4)]
                if block[0] == b"":
                    return  # 恰好结束在记录边界
                self.records += 1
                for line in block:
                    self.hasher.update(line)
                if block[1] == b"" or block[2] == b"" or block[3] == b"":
                    raise FastqFormatError(
                        f"{self.label}: 第 {self.records} 条记录在文件末尾被截断"
                    )
                try:
                    lines = [line.decode("ascii").rstrip("\r\n") for line in block]
                except UnicodeDecodeError as exc:
                    raise FastqFormatError(
                        f"{self.label}: 第 {self.records} 条附近不是纯文本"
                        f"(需要未压缩 FASTQ): {exc}"
                    ) from exc
                header, seq, plus, qual = lines
                if not header.startswith("@"):
                    raise FastqFormatError(
                        f"{self.label}: 第 {self.records} 条记录头缺少 '@'"
                    )
                if not plus.startswith("+"):
                    raise FastqFormatError(
                        f"{self.label}: 第 {self.records} 条记录分隔行缺少 '+'"
                    )
                if len(seq) != len(qual):
                    raise FastqFormatError(
                        f"{self.label}: 第 {self.records} 条({header})序列长度 "
                        f"{len(seq)} 与质量长度 {len(qual)} 不一致"
                    )
                yield header[1:], seq, plus, qual


def _signature(seq1, qual1, seq2, qual2) -> bytes:
    """一对读段的内容指纹:只含 R1/R2 的序列与质量,带长度前缀避免拼接歧义。"""
    digest = hashlib.sha256()
    for part in (seq1, qual1, seq2, qual2):
        raw = part.encode("ascii")
        digest.update(len(raw).to_bytes(8, "big"))
        digest.update(raw)
    return digest.digest()


def deduplicate(r1_path, r2_path, expected_pairs, out_r1, out_r2, seen_db_path):
    """逐对去重并写入 out_r1/out_r2;seen_db_path 是磁盘索引数据库。

    返回 {'input_pairs', 'kept_pairs', 'dropped_pairs', 'r1_sha256', 'r2_sha256'}。
    任何异常都不会自行删除半成品,由上层统一清理,保证不发布部分结果。
    """
    if expected_pairs < 0:
        raise ValueError("expected_pairs 不能为负")
    stream1 = _HashedFastq(r1_path, "R1")
    stream2 = _HashedFastq(r2_path, "R2")
    conn = sqlite3.connect(str(seen_db_path))
    total = kept = dropped = 0
    try:
        # 索引数据库只服务于本次任务,崩溃后可整体丢弃,因此不必保证其自身持久性;
        # 页缓存限制在 16 MiB,大样本的读段号由 SQLite 落到磁盘。
        conn.execute("PRAGMA journal_mode = OFF")
        conn.execute("PRAGMA synchronous = OFF")
        conn.execute("PRAGMA temp_store = FILE")
        conn.execute("PRAGMA cache_size = -16384")
        # 同一暂存目录里若残留上次的索引(例如失败重试),必须以空索引重新开始
        conn.execute("DROP TABLE IF EXISTS seen")
        conn.execute(
            "CREATE TABLE seen (rid TEXT PRIMARY KEY, sig BLOB NOT NULL) WITHOUT ROWID"
        )
        conn.commit()
        with open(out_r1, "w", encoding="ascii", newline="\n",
                  buffering=1 << 20) as w1, \
             open(out_r2, "w", encoding="ascii", newline="\n",
                  buffering=1 << 20) as w2:
            for rec1, rec2 in zip_longest(stream1, stream2, fillvalue=_MISSING):
                if rec1 is _MISSING or rec2 is _MISSING:
                    raise PairingError("R1 与 R2 记录数不一致,双端无法配对")
                total += 1
                if total > expected_pairs:
                    raise PairCountMismatch(
                        f"读段对数超过原入库对数 {expected_pairs},源文件可能已变化"
                    )
                h1, s1, p1, q1 = rec1
                h2, s2, p2, q2 = rec2
                rid1, rid2 = read_id(h1), read_id(h2)
                if not rid1 or rid1 != rid2:
                    raise PairingError(
                        f"第 {total} 对读段号不一致: R1={h1!r} R2={h2!r}"
                    )
                sig = _signature(s1, q1, s2, q2)
                inserted = conn.execute(
                    "INSERT OR IGNORE INTO seen(rid, sig) VALUES (?, ?)", (rid1, sig)
                ).rowcount
                if inserted == 1:  # 首次出现:保留并写出
                    kept += 1
                    w1.write(f"@{h1}\n{s1}\n{p1}\n{q1}\n")
                    w2.write(f"@{h2}\n{s2}\n{p2}\n{q2}\n")
                else:
                    row = conn.execute(
                        "SELECT sig FROM seen WHERE rid = ?", (rid1,)
                    ).fetchone()
                    if row is None:  # 防御:刚插入过的键不应查不到
                        raise DedupError(f"第 {total} 对读段号 {rid1!r} 索引异常")
                    if row[0] != sig:
                        raise DuplicateContentConflict(
                            f"第 {total} 对读段号 {rid1!r} 重复但内容不同"
                            "(同号不同内容),整次去重失败"
                        )
                    dropped += 1
                if total % _COMMIT_EVERY == 0:
                    conn.commit()
        conn.commit()
    finally:
        conn.close()

    if total != expected_pairs:
        raise PairCountMismatch(
            f"读段对数 {total} 与原入库对数 {expected_pairs} 不一致,"
            "源文件可能损坏或已变化"
        )
    if kept + dropped != total:
        raise DedupError(f"计数不自洽: 保留 {kept} + 剔除 {dropped} != 输入 {total}")
    return {
        "input_pairs": total,
        "kept_pairs": kept,
        "dropped_pairs": dropped,
        "r1_sha256": stream1.hasher.hexdigest(),
        "r2_sha256": stream2.hasher.hexdigest(),
    }


def stat_fingerprint(path):
    """文件身份快照:设备、inode、大小、mtime(纳秒)。"""
    st = os.stat(path)
    return (st.st_dev, st.st_ino, st.st_size, st.st_mtime_ns)


def file_sha256(path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        while chunk := fh.read(1 << 20):
            digest.update(chunk)
    return digest.hexdigest()


def verify_source_unchanged(path, fingerprint, expected_sha256):
    """发布前复核源文件未被替换/修改。

    先比 stat,再完整重读算哈希(即使大小与 mtime 被刻意保持一致也能发现内容
    变化),最后再比一次 stat,缩小发布前的竞态窗口。
    """
    if stat_fingerprint(path) != fingerprint:
        raise SourceChangedError(f"源文件在处理期间发生变化: {os.path.basename(path)}")
    if file_sha256(path) != expected_sha256:
        raise SourceChangedError(
            f"源文件在处理期间内容变化(哈希不符): {os.path.basename(path)}"
        )
    if stat_fingerprint(path) != fingerprint:
        raise SourceChangedError(f"源文件在处理期间发生变化: {os.path.basename(path)}")
