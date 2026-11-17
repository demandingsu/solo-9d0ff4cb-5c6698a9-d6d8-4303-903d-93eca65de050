"""读段号去重关键情形自测:规则、配对/计数校验、原子发布、源变化与恢复。"""
import os
import sqlite3
import sys
import time
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from app import dedup as dedup_mod  # noqa: E402
from app.dedup import (  # noqa: E402
    DedupError,
    DuplicateContentConflict,
    FastqFormatError,
    PairCountMismatch,
    PairingError,
    deduplicate,
    verify_source_unchanged,
)

SAMPLES_CSV = "sample_id,index\nS1,AAAACCCC\nS2,AAAACCGG\n"
SEQ1, SEQ2 = "ACGT" * 10, "TGCA" * 10
QUAL = "I" * 40


def rec(head, seq=SEQ1, qual=QUAL):
    return f"@{head} 1:N:0:0\n{seq}\n+\n{qual}\n"


def write_fastq(path, records):
    Path(path).write_text("".join(records), encoding="ascii")


def make_pair(r1_path, r2_path, entries):
    """entries: [(读段号, r1序列, r1质量, r2序列, r2质量)]。"""
    write_fastq(r1_path, [f"@{rid} 1:N:0:0\n{s1}\n+\n{q1}\n"
                          for rid, s1, q1, _s2, _q2 in entries])
    write_fastq(r2_path, [f"@{rid} 2:N:0:0\n{s2}\n+\n{q2}\n"
                          for rid, _s1, _q1, s2, q2 in entries])


def run_dedup(tmp_path, entries, expected=None, r1_name="r1.fastq",
              r2_name="r2.fastq"):
    r1, r2 = tmp_path / r1_name, tmp_path / r2_name
    make_pair(r1, r2, entries)
    out1, out2 = tmp_path / "out_R1.fastq", tmp_path / "out_R2.fastq"
    stat = deduplicate(
        str(r1), str(r2), len(entries) if expected is None else expected,
        str(out1), str(out2), str(tmp_path / "seen.db"),
    )
    return stat, out1, out2


def ids_of(path):
    lines = Path(path).read_text().splitlines()
    return [lines[i][1:].split()[0] for i in range(0, len(lines), 4)]


def e(rid, s1=SEQ1, q1=QUAL, s2=SEQ2, q2=QUAL):
    return (rid, s1, q1, s2, q2)


# ------------------------------------------------------------ 单元级

def test_keep_first_drop_identical_duplicates(tmp_path):
    stat, out1, out2 = run_dedup(tmp_path, [
        e("a"), e("b"), e("a"), e("c"), e("b"),
    ])
    assert stat["input_pairs"] == 5
    assert stat["kept_pairs"] == 3 and stat["dropped_pairs"] == 2
    assert ids_of(out1) == ["a", "b", "c"]  # 原顺序,保留首次出现
    assert ids_of(out2) == ["a", "b", "c"]
    # 输出是完整 FASTQ(4 行/条),内容可取回
    assert out1.read_text().count("\n") == 12
    assert SEQ1 in out1.read_text() and SEQ2 in out2.read_text()


def test_conflict_aborts_whole_run(tmp_path):
    with pytest.raises(DuplicateContentConflict, match="重复但内容不同"):
        run_dedup(tmp_path, [e("a"), e("b"), e("a", s1="AAAA", q1="IIII")])


def test_conflict_on_quality_only(tmp_path):
    with pytest.raises(DuplicateContentConflict):
        run_dedup(tmp_path, [e("a"), e("a", q1="J" * 40)])


def test_conflict_on_r2_only(tmp_path):
    with pytest.raises(DuplicateContentConflict):
        run_dedup(tmp_path, [e("a"), e("a", s2="AAAA", q2="IIII")])


def test_same_seq_different_id_kept(tmp_path):
    stat, out1, _ = run_dedup(tmp_path, [e("a"), e("b")])  # 内容相同但号不同
    assert stat["kept_pairs"] == 2 and stat["dropped_pairs"] == 0
    assert ids_of(out1) == ["a", "b"]


def test_read_id_ignores_suffix(tmp_path):
    r1 = tmp_path / "r1.fastq"
    r2 = tmp_path / "r2.fastq"
    write_fastq(r1, [rec("a"), rec("a")])
    write_fastq(r2, [rec("a"), rec("a")])
    stat = deduplicate(str(r1), str(r2), 2, str(tmp_path / "o1"), str(tmp_path / "o2"),
                       str(tmp_path / "seen.db"))
    assert stat["dropped_pairs"] == 1


def test_signature_length_prefix(tmp_path):
    """不同字段切分不产生同一指纹(长度前缀)。"""
    a = dedup_mod._signature("AC", "II", "GT", "II")
    b = dedup_mod._signature("ACG", "TII", "GT", "II")
    assert a != b


def test_pairing_mismatch(tmp_path):
    make_pair(tmp_path / "r1.fastq", tmp_path / "r2.fastq", [e("a"), e("b")])
    write_fastq(tmp_path / "r2.fastq", [rec("a")])  # R2 少一条
    with pytest.raises(PairingError, match="记录数不一致"):
        deduplicate(str(tmp_path / "r1.fastq"), str(tmp_path / "r2.fastq"), 2,
                    str(tmp_path / "o1"), str(tmp_path / "o2"),
                    str(tmp_path / "seen.db"))

    make_pair(tmp_path / "r1.fastq", tmp_path / "r2.fastq", [e("a")])
    write_fastq(tmp_path / "r2.fastq", [rec("b")])  # 号不同
    with pytest.raises(PairingError, match="读段号不一致"):
        deduplicate(str(tmp_path / "r1.fastq"), str(tmp_path / "r2.fastq"), 1,
                    str(tmp_path / "o1"), str(tmp_path / "o2"),
                    str(tmp_path / "seen.db"))


def test_expected_pairs_mismatch(tmp_path):
    with pytest.raises(PairCountMismatch, match="不一致"):
        run_dedup(tmp_path, [e("a"), e("b")], expected=3)


def test_corrupt_input(tmp_path):
    r1, r2 = tmp_path / "r1.fastq", tmp_path / "r2.fastq"
    write_fastq(r1, ["\n".join(rec("a").splitlines()[:3]) + "\n"])  # 缺质量行
    write_fastq(r2, [rec("a")])
    with pytest.raises(FastqFormatError, match="截断"):
        deduplicate(str(r1), str(r2), 1, str(tmp_path / "o1"),
                    str(tmp_path / "o2"), str(tmp_path / "seen.db"))


def test_non_text_rejected(tmp_path):
    r1, r2 = tmp_path / "r1.fastq", tmp_path / "r2.fastq"
    r1.write_bytes(b"@a\n\x1f\x8b\x08\n+\nIIII\n")
    write_fastq(r2, [rec("a")])
    with pytest.raises(FastqFormatError, match="纯文本"):
        deduplicate(str(r1), str(r2), 1, str(tmp_path / "o1"),
                    str(tmp_path / "o2"), str(tmp_path / "seen.db"))


def test_source_change_detected(tmp_path):
    r1, r2 = tmp_path / "r1.fastq", tmp_path / "r2.fastq"
    make_pair(r1, r2, [e("a")])
    fp = dedup_mod.stat_fingerprint(r1)
    stat = deduplicate(str(r1), str(r2), 1, str(tmp_path / "o1"),
                       str(tmp_path / "o2"), str(tmp_path / "seen.db"))
    verify_source_unchanged(r1, fp, stat["r1_sha256"])  # 未变:通过

    # 保持大小与 mtime 不变,只改内容 -> 哈希复核仍能发现
    before = r1.stat()
    r1.write_text(r1.read_text().replace("ACGT", "AAAA"))
    os.utime(r1, ns=(before.st_atime_ns, before.st_mtime_ns))
    with pytest.raises(DedupError, match="内容变化"):
        verify_source_unchanged(r1, fp, stat["r1_sha256"])


def test_index_on_disk_not_in_memory(tmp_path):
    """读段号索引必须落在 SQLite 文件里,而不是 Python 集合/字典。"""
    stat, _, _ = run_dedup(tmp_path, [e(f"r{i}") for i in range(50)])
    db = tmp_path / "seen.db"
    assert db.is_file()
    conn = sqlite3.connect(str(db))
    assert conn.execute("SELECT count(*) FROM seen").fetchone()[0] == 50
    conn.close()


def test_large_stream_bounded(tmp_path, monkeypatch):
    """大量读段:索引库有磁盘占用,同时流式输出正确。"""
    n = 5000
    entries = [e(f"r{i % 3000}") for i in range(n)]  # 3000 个唯一号,2000 条重复
    stat, out1, _ = run_dedup(tmp_path, entries)
    assert stat["kept_pairs"] == 3000 and stat["dropped_pairs"] == 2000
    assert len(ids_of(out1)) == 3000
    assert ids_of(out1)[:3] == ["r0", "r1", "r2"]  # 保留首次出现


# ------------------------------------------------------------ 端到端

@pytest.fixture()
def client(tmp_path, monkeypatch):
    monkeypatch.setenv("DATA_DIR", str(tmp_path / "data"))
    from fastapi.testclient import TestClient
    from app.main import app
    with TestClient(app) as c:
        c.data_dir = tmp_path / "data"
        yield c


def submit_demux(client, run_number="RUN1"):
    reads = []
    for rid, idx, seq1, seq2 in [
        ("r001", "AAAACCCC", SEQ1, SEQ2),
        ("r002", "AAAACCCC", SEQ1, SEQ2),   # 与 r001 内容相同的重复号
        ("r003", "AAAACCGG", SEQ1, SEQ2),
        ("r001", "AAAACCCC", SEQ1, SEQ2),   # 再一次重复
        ("r004", "AAAACCCC", SEQ1, SEQ2),   # S1 中的唯一条
    ]:
        reads.append((rid, idx, seq1, seq2))
    r1 = "".join(f"@{rid} 1:N:0:0\n{s1}\n+\n{QUAL}\n" for rid, _i, s1, _s2 in reads)
    r2 = "".join(f"@{rid} 2:N:0:0\n{s2}\n+\n{QUAL}\n" for rid, _i, _s1, s2 in reads)
    i1 = "".join(f"@{rid} 1:N:0:0\n{idx}\n+\n{'I' * len(idx)}\n"
                 for rid, idx, _s1, _s2 in reads)
    return client.post(
        "/runs", data={"run_number": run_number},
        files={
            "r1": ("r1.fastq", r1.encode()),
            "r2": ("r2.fastq", r2.encode()),
            "i1": ("i1.fastq", i1.encode()),
            "samples": ("samples.csv", SAMPLES_CSV.encode()),
        },
        follow_redirects=False,
    )


def wait_done(client, url, timeout=15):
    deadline = time.time() + timeout
    while time.time() < deadline:
        st = client.get(url).json()
        if st["status"] != "processing":
            return st
        time.sleep(0.05)
    raise TimeoutError("处理超时")


def publish_demux(client, run_number="RUN1"):
    resp = submit_demux(client, run_number)
    assert resp.status_code == 303, resp.text
    rid = int(resp.headers["location"].rsplit("/", 1)[1])
    assert wait_done(client, f"/api/runs/{rid}")["status"] == "published"
    return rid


def start_dedup(client, run_id, sample_id="S1"):
    resp = client.post(f"/runs/{run_id}/dedup", data={"sample_id": sample_id},
                       follow_redirects=False)
    assert resp.status_code == 303, resp.text
    return int(resp.headers["location"].rsplit("/", 1)[1])


def test_end_to_end_dedup(client):
    rid = publish_demux(client)
    # 原文件内容作为对照
    before = client.get(f"/runs/{rid}/download/S1_R1.fastq").text
    jid = start_dedup(client, rid, "S1")
    st = wait_done(client, f"/api/dedup/{jid}")
    assert st["status"] == "published", st["error"]
    assert st["expected_pairs"] == 4 and st["input_pairs"] == 4
    assert st["kept_pairs"] == 3 and st["dropped_pairs"] == 1

    o1 = client.get(f"/dedup/{jid}/download/S1_dedup_R1.fastq")
    o2 = client.get(f"/dedup/{jid}/download/S1_dedup_R2.fastq")
    assert o1.status_code == 200 and o2.status_code == 200
    lines1 = o1.text.splitlines()
    assert [lines1[i][1:].split()[0] for i in range(0, len(lines1), 4)] == [
        "r001", "r002", "r004"]
    assert SEQ2 in o2.text  # R2 内容正确

    log = client.get(f"/dedup/{jid}/download/dedup.log")
    assert "保留对数: 3" in log.text and "剔除对数: 1" in log.text
    assert "原入库对数: 4" in log.text

    # 原文件不变
    assert client.get(f"/runs/{rid}/download/S1_R1.fastq").text == before

    page = client.get(f"/dedup/{jid}")
    assert page.status_code == 200 and "保留对数" in page.text
    run_page = client.get(f"/runs/{rid}")
    assert "查看去重结果" in run_page.text


def test_dedup_conflict_fails_without_partial_publish(client):
    rid = publish_demux(client, "RUN2")
    # 直接改写已发布源文件,构造"同号不同内容"(模拟源数据本身冲突)
    pub = client.data_dir / "published" / "RUN2" / "S1_R1.fastq"
    pub2 = client.data_dir / "published" / "RUN2" / "S1_R2.fastq"
    l1, l2 = pub.read_text().splitlines(), pub2.read_text().splitlines()
    first = l1[0][1:].split()[0]
    assert l1[4][1:].split()[0] != first
    # 把第 2 条改成与第 1 条同号(内容不同),总对数不变
    l1[4] = f"@{first} 1:N:0:0"
    l1[5] = "A" * len(l1[5])  # 同号但 R1 序列不同
    l2[4] = f"@{first} 2:N:0:0"
    pub.write_text("\n".join(l1) + "\n")
    pub2.write_text("\n".join(l2) + "\n")

    jid = start_dedup(client, rid, "S1")
    st = wait_done(client, f"/api/dedup/{jid}")
    assert st["status"] == "failed" and "重复但内容不同" in st["error"]
    assert st["kept_pairs"] is None and st["dropped_pairs"] is None  # 不发布计数
    assert st["input_pairs"] is None
    assert not (client.data_dir / "published" / "_dedup" / "RUN2").exists()
    assert client.get(f"/dedup/{jid}/download/S1_dedup_R1.fastq").status_code == 409


def test_dedup_source_missing_fails(client):
    rid = publish_demux(client, "RUN3")
    (client.data_dir / "published" / "RUN3" / "S1_R2.fastq").unlink()
    jid = start_dedup(client, rid, "S1")
    st = wait_done(client, f"/api/dedup/{jid}")
    assert st["status"] == "failed" and "缺失" in st["error"]
    assert not (client.data_dir / "published" / "_dedup" / "RUN3").exists()


def test_dedup_source_changed_during_processing(client, monkeypatch):
    rid = publish_demux(client, "RUN4")
    real = dedup_mod.deduplicate

    def tamper(*args, **kwargs):
        stat = real(*args, **kwargs)
        src = Path(args[0])
        src.write_text(src.read_text())  # 改写(内容相同但 mtime/inode 可能变)
        return stat

    monkeypatch.setattr("app.main.deduplicate", tamper)
    jid = start_dedup(client, rid, "S1")
    st = wait_done(client, f"/api/dedup/{jid}")
    assert st["status"] == "failed" and "变化" in st["error"]
    assert not (client.data_dir / "published" / "_dedup" / "RUN4").exists()


def test_dedup_requires_published_run(client):
    resp = client.post("/runs/999/dedup", data={"sample_id": "S1"})
    assert resp.status_code == 404
    rid = publish_demux(client, "RUN5")
    assert client.post(f"/runs/{rid}/dedup", data={"sample_id": "NOPE"}).status_code == 404
    assert client.post(f"/runs/{rid}/dedup", data={"sample_id": "undetermined"}).status_code == 404


def test_dedup_published_immutable_and_retry(client):
    rid = publish_demux(client, "RUN6")
    jid = start_dedup(client, rid, "S1")
    assert wait_done(client, f"/api/dedup/{jid}")["status"] == "published"
    before = client.get(f"/dedup/{jid}/download/S1_dedup_R1.fastq").text
    resp = client.post(f"/runs/{rid}/dedup", data={"sample_id": "S1"},
                       follow_redirects=False)
    assert resp.status_code == 409 and "不可改写" in resp.json()["detail"]
    after = client.get(f"/dedup/{jid}/download/S1_dedup_R1.fastq").text
    assert before == after



def test_dedup_failed_job_retry(client):
    rid = publish_demux(client, "RUN7")
    src = client.data_dir / "published" / "RUN7" / "S1_R1.fastq"
    original = src.read_text()
    src.unlink()
    jid = start_dedup(client, rid, "S1")
    assert wait_done(client, f"/api/dedup/{jid}")["status"] == "failed"
    src.write_text(original)  # 恢复源文件
    jid2 = start_dedup(client, rid, "S1")
    st = wait_done(client, f"/api/dedup/{jid2}")
    assert st["status"] == "published", st["error"]
    assert client.get(f"/dedup/{jid2}/download/S1_dedup_R1.fastq").status_code == 200


def test_dedup_recovered_on_startup(client):
    from app.main import recover_interrupted
    from app.db import db_conn
    from app.config import db_path, staging_dir, published_dir

    rid = publish_demux(client, "RUN8")
    staging = staging_dir() / "deadcafe"
    (staging / "out").mkdir(parents=True)
    (staging / "out" / "S1_dedup_R1.fastq").write_text("半成品")
    stray = published_dir() / "dedup" / "RUN8" / "S1"
    stray.mkdir(parents=True)
    (stray / "S1_dedup_R1.fastq").write_text("改名后未落库")
    with db_conn(db_path()) as conn:
        conn.execute(
            "INSERT INTO dedup_jobs(run_id, sample_id, status, staging_dir, source_r1, "
            "source_r2, source_r1_name, source_r2_name, expected_pairs, created_at, "
            "updated_at) VALUES (?,?,?,?,?,?,?,?,?,?,?)",
            (rid, "S1", "processing", str(staging), "x", "y", "x", "y", 4, "t", "t"),
        )

    recover_interrupted()

    with db_conn(db_path()) as conn:
        row = conn.execute("SELECT * FROM dedup_jobs WHERE run_id=?", (rid,)).fetchone()
    assert row["status"] == "failed" and "中断" in row["error"]
    assert not staging.exists() and not stray.exists()
    jid = start_dedup(client, rid, "S1")
    assert wait_done(client, f"/api/dedup/{jid}")["status"] == "published"


def test_demux_pages_still_work(client):
    resp = client.get("/")
    assert resp.status_code == 200 and "新建拆样运行" in resp.text
    rid = publish_demux(client, "RUN9")
    page = client.get(f"/runs/{rid}")
    assert page.status_code == 200 and "已发布" in page.text


def test_dedup_source_run_invalidated_during_processing(client, monkeypatch):
    rid = publish_demux(client, "RUN10")
    real = dedup_mod.deduplicate
    from app.db import db_conn
    from app.config import db_path

    def invalidate(*args, **kwargs):
        stat = real(*args, **kwargs)
        with db_conn(db_path()) as conn:  # 处理期间源运行失效
            conn.execute("UPDATE runs SET status='failed' WHERE id=?", (rid,))
        return stat

    monkeypatch.setattr("app.main.deduplicate", invalidate)
    jid = start_dedup(client, rid, "S1")
    st = wait_done(client, f"/api/dedup/{jid}")
    assert st["status"] == "failed" and "已变化" in st["error"]
    assert st["kept_pairs"] is None
    assert not (client.data_dir / "published" / "_dedup" / "RUN10").exists()


def test_dedup_failed_page(client):
    rid = publish_demux(client, "RUN11")
    (client.data_dir / "published" / "RUN11" / "S1_R1.fastq").unlink()
    jid = start_dedup(client, rid, "S1")
    assert wait_done(client, f"/api/dedup/{jid}")["status"] == "failed"
    page = client.get(f"/dedup/{jid}")
    assert page.status_code == 200
    assert "未发布任何文件" in page.text and "失败原因" in page.text
    assert "保留对数" not in page.text
    home = client.get("/")
    assert home.status_code == 200 and "读段号去重记录" in home.text
    run_page = client.get(f"/runs/{rid}")
    assert "上次去重失败" in run_page.text and "重新去重" in run_page.text
