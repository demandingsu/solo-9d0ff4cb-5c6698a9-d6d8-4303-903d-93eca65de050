"""泳道拆样台:按样本索引拆分 R1/R2/I1,并可对已发布样本做读段号去重。

去重发布模型与拆样一致:
- 处理中:去重索引与半成品输出都在 staging/<token>/,不对外可见。
- 成功:输出目录整体 os.rename 到 published/_dedup/<运行号>/<样本>/,同事务
  写入保留/剔除计数并把任务置为 published;之后才允许下载。
- 失败:任务标 failed,暂存与半成品目录清理;不留下部分文件或计数。
- 源运行失效或源文件变化:同样整次失败,不发布。
- 再次去重同一运行同一 Read1/Read2:返回 409,已发布结果保持原内容。
"""
import os
import re
import shutil
import sqlite3
import threading
import uuid
from contextlib import asynccontextmanager
from datetime import datetime, timezone
from pathlib import Path

from fastapi import FastAPI, File, Form, HTTPException, Request, UploadFile
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse, RedirectResponse
from fastapi.templating import Jinja2Templates

from .config import db_path, published_dir, staging_dir
from .db import db_conn, init_db
from .dedup import (
    DedupError,
    deduplicate,
    stat_fingerprint,
    verify_source_unchanged,
)
from .demux import demultiplex
from .samples import load_samples

RUN_NUMBER_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,63}$")
STATUS_LABELS = {"processing": "处理中", "published": "已发布", "failed": "失败"}
UNDETERMINED = "undetermined"
# 去重发布区放在发布区根下的保留子目录里。运行号必须以字母/数字开头,
# 因此 "_dedup" 不可能与任何运行的发布目录重名。
DEDUP_DIRNAME = "_dedup"

BASE_DIR = Path(__file__).resolve().parent
templates = Jinja2Templates(directory=str(BASE_DIR / "templates"))


def utcnow() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _remove(path) -> None:
    p = Path(path)
    if p.is_dir() and not p.is_symlink():
        shutil.rmtree(p, ignore_errors=True)
    else:
        try:
            p.unlink()
        except (FileNotFoundError, NotADirectoryError):
            pass


def dedup_published_dir(run_number: str, sample_id: str) -> Path:
    return published_dir() / DEDUP_DIRNAME / run_number / sample_id


def recover_interrupted() -> None:
    """启动恢复:凡未走到 published 的运行/去重任务一律标失败并清掉半成品,
    保证发布区里只存在完整运行的完整文件。"""
    with db_conn(db_path()) as conn:
        stale = conn.execute(
            "SELECT id, run_number, staging_dir FROM runs WHERE status = 'processing'"
        ).fetchall()
        for row in stale:
            conn.execute(
                "UPDATE runs SET status='failed', error=?, updated_at=? WHERE id=?",
                ("进程中断,未完成即停止;可用同一运行号重新提交", utcnow(), row["id"]),
            )
        published = {
            r["run_number"]
            for r in conn.execute("SELECT run_number FROM runs WHERE status='published'")
        }
        stale_jobs = conn.execute(
            "SELECT j.id, j.staging_dir, r.run_number, j.sample_id "
            "FROM dedup_jobs j JOIN runs r ON r.id = j.run_id "
            "WHERE j.status = 'processing'"
        ).fetchall()
        for row in stale_jobs:
            conn.execute(
                "UPDATE dedup_jobs SET status='failed', error=?, updated_at=? WHERE id=?",
                ("进程中断,未完成即停止;可再次发起同一去重", utcnow(), row["id"]),
            )
        published_jobs = {
            (r["run_number"], r["sample_id"])
            for r in conn.execute(
                "SELECT r.run_number, j.sample_id FROM dedup_jobs j "
                "JOIN runs r ON r.id = j.run_id WHERE j.status='published'"
            )
        }
    for row in stale:
        _remove(row["staging_dir"])
        _remove(published_dir() / row["run_number"])
    for row in stale_jobs:
        _remove(row["staging_dir"])
        _remove(dedup_published_dir(row["run_number"], row["sample_id"]))
    if staging_dir().exists():
        for child in staging_dir().iterdir():  # 无主的暂存目录一并清理
            _remove(child)
    if published_dir().exists():
        for child in published_dir().iterdir():
            if child.name == DEDUP_DIRNAME:
                continue
            if child.name not in published:
                _remove(child)
    dedup_root = published_dir() / DEDUP_DIRNAME
    if dedup_root.exists():
        known_runs = published | {r for r, _s in published_jobs}
        for run_dir in dedup_root.iterdir():
            if not run_dir.is_dir():
                continue
            if run_dir.name not in known_runs:
                _remove(run_dir)  # 运行已不存在/已失效:整目录都不该提供下载
                continue
            for sample_dir in run_dir.iterdir():  # 无主/未落库的去重目录
                if (run_dir.name, sample_dir.name) not in published_jobs:
                    _remove(sample_dir)
            if not any(run_dir.iterdir()):
                _remove(run_dir)
        if not any(dedup_root.iterdir()):
            _remove(dedup_root)  # 空的保留目录本身也清掉
    with db_conn(db_path()) as conn:  # 已发布但文件丢失:不允许提供残缺下载
        for num in published:
            if not (published_dir() / num).is_dir():
                conn.execute(
                    "UPDATE runs SET status='failed', error=?, updated_at=? "
                    "WHERE run_number=? AND status='published'",
                    ("已发布文件缺失,运行标记为失败", utcnow(), num),
                )
        for run_number, sample_id in published_jobs:
            needed = [f"{sample_id}_dedup_R1.fastq", f"{sample_id}_dedup_R2.fastq",
                      "dedup.log"]
            d = dedup_published_dir(run_number, sample_id)
            if not d.is_dir() or not all((d / n).is_file() for n in needed):
                conn.execute(
                    "UPDATE dedup_jobs SET status='failed', error=?, updated_at=? "
                    "WHERE sample_id=? AND status='published' AND run_id IN "
                    "(SELECT id FROM runs WHERE run_number=?)",
                    ("已发布去重文件缺失,任务标记为失败", utcnow(), sample_id,
                     run_number),
                )


def setup() -> None:
    staging_dir().mkdir(parents=True, exist_ok=True)
    published_dir().mkdir(parents=True, exist_ok=True)
    init_db(db_path())
    recover_interrupted()


@asynccontextmanager
async def lifespan(_app: FastAPI):
    setup()
    yield


app = FastAPI(title="泳道拆样台", lifespan=lifespan)


# ---------------------------------------------------------------- 去重后台处理

def write_dedup_log(path: Path, job, run, stat) -> None:
    lines = [
        f"运行号: {run['run_number']}",
        f"样本: {job['sample_id']}",
        f"源文件: R1={job['source_r1_name']} R2={job['source_r2_name']}",
        f"发起时间: {job['created_at']}",
        f"完成时间: {utcnow()}",
        "规则: 按原顺序保留首次出现的读段对;同号且 R1/R2 序列与质量均相同的对剔除;",
        "      同号内容不同则整次失败。",
        f"原入库对数: {job['expected_pairs']}",
        f"输入对数: {stat['input_pairs']}",
        f"保留对数: {stat['kept_pairs']}",
        f"剔除对数: {stat['dropped_pairs']}",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def process_dedup(job_id: int) -> None:
    with db_conn(db_path()) as conn:
        job = conn.execute("SELECT * FROM dedup_jobs WHERE id=?", (job_id,)).fetchone()
        run = None if job is None else conn.execute(
            "SELECT * FROM runs WHERE id=?", (job["run_id"],)
        ).fetchone()
    if job is None:
        return
    staging = Path(job["staging_dir"])
    out = staging / "out"
    final = None
    try:
        if run is None or run["status"] != "published":
            raise DedupError("源运行已失效(不是已发布状态),不能去重")
        source_r1, source_r2 = Path(job["source_r1"]), Path(job["source_r2"])
        if not source_r1.is_file() or not source_r2.is_file():
            raise DedupError("源运行文件缺失,不能去重")
        # 处理前锁定源文件身份,处理结束后复核,防止期间被替换或改写。
        fp1, fp2 = stat_fingerprint(source_r1), stat_fingerprint(source_r2)
        out.mkdir(parents=True)
        _remove(staging / "seen.db")  # 失败重试不得复用上次的读段号索引
        stat = deduplicate(
            str(source_r1), str(source_r2), job["expected_pairs"],
            str(out / f"{job['sample_id']}_dedup_R1.fastq"),
            str(out / f"{job['sample_id']}_dedup_R2.fastq"),
            str(staging / "seen.db"),
        )
        verify_source_unchanged(source_r1, fp1, stat["r1_sha256"])
        verify_source_unchanged(source_r2, fp2, stat["r2_sha256"])
        final = dedup_published_dir(run["run_number"], job["sample_id"])
        write_dedup_log(out / "dedup.log", job, run, stat)
        if final.exists():
            raise DedupError(f"发布目录已存在: {final}")
        final.parent.mkdir(parents=True, exist_ok=True)
        os.rename(out, final)  # 同一数据卷内原子移动:先完整落盘,再发布
        # 计数只在此刻(文件已就位)与 published 状态同一事务落库:
        # 之前任何失败都不会让人看到"有计数却没有完整文件"。
        with db_conn(db_path()) as conn:
            updated = conn.execute(
                "UPDATE dedup_jobs SET status='published', input_pairs=?, kept_pairs=?, "
                "dropped_pairs=?, published_at=?, updated_at=? "
                "WHERE id=? AND status='processing' AND EXISTS ("
                "SELECT 1 FROM runs WHERE id=? AND status='published')",
                (stat["input_pairs"], stat["kept_pairs"], stat["dropped_pairs"],
                 utcnow(), utcnow(), job_id, job["run_id"]),
            ).rowcount
            if updated != 1:
                raise DedupError("任务状态或源运行状态已变化,放弃发布")
        _remove(staging)
    except Exception as exc:  # 任何失败:标失败 + 清半成品,不发布文件或计数
        with db_conn(db_path()) as conn:
            conn.execute(
                "UPDATE dedup_jobs SET status='failed', error=?, updated_at=? "
                "WHERE id=? AND status='processing'",
                (str(exc), utcnow(), job_id),
            )
        _remove(staging)
        if final is not None:
            with db_conn(db_path()) as conn:
                still_published = conn.execute(
                    "SELECT 1 FROM dedup_jobs WHERE id=? AND status='published'", (job_id,)
                ).fetchone()
            if not still_published:  # 已发布结果保持原内容,不被清理
                _remove(final)
                parent = final.parent
                if parent.is_dir() and not any(parent.iterdir()):
                    _remove(parent)


# ---------------------------------------------------------------- 页面与接口

@app.exception_handler(HTTPException)
async def http_exception_handler(request: Request, exc: HTTPException):
    if "text/html" in request.headers.get("accept", ""):
        return templates.TemplateResponse(
            request, "error.html",
            {"status": exc.status_code, "detail": exc.detail},
            status_code=exc.status_code,
        )
    return JSONResponse({"detail": exc.detail}, status_code=exc.status_code)


@app.get("/", response_class=HTMLResponse)
def index(request: Request):
    with db_conn(db_path()) as conn:
        runs = conn.execute("SELECT * FROM runs ORDER BY id DESC").fetchall()
        jobs = conn.execute(
            "SELECT j.*, r.run_number FROM dedup_jobs j JOIN runs r ON r.id = j.run_id "
            "ORDER BY j.id DESC"
        ).fetchall()
    return templates.TemplateResponse(
        request, "index.html",
        {"runs": runs, "jobs": jobs, "labels": STATUS_LABELS},
    )


async def save_upload(upload: UploadFile, dest: Path) -> None:
    with open(dest, "wb") as fh:
        while chunk := await upload.read(1 << 20):
            fh.write(chunk)


@app.post("/runs")
async def create_run(
    run_number: str = Form(...),
    r1: UploadFile = File(...),
    r2: UploadFile = File(...),
    i1: UploadFile = File(...),
    samples: UploadFile = File(...),
):
    run_number = run_number.strip()
    if not RUN_NUMBER_RE.fullmatch(run_number):
        raise HTTPException(400, "运行号格式不正确:字母或数字开头,可含 - _ .,最长 64 字符")
    for label, up in (("R1", r1), ("R2", r2), ("I1", i1), ("样本表", samples)):
        if not up.filename:
            raise HTTPException(400, f"请上传{label}文件")

    token = uuid.uuid4().hex
    staging = staging_dir() / token
    uploads = staging / "uploads"
    names = {k: Path(u.filename or "").name for k, u in
             (("r1", r1), ("r2", r2), ("i1", i1), ("samples", samples))}
    try:
        with db_conn(db_path()) as conn:
            row = conn.execute(
                "SELECT id, status FROM runs WHERE run_number=?", (run_number,)
            ).fetchone()
            if row and row["status"] == "published":
                raise HTTPException(409, f"运行 {run_number} 已发布,不可改写")
            if row and row["status"] == "processing":
                raise HTTPException(409, f"运行 {run_number} 正在处理中,请稍后查看")
            if row:  # failed:允许同号重试
                conn.execute("DELETE FROM runs WHERE id=?", (row["id"],))
            cur = conn.execute(
                "INSERT INTO runs(run_number, status, staging_dir, r1_name, r2_name, "
                "i1_name, samples_name, created_at, updated_at) "
                "VALUES (?,?,?,?,?,?,?,?,?)",
                (run_number, "processing", str(staging), names["r1"], names["r2"],
                 names["i1"], names["samples"], utcnow(), utcnow()),
            )
            run_id = cur.lastrowid
    except sqlite3.IntegrityError:
        raise HTTPException(409, f"运行 {run_number} 已存在,请刷新后查看")

    try:
        uploads.mkdir(parents=True)
        await save_upload(r1, uploads / "r1.fastq")
        await save_upload(r2, uploads / "r2.fastq")
        await save_upload(i1, uploads / "i1.fastq")
        await save_upload(samples, uploads / "samples.csv")
        load_samples(uploads / "samples.csv")  # 先校验样本表,尽早报错
    except HTTPException:
        raise
    except Exception as exc:
        with db_conn(db_path()) as conn:
            conn.execute(
                "UPDATE runs SET status='failed', error=?, updated_at=? WHERE id=?",
                (f"上传或样本表错误: {exc}", utcnow(), run_id),
            )
        _remove(staging)
        raise HTTPException(400, f"上传或样本表错误: {exc}")

    threading.Thread(target=process_run, args=(run_id,), daemon=True).start()
    return RedirectResponse(f"/runs/{run_id}", status_code=303)


def write_log(path: Path, run, sample_list, stats) -> None:
    lines = [
        f"运行号: {run['run_number']}",
        f"提交时间: {run['created_at']}",
        f"完成时间: {utcnow()}",
        f"输入文件: R1={run['r1_name']} R2={run['r2_name']} "
        f"I1={run['i1_name']} 样本表={run['samples_name']}",
        "样本索引:",
    ]
    lines += [f"  {sid}\t{idx}" for sid, idx in sample_list]
    lines.append(f"总读段对: {stats['total_pairs']}")
    lines.append("分组对数(对数\t其中单碱基纠错):")
    for sid in stats["order"]:
        pairs, corrected = stats["counts"][sid]
        lines.append(f"  {sid}\t{pairs}\t{corrected}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def process_run(run_id: int) -> None:
    with db_conn(db_path()) as conn:
        run = conn.execute("SELECT * FROM runs WHERE id=?", (run_id,)).fetchone()
    if run is None:
        return
    run_number = run["run_number"]
    staging = Path(run["staging_dir"])
    uploads = staging / "uploads"
    out = staging / "out"
    final = published_dir() / run_number
    try:
        sample_list = load_samples(uploads / "samples.csv")
        out.mkdir(parents=True)
        stats = demultiplex(
            str(uploads / "r1.fastq"), str(uploads / "r2.fastq"),
            str(uploads / "i1.fastq"), sample_list, str(out),
        )
        write_log(out / "demux.log", run, sample_list, stats)
        if final.exists():
            raise RuntimeError(f"发布目录已存在: {final}")
        os.rename(out, final)  # 同一数据卷内原子移动:先完整落盘,再发布
        with db_conn(db_path()) as conn:
            updated = conn.execute(
                "UPDATE runs SET status='published', total_pairs=?, published_at=?, "
                "updated_at=? WHERE id=? AND status='processing'",
                (stats["total_pairs"], utcnow(), utcnow(), run_id),
            ).rowcount
            if updated != 1:
                raise RuntimeError("运行状态已变化,放弃发布")
            index_of = dict(sample_list)
            for sid in stats["order"]:
                pairs, corrected = stats["counts"][sid]
                conn.execute(
                    "INSERT INTO sample_counts(run_id, sample_id, index_seq, pairs, "
                    "corrected) VALUES (?,?,?,?,?)",
                    (run_id, sid, index_of.get(sid, ""), pairs, corrected),
                )
        _remove(staging)
    except Exception as exc:  # 任何失败:标失败 + 清半成品,不发布
        with db_conn(db_path()) as conn:
            conn.execute(
                "UPDATE runs SET status='failed', error=?, updated_at=? WHERE id=?",
                (str(exc), utcnow(), run_id),
            )
        _remove(staging)
        _remove(final)


def get_run_or_404(run_id: int):
    with db_conn(db_path()) as conn:
        run = conn.execute("SELECT * FROM runs WHERE id=?", (run_id,)).fetchone()
    if run is None:
        raise HTTPException(404, "运行不存在")
    return run


def get_counts(run_id: int):
    with db_conn(db_path()) as conn:
        return conn.execute(
            "SELECT * FROM sample_counts WHERE run_id=? "
            "ORDER BY sample_id = 'undetermined', sample_id",
            (run_id,),
        ).fetchall()


def get_published_jobs(run_id: int):
    with db_conn(db_path()) as conn:
        return conn.execute(
            "SELECT * FROM dedup_jobs WHERE run_id=? ORDER BY id DESC", (run_id,)
        ).fetchall()


def get_job_or_404(job_id: int):
    with db_conn(db_path()) as conn:
        job = conn.execute("SELECT * FROM dedup_jobs WHERE id=?", (job_id,)).fetchone()
    if job is None:
        raise HTTPException(404, "去重任务不存在")
    return job


@app.get("/runs/{run_id}", response_class=HTMLResponse)
def run_detail(request: Request, run_id: int):
    run = get_run_or_404(run_id)
    groups = [
        {
            "sample_id": c["sample_id"],
            "index_seq": c["index_seq"],
            "pairs": c["pairs"],
            "corrected": c["corrected"],
            "r1": f"{c['sample_id']}_R1.fastq",
            "r2": f"{c['sample_id']}_R2.fastq",
        }
        for c in get_counts(run_id)
    ]
    jobs = get_published_jobs(run_id)
    dedup = {}
    for j in jobs:  # jobs 按 id 倒序,每个样本取最新一条
        dedup.setdefault(j["sample_id"], j)
    return templates.TemplateResponse(
        request, "run.html",
        {"run": run, "groups": groups, "jobs": jobs, "dedup": dedup,
         "labels": STATUS_LABELS},
    )


@app.post("/runs/{run_id}/dedup")
def create_dedup(run_id: int, sample_id: str = Form(...)):
    run = get_run_or_404(run_id)
    sample_id = sample_id.strip()
    if "/" in sample_id or "\\" in sample_id or "\x00" in sample_id or not sample_id:
        raise HTTPException(400, "样本号无效")
    if run["status"] != "published":
        raise HTTPException(409, "只有已发布的运行才能去重")
    accepted = {
        c["sample_id"] for c in get_counts(run_id)
        if c["sample_id"] != UNDETERMINED
    }
    if sample_id not in accepted:
        raise HTTPException(404, "该运行里没有这个样本")

    token = uuid.uuid4().hex
    staging = staging_dir() / token
    run_dir = published_dir() / run["run_number"]
    try:
        with db_conn(db_path()) as conn:
            cur = conn.execute(
                "SELECT id, status FROM dedup_jobs WHERE run_id=? AND sample_id=?",
                (run_id, sample_id),
            ).fetchone()
            if cur and cur["status"] == "published":
                raise HTTPException(409, "该样本已去重发布,结果不可改写")
            if cur and cur["status"] == "processing":
                raise HTTPException(409, "该样本正在去重,请稍后查看")
            if cur:  # failed:允许再次发起
                conn.execute("DELETE FROM dedup_jobs WHERE id=?", (cur["id"],))
            job_id = conn.execute(
                "INSERT INTO dedup_jobs(run_id, sample_id, status, staging_dir, "
                "source_r1, source_r2, source_r1_name, source_r2_name, "
                "expected_pairs, created_at, updated_at) "
                "VALUES (?,?,?,?,?,?,?,?,?,?,?)",
                (run_id, sample_id, "processing", str(staging),
                 str(run_dir / f"{sample_id}_R1.fastq"),
                 str(run_dir / f"{sample_id}_R2.fastq"),
                 f"{sample_id}_R1.fastq", f"{sample_id}_R2.fastq",
                 next(c["pairs"] for c in get_counts(run_id)
                      if c["sample_id"] == sample_id),
                 utcnow(), utcnow()),
            ).lastrowid
    except sqlite3.IntegrityError:
        raise HTTPException(409, "该样本正在去重,请刷新后查看")

    threading.Thread(target=process_dedup, args=(job_id,), daemon=True).start()
    return RedirectResponse(f"/dedup/{job_id}", status_code=303)


@app.get("/dedup/{job_id}", response_class=HTMLResponse)
def dedup_detail(request: Request, job_id: int):
    job = get_job_or_404(job_id)
    run = get_run_or_404(job["run_id"])
    return templates.TemplateResponse(
        request, "dedup.html",
        {"job": job, "run": run, "labels": STATUS_LABELS},
    )


@app.get("/api/dedup/{job_id}")
def dedup_status(job_id: int):
    job = get_job_or_404(job_id)
    published = job["status"] == "published"
    return {
        "run_number": get_run_or_404(job["run_id"])["run_number"],
        "sample_id": job["sample_id"],
        "status": job["status"],
        "error": job["error"],
        "expected_pairs": job["expected_pairs"],
        # 只有已发布的完整结果才对外给出计数
        "input_pairs": job["input_pairs"] if published else None,
        "kept_pairs": job["kept_pairs"] if published else None,
        "dropped_pairs": job["dropped_pairs"] if published else None,
    }


@app.get("/dedup/{job_id}/download/{filename}")
def dedup_download(job_id: int, filename: str):
    job = get_job_or_404(job_id)
    if job["status"] != "published":
        raise HTTPException(409, "去重任务未发布,没有完整文件可下载")
    run = get_run_or_404(job["run_id"])
    valid = {
        f"{job['sample_id']}_dedup_R1.fastq",
        f"{job['sample_id']}_dedup_R2.fastq",
        "dedup.log",
    }
    if filename not in valid:
        raise HTTPException(404, "文件不存在")
    path = dedup_published_dir(run["run_number"], job["sample_id"]) / filename
    if not path.is_file():
        raise HTTPException(404, "文件不存在")
    return FileResponse(path, filename=filename, media_type="application/octet-stream")


@app.get("/api/runs/{run_id}")
def run_status(run_id: int):
    run = get_run_or_404(run_id)
    counts = get_counts(run_id)
    return {
        "run_number": run["run_number"],
        "status": run["status"],
        "error": run["error"],
        "total_pairs": run["total_pairs"],
        "counts": [
            {"sample_id": c["sample_id"], "index_seq": c["index_seq"],
             "pairs": c["pairs"], "corrected": c["corrected"]}
            for c in counts
        ],
    }


@app.get("/runs/{run_id}/download/{filename}")
def download(run_id: int, filename: str):
    run = get_run_or_404(run_id)
    if run["status"] != "published":
        raise HTTPException(409, "运行未发布,没有完整文件可下载")
    valid = {"demux.log"}
    for c in get_counts(run_id):
        valid.add(f"{c['sample_id']}_R1.fastq")
        valid.add(f"{c['sample_id']}_R2.fastq")
    if filename not in valid:
        raise HTTPException(404, "文件不存在")
    path = published_dir() / run["run_number"] / filename
    if not path.is_file():
        raise HTTPException(404, "文件不存在")
    return FileResponse(path, filename=filename, media_type="application/octet-stream")
