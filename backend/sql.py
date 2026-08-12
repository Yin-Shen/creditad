import os
import json
import datetime
from sqlalchemy import create_engine, Column, String, DateTime, Text, Integer, Float, text
from sqlalchemy.orm import declarative_base, sessionmaker

Base = declarative_base()
DB_PATH = os.getenv("TAD_DB_PATH", "sqlite:///extra_mode/tad_tokens.db")
engine = create_engine(DB_PATH, echo=False, future=True)
SessionLocal = sessionmaker(bind=engine, autoflush=False, autocommit=False)

class TokenRecord(Base):
    __tablename__ = "token_records"

    token      = Column(String(36), primary_key=True, index=True)
    created_at = Column(DateTime, default=datetime.datetime.utcnow)
    mcool_path = Column(Text, nullable=False)

    resolution   = Column(Integer,    nullable=True)
    is_threshold = Column(Float,      nullable=True)
    overlap      = Column(Float,      nullable=True)
    hash_tag     = Column(String(64), nullable=True)
    fast_mode    = Column(String(8),  nullable=True)

    status      = Column(String(16), default="pending")
    finished_at = Column(DateTime,   nullable=True)
    error_msg   = Column(Text,       nullable=True)

    # Legacy columns (kept for read-back of pre-2026-05-14 records).
    # Under method_set=v2 these are reused: di=insulation, sdi=topdom_like,
    # is=contact_contrast (path strings only — semantics carried by
    # ``output_paths_json``).
    output_di        = Column(Text, nullable=True)
    output_sdi       = Column(Text, nullable=True)
    output_is        = Column(Text, nullable=True)
    output_tad_final = Column(Text, nullable=True)
    # Full per-caller output dict (JSON-serialized).  Source of truth post-v2.
    output_paths_json = Column(Text, nullable=True)
    prediction_output = Column(Text, nullable=True)  # 预测结果路径


def _ensure_output_paths_json_column() -> None:
    """Idempotent migration: add ``output_paths_json`` column on existing DBs."""
    try:
        with engine.begin() as conn:
            rows = conn.execute(text("PRAGMA table_info(token_records)")).fetchall()
            cols = {row[1] for row in rows}
            if "output_paths_json" not in cols:
                conn.execute(text("ALTER TABLE token_records ADD COLUMN output_paths_json TEXT"))
    except Exception as exc:
        # Don't crash backend startup on migration glitches — log only.
        print(f"[sql] output_paths_json migration skipped: {exc}")


def init_db():
    Base.metadata.create_all(bind=engine)
    _ensure_output_paths_json_column()

def add_token_record(
    token: str,
    mcool_path: str,
    resolution: int       = None,
    is_threshold: float   = None,
    overlap: float        = None,
    hash_tag: str         = None,
    fast_mode: bool       = None,
    prediction_output: str = None,  # 新增这个参数
):
    """注册 token，写入 mcool 路径和 TAD 检测参数快照，状态置为 pending。

    幂等：同一 token 重复注册（manager 对同一文件路径返回相同 token）时用 merge 覆盖
    更新参数快照，而不是抛 PRIMARY KEY 冲突——避免旧实现里 IntegrityError 被上层静默吞掉、
    且 DB 残留首次注册的旧 resolution。
    """
    session = SessionLocal()
    try:
        record = TokenRecord(
            token             = token,
            mcool_path        = mcool_path,
            resolution        = resolution,
            is_threshold      = is_threshold,
            overlap           = overlap,
            hash_tag          = hash_tag,
            fast_mode         = str(fast_mode) if fast_mode is not None else None,
            prediction_output = prediction_output,  # 写入这个字段
            status            = "pending",
        )
        session.merge(record)   # insert-or-update by primary key (token)
        session.commit()
    except Exception:
        session.rollback()
        raise
    finally:
        session.close()

def update_tad_outputs(token: str, output_paths: dict, prediction_output: str = None):
    """Persist caller output paths to the row identified by ``token``.

    For backward compatibility, the legacy four columns (output_di, output_sdi,
    output_is, output_tad_final) get populated with the first three v2 callers
    plus TAD_final when running under method_set=v2.  The full dict (which
    includes spectral_profile / Consensus / Disputed and any other keys) is
    JSON-serialized into ``output_paths_json``.
    """
    if not isinstance(output_paths, dict):
        output_paths = {}
    session = SessionLocal()
    try:
        record = session.query(TokenRecord).filter(TokenRecord.token == token).first()
        if record:
            # v2 first (insulation/topdom_like/contact_contrast), fall back to
            # legacy DI/TopDom/IS so the old path also keeps writing the same
            # columns it always did.
            record.output_di = (
                output_paths.get("insulation")
                or output_paths.get("DI")
            )
            record.output_sdi = (
                output_paths.get("topdom_like")
                or output_paths.get("TopDom")
                or output_paths.get("SDI")
            )
            record.output_is = (
                output_paths.get("contact_contrast")
                or output_paths.get("IS")
            )
            record.output_tad_final = output_paths.get("TAD_final")
            try:
                record.output_paths_json = json.dumps(output_paths)
            except (TypeError, ValueError):
                record.output_paths_json = None
            if prediction_output:
                record.prediction_output = prediction_output
            record.status           = "done"
            record.finished_at      = datetime.datetime.utcnow()
            session.commit()
    finally:
        session.close()

def mark_running(token: str):
    session = SessionLocal()
    try:
        record = session.query(TokenRecord).filter(TokenRecord.token == token).first()
        if record:
            record.status = "running"
            session.commit()
    finally:
        session.close()

def mark_failed(token: str, error_msg: str):
    session = SessionLocal()
    try:
        record = session.query(TokenRecord).filter(TokenRecord.token == token).first()
        if record:
            record.status      = "failed"
            record.error_msg   = error_msg
            record.finished_at = datetime.datetime.utcnow()
            session.commit()
    finally:
        session.close()

def _decode_output_paths(record: "TokenRecord") -> dict:
    """Prefer the JSON dict (v2 + future).  Fall back to the legacy 4 columns
    for rows persisted before the 2026-05-14 v2 switch."""
    raw = getattr(record, "output_paths_json", None)
    if raw:
        try:
            payload = json.loads(raw)
            if isinstance(payload, dict):
                return payload
        except (TypeError, ValueError):
            pass
    return {
        "DI":        record.output_di,
        "TopDom":    record.output_sdi,
        "IS":        record.output_is,
        "TAD_final": record.output_tad_final,
    }


def get_record(token: str) -> dict | None:
    session = SessionLocal()
    try:
        r = session.query(TokenRecord).filter(TokenRecord.token == token).first()
        if not r:
            return None
        return {
            "token":             r.token,
            "created_at":        r.created_at.isoformat() if r.created_at else None,
            "finished_at":       r.finished_at.isoformat() if r.finished_at else None,
            "mcool_path":        r.mcool_path,
            "resolution":        r.resolution,
            "is_threshold":      r.is_threshold,
            "overlap":           r.overlap,
            "hash_tag":          r.hash_tag,
            "fast_mode":         r.fast_mode,
            "status":            r.status,
            "error_msg":         r.error_msg,
            "prediction_output": r.prediction_output,
            "output_paths":      _decode_output_paths(r),
        }
    finally:
        session.close()

def get_all_records() -> list[dict]:
    session = SessionLocal()
    try:
        records = session.query(TokenRecord).all()
        result = []
        for r in records:
            result.append({
                "token":             r.token,
                "created_at":        r.created_at.isoformat() if r.created_at else None,
                "finished_at":       r.finished_at.isoformat() if r.finished_at else None,
                "mcool_path":        r.mcool_path,
                "resolution":        r.resolution,
                "is_threshold":      r.is_threshold,
                "overlap":           r.overlap,
                "hash_tag":          r.hash_tag,
                "fast_mode":         r.fast_mode,
                "status":            r.status,
                "error_msg":         r.error_msg,
                "prediction_output": r.prediction_output,
                "output_paths":      _decode_output_paths(r),
            })
        return result
    finally:
        session.close()

def delete_record(token: str) -> bool:
    session = SessionLocal()
    try:
        record = session.query(TokenRecord).filter(TokenRecord.token == token).first()
        if not record:
            return False
        session.delete(record)
        session.commit()
        return True
    except Exception:
        session.rollback()
        raise
    finally:
        session.close()

def delete_all_records() -> int:
    session = SessionLocal()
    try:
        count = session.query(TokenRecord).delete()
        session.commit()
        return count
    except Exception:
        session.rollback()
        raise
    finally:
        session.close()

# 新增：单独更新预测结果路径的函数（备用）
def update_prediction_output(token: str, prediction_output: str):
    session = SessionLocal()
    try:
        record = session.query(TokenRecord).filter(TokenRecord.token == token).first()
        if record:
            record.prediction_output = prediction_output
            session.commit()
    finally:
        session.close()