"""Evidence cards and application-level versioned decision records.

Server-side and UI-agnostic, so it is fully unit-testable without a front end.
The application appends one JSON object per event and does not rewrite prior lines
through this API. Ordinary filesystem users can still edit, delete or reorder the
file. There is no identity authentication, hash chain, external anchor, locking or
immutable storage, so this module must not be described as signed or tamper-evident.
The checksum only detects a mismatch between the fields in one record and its saved
digest.
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Optional

from .graded_evidence import CRITERIA

_CARD_COLS = {
    "BD1": "BD1_caller_support", "BD2": "BD2_ctcf", "BD3": "BD3_rad21",
    "BD4": "BD4_loop_anchor",
}  # BD5 (CTCF motif) removed 2026-06-04


def _tier_why(grades: dict) -> str:
    """One-line, human-readable derivation of max_support_v2 for this boundary."""
    from tad_vci.graded_evidence import combine_tier, _rank, _RANK
    bd1 = grades.get("BD1", "not_assessable")
    try:
        tier, meta = combine_tier(grades)
    except ValueError:
        return "BD1 not assessable — cannot form a tier"
    assessed_sup = [c for c in meta.get("assessed", []) if c != "BD1"]
    best_name = "none"
    best_r = 0
    for c in assessed_sup:
        g = grades.get(c, "not_assessable")
        if g != "not_assessable" and _rank(g) >= best_r:
            best_r = _rank(g)
            best_name = f"{c}={g}"
    # Mirror combine_tier branches in plain language
    n_assessed_sup = len(assessed_sup)
    if bd1 == "strong" and best_r >= _RANK["strong"]:
        rule = "BD1 strong + supporting strong → D5"
    elif bd1 == "strong" and n_assessed_sup and best_r == 0:
        # max_support_v2 depletion clause. Without this branch the string said "→ D4"
        # for a boundary the rule actually puts at D3, and that string is written into
        # the signed decision record.
        # Wording corrected 2026-07-29 with the v2 band ladder: on a ChIP axis `none` now
        # means "below the p75 of this cell's random-window null", i.e. not enriched — not
        # "depleted below background", which is what this string used to assert.
        rule = (f"BD1 strong but no measured supporting axis reaches its weak band "
                f"({', '.join(assessed_sup)} all none) → D3")
    elif bd1 == "strong":
        rule = ("BD1 strong, no supporting axis measured → D4"
                if not n_assessed_sup else
                "BD1 strong (no supporting strong) → D4")
    elif bd1 == "moderate" and best_r >= _RANK["moderate"]:
        rule = "BD1 moderate + supporting ≥moderate → D4"
    elif bd1 == "moderate":
        rule = "BD1 moderate (supporting < moderate) → D3"
    elif bd1 == "weak" and best_r >= _RANK["moderate"]:
        rule = "BD1 weak + supporting ≥moderate → D3"
    elif bd1 == "weak" and best_r >= _RANK["weak"]:
        rule = "BD1 weak + supporting weak → D2"
    else:
        rule = "BD1 weak + no positive supporting → D1"
    sup = best_name if best_r > 0 else "no supporting signal"
    return f"{rule} (best supporting: {sup})"


def _num(row: dict, key: str):
    """Float from an annotation row, or None for missing / NaN / unparseable."""
    v = row.get(key)
    if v is None:
        return None
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return None if f != f else f


def _measurement(row: dict, fold_key: str, bg_key: str, pct_key: str,
                 bands) -> Optional[dict]:
    """The numbers behind a BD2/BD3 grade: observed fold, background, cutpoints, null rank.

    A grade word alone is not inspectable — `weak` covered folds from 1.0 to the moderate
    cut, and under the retired v1 ladder its bottom half was indistinguishable from a
    random window. The card must therefore carry what was measured, what it was divided
    by, where the cutpoints are, and (when the cell's null is on file) what fraction of
    random +/-25 kb windows in this cell reach the observed fold."""
    fold = _num(row, fold_key)
    if fold is None:
        return None
    bd = dict(bands or ())
    pct = _num(row, pct_key)
    out = {
        "statistic": "mean signal over +/-25 kb / random-window background",
        "fold": round(fold, 4),
        "background_mean": (round(bg, 5) if (bg := _num(row, bg_key)) is not None else None),
        "window_bp": int(_num(row, "chip_window_bp") or 50000),
        "cutpoints": {name: round(float(cut), 4) for name, cut in bd.items()},
        "cutpoint_percentiles": {"weak": 75, "moderate": 90, "strong": 97.5},
    }
    if pct is not None:
        out["null_percentile"] = round(pct, 2)
        out["random_windows_at_or_above"] = round(100.0 - pct, 2)
    return out


def evidence_card(row: dict) -> dict:
    """Structured, human+machine readable evidence card for one boundary."""
    import math as _math
    criteria = {}
    for code, col in _CARD_COLS.items():
        c = CRITERIA[code]
        criteria[code] = {
            "name": c.name, "availability": c.availability, "role": c.role,
            "requires": c.assay, "grading": c.grading,
            "grade": row.get(col, "not_assessable"),
        }
    # BD1 continuous companion: the insulation boundary-strength (depth of the insulation
    # minimum from the Hi-C map). BD1's vote count is discrete; this quantifies HOW sharp the
    # boundary is. Informational refinement of BD1, not a separate criterion / not in the tier.
    _ins = row.get("insulation_strength")
    if _ins is not None and not (isinstance(_ins, float) and _ins != _ins):
        try:
            criteria["BD1"]["insulation_strength"] = round(float(_ins), 3)
        except (TypeError, ValueError):
            pass
    _anch = _num(row, "loop_anchor_count")
    if _anch is not None:
        # BD4's window is its own parameter (published as 50 kb). Older tables predate the
        # split and carry only vote_tol_bp, which is what they actually used.
        criteria["BD4"]["measured"] = {
            "loop_anchor_count": int(_anch),
            "window_bp": int(_num(row, "loop_window_bp") or _num(row, "vote_tol_bp") or 50000),
            "cutpoints": {"strong": 2, "moderate": 1, "none": 0},
        }
    assessed = [k for k, v in criteria.items() if v["grade"] != "not_assessable"]
    bid = f"{row['chrom']}:{int(row['pos'])}"
    grade_summary = ", ".join(
        f"{criterion}={criteria[criterion]['grade']}" for criterion in assessed
    )
    summary = (f"{bid}: evidence-tier {row.get('D_tier')} from "
               f"{grade_summary}"
               f"; not assessable: {', '.join(k for k in criteria if k not in assessed) or 'none'}")
    ub = row.get("user_bed")
    # BD2/BD3 per-cell band-calibration disclosure (absent column -> treat as calibrated)
    bcal = row.get("bands_calibrated")
    bcal_ok = True if bcal is None or (isinstance(bcal, float) and bcal != bcal) else bool(bcal)
    dcell = row.get("data_cell_line")
    dcell = dcell if isinstance(dcell, str) and dcell and dcell != "nan" else None
    bands_disclosure = (None if bcal_ok else
                        f"BD2/BD3 thresholds use GM12878 calibration — NOT calibrated for "
                        f"{dcell or 'this cell line'}; CTCF/RAD21 grades are approximate.")
    # Per-cell BD2/BD3 grading TEXT: the static registry string quotes GM12878's thresholds;
    # rewrite the numbers with THIS cell's calibrated cutpoints so a reviewer reading e.g. a K562
    # card sees K562's thresholds (the grade itself already uses the correct per-cell bands).
    #
    # Two silent failures used to hide here. (a) The whole block was wrapped in
    # `except Exception: pass`, so a band-table failure left the static GM12878 numbers on
    # the card with nothing saying so. (b) When load_bands could not calibrate this cell it
    # returns the GM12878 placeholder with calibrated=False, yet the text was still written
    # as "per-cell GM12878" — a card for an uncalibrated cell line claimed a per-cell
    # threshold it never had, while bands_disclosure stayed None because the annotation had
    # no bands_calibrated column. State the provenance of the numbers instead of implying it.
    try:
        from tad_vci.graded_evidence import load_bands as _lb
        _ctcf_b, _rad21_b, _cal, _used = _lb(dcell)
        _prov = (f"per-cell {_used}" if _cal else
                 f"{_used} reference thresholds, NOT calibrated for "
                 f"{dcell or 'this cell line'}")

        def _fmt(axis, bands):
            bd = dict(bands)
            return (f"+/-25 kb mean {axis} fold vs background ({_prov}): "
                    f">={bd.get('strong'):.3f} strong (p97.5 of random windows), "
                    f">={bd.get('moderate'):.3f} moderate (p90), "
                    f">={bd.get('weak'):.3f} weak (p75), below that none")
        # The measurement block is omitted, not set to null, when the annotation carries no
        # fold: an empty block on the card would read as "measured, and it was nothing".
        for code, fold_key, bg_key, pct_key, bands, axis in (
            ("BD2", "ctcf_fold", "bg_ctcf", "ctcf_null_pct", _ctcf_b, "CTCF"),
            ("BD3", "rad21_fold", "bg_rad21", "rad21_null_pct", _rad21_b, "RAD21"),
        ):
            if code not in criteria:
                continue
            criteria[code]["grading"] = _fmt(axis, bands)
            meas = _measurement(row, fold_key, bg_key, pct_key, bands)
            if meas is not None:
                criteria[code]["measured"] = meas
        if not _cal and bands_disclosure is None:
            bands_disclosure = (
                f"BD2/BD3 thresholds use {_used} calibration — NOT calibrated for "
                f"{dcell or 'this cell line'}; CTCF/RAD21 grades are approximate.")
    except Exception as exc:
        bands_disclosure = (
            "BD2/BD3 threshold text could not be resolved for this cell line; the numbers "
            "shown are the GM12878 reference thresholds from the criteria registry.")
        print(f"[evidence-card] per-cell band text unavailable for "
              f"{dcell or 'unnamed cell line'}: {type(exc).__name__}: {exc}")

    # BD1 panel context: votes/n_methods + cutpoints (the usual "3 methods but still weak" trap)
    votes = None
    n_methods = None
    try:
        if row.get("votes") is not None and not (isinstance(row.get("votes"), float) and row["votes"] != row["votes"]):
            votes = int(row["votes"])
    except (TypeError, ValueError):
        votes = None
    try:
        if row.get("n_methods") is not None and not (isinstance(row.get("n_methods"), float) and row["n_methods"] != row["n_methods"]):
            n_methods = int(row["n_methods"])
    except (TypeError, ValueError):
        n_methods = None
    if n_methods is None or n_methods < 1:
        n_methods = 5  # built-in panel default when annotation lacks the column
    mod_cut = int(_math.ceil(0.6 * n_methods))
    str_cut = int(_math.ceil(0.8 * n_methods))
    if n_methods <= 1:
        bd1_threshold = (
            "BD1 = cross-caller support. With a single-method panel (n=1) "
            "cross-caller agreement cannot be assessed, so BD1 is floored to weak "
            "(the boundary is supported only by the one caller that called it)."
        )
    else:
        bd1_threshold = (
            f"BD1 uses fraction of panel n={n_methods}: "
            f"strong if votes≥{str_cut} (≥80%), moderate if votes≥{mod_cut} (≥60%), else weak. "
            f"Seeing 2–3 method tracks is still weak when n={n_methods}."
        )
    sc = row.get("supporting_callers")
    if isinstance(sc, float) and sc != sc:
        sc = None
    supporting_callers = str(sc) if sc else None

    # BD1 POSITION-MATCHING TOLERANCE. A method votes when ANY of its calls falls within
    # vote_tol_bp of the candidate, so "2/6 methods" does NOT mean two callers drew a block
    # at this position. Until 2026-07-29 the card said nothing about this while it spelled
    # out the +/-25 kb window for BD2/BD3 and the 50 kb window for BD4 — the one criterion
    # with the largest influence on the tier was the one with no window on the card.
    tol_bp = _num(row, "vote_tol_bp")
    res_bp = _num(row, "resolution_bp")
    vote_model = row.get("bd1_vote_model")
    vote_model = str(vote_model) if isinstance(vote_model, str) and vote_model else None
    offs_raw = row.get("supporting_offsets_bp")
    if isinstance(offs_raw, float) and offs_raw != offs_raw:
        offs_raw = None
    supporters: list[dict] = []
    if supporting_callers:
        names = supporting_callers.split(";")
        offs = str(offs_raw).split(";") if offs_raw not in (None, "") else []
        for i, nm in enumerate(names):
            entry = {"caller": nm}
            if i < len(offs):
                try:
                    o = int(float(offs[i]))
                except (TypeError, ValueError):
                    o = None
                if o is not None:
                    entry["offset_bp"] = o
                    if res_bp:
                        entry["offset_bins"] = round(o / res_bp, 2)
            supporters.append(entry)
    exact = [s for s in supporters if s.get("offset_bp") == 0]
    bd1_match = None
    if tol_bp:
        bins = f" (+/-{tol_bp / res_bp:g} bin{'s' if res_bp and tol_bp / res_bp != 1 else ''})" if res_bp else ""
        bd1_match = (
            f"A method votes if any of its calls lies within +/-{int(tol_bp):,} bp{bins} "
            f"of this position; a vote is not a call at this position."
        )
        if vote_model == "nearest_within_tol":
            bd1_match += (" Candidates are the raw union of caller calls, so one call inside "
                          "that window supports every candidate in it: votes at neighbouring "
                          "candidates are correlated, not independent observations.")
        elif vote_model == "cluster_consensus":
            bd1_match += (" Positions were clustered within the same window first, so each "
                          "method contributes at most one vote per candidate.")
    criteria["BD1"]["measured"] = {
        "votes": votes,
        "n_methods": n_methods,
        "vote_tol_bp": int(tol_bp) if tol_bp else None,
        "vote_tol_bins": (round(tol_bp / res_bp, 2) if (tol_bp and res_bp) else None),
        "vote_model": vote_model,
        "supporters": supporters or None,
        "n_supporters_at_this_position": len(exact) if supporters else None,
    }

    grades_now = {
        "BD1": criteria["BD1"]["grade"],
        "BD2": criteria["BD2"]["grade"],
        "BD3": criteria["BD3"]["grade"],
        "BD4": criteria["BD4"]["grade"],
    }
    tier_why = _tier_why(grades_now)

    return {
        "boundary_id": bid,
        "auto_tier": row.get("D_tier"),
        "bands_disclosure": bands_disclosure,   # None if BD2/BD3 bands calibrated for the data's cell line
        # True when this candidate came from a user-uploaded BED (other algorithm),
        # folded into the tier and graded by the same evidence framework.
        "user_bed": bool(ub) if ub is not None and not (isinstance(ub, float) and ub != ub) else False,
        "criteria": criteria,
        "assessed": assessed,
        "summary": summary,
        "votes": votes,
        "n_methods": n_methods,
        "bd1_strong_at": str_cut,
        "bd1_moderate_at": mod_cut,
        "bd1_threshold": bd1_threshold,
        "bd1_match_rule": bd1_match,
        "vote_tol_bp": int(tol_bp) if tol_bp else None,
        "supporting_callers": supporting_callers,
        "supporters": supporters or None,
        "tier_why": tier_why,
        "tier_rule": "max_support_v2",
    }


def _canonical(obj) -> str:
    return json.dumps(obj, sort_keys=True, separators=(",", ":"))


def content_hash(payload: dict) -> str:
    """Full SHA-256 checksum for one canonical record payload."""
    return hashlib.sha256(_canonical(payload).encode("utf-8")).hexdigest()


def _legacy_content_hash(payload: dict) -> str:
    """Verify pre-July records that used a 64-bit truncated checksum."""
    return hashlib.sha256(_canonical(payload).encode("utf-8")).hexdigest()[:16]


@dataclass(frozen=True)
class CurationSnapshot:
    boundary_id: str
    assembly: str
    version: int
    auto_tier: Optional[str]
    final_tier: Optional[str]          # reviewed tier, else auto_tier
    evidence: dict
    review_decision: Optional[dict]    # {"final_tier","verdict","rationale"} or None
    actor_label: str                   # caller-supplied label; not authenticated identity
    created_at: str                    # ISO timestamp, supplied by caller
    content_hash: str
    hash_schema: str = "creditad-decision-checksum-v3"
    resolution: Optional[int] = None
    cell_line: Optional[str] = None

    def to_json(self) -> str:
        return _canonical(asdict(self))


def make_snapshot(card: dict, assembly: str, actor_label: str, created_at: str,
                  decision: Optional[dict] = None, prior_version: int = 0,
                  resolution: Optional[int] = None,
                  cell_line: Optional[str] = None) -> CurationSnapshot:
    """Build a versioned record from an unauthenticated local review label."""
    actor_label = str(actor_label).strip()
    if not actor_label:
        raise ValueError("actor_label is required (it is not an authenticated identity)")
    if decision is not None:
        if not decision.get("rationale", "").strip():
            raise ValueError("review decision requires a non-empty rationale")
        if not decision.get("final_tier"):
            raise ValueError("review decision requires a final_tier")
        final_tier = decision["final_tier"]
    else:
        final_tier = card["auto_tier"]
    version = prior_version + 1
    payload = {
        "boundary_id": card["boundary_id"], "assembly": assembly,
        "auto_tier": card["auto_tier"], "final_tier": final_tier,
        "evidence": card["criteria"], "review_decision": decision,
        "version": version, "actor_label": actor_label, "created_at": created_at,
        "resolution": resolution, "cell_line": cell_line,
        "hash_schema": "creditad-decision-checksum-v3",
    }
    return CurationSnapshot(
        boundary_id=card["boundary_id"], assembly=assembly, version=version,
        auto_tier=card["auto_tier"], final_tier=final_tier, evidence=card["criteria"],
        review_decision=decision, actor_label=actor_label, created_at=created_at,
        content_hash=content_hash(payload),
        resolution=resolution, cell_line=cell_line,
    )


def append_snapshot(snap: CurationSnapshot, log_path: str | Path) -> None:
    """Append through the application API; this is not an immutable-storage guarantee."""
    p = Path(log_path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("a", encoding="utf-8") as f:
        f.write(snap.to_json() + "\n")


_ASSEMBLY_ALIASES = {"grch38": "hg38", "hg38": "hg38", "grch37": "hg19", "hg19": "hg19"}


def _norm_cell(value) -> Optional[str]:
    if not isinstance(value, str):
        return None
    v = value.strip()
    return v.upper() if v and v.lower() != "nan" else None


def _norm_asm(value) -> Optional[str]:
    if not isinstance(value, str):
        return None
    v = value.strip().lower()
    if not v or v == "nan":
        return None
    return _ASSEMBLY_ALIASES.get(v, v)


def _norm_res(value) -> Optional[int]:
    try:
        if value is None or value == "":
            return None
        return int(value)
    except (TypeError, ValueError):
        return None


def record_scope(rec: dict) -> dict:
    """Normalised dataset scope of one record: which map it was actually reviewed on."""
    return {
        "cell_line": _norm_cell(rec.get("cell_line")),
        "assembly": _norm_asm(rec.get("assembly")),
        "resolution": _norm_res(rec.get("resolution")),
    }


def parse_boundary_selector(selector: str) -> tuple[str, Optional[dict]]:
    """Split "chr7:1000000@GM12878/hg38/10000" into (boundary_id, scope).

    A boundary id is only chrom:pos, and that string collides across datasets: the
    same coordinate exists in GM12878 hg38 10 kb and in K562 hg38 25 kb. The optional
    "@cell/assembly/resolution" suffix names WHICH dataset the caller is asking about.
    A bare id returns scope None, which means "every scope" (used by the export /
    all-records view and by pre-scope callers).
    """
    if "@" not in selector:
        return selector, None
    boundary_id, _, tail = selector.partition("@")
    parts = tail.split("/")
    if len(parts) != 3:
        # Not a scope suffix we understand — treat the whole thing as an id rather
        # than silently widening the query to every cell line.
        return selector, None
    return boundary_id, {
        "cell_line": _norm_cell(parts[0]),
        "assembly": _norm_asm(parts[1]),
        "resolution": _norm_res(parts[2]),
    }


def scope_matches(rec: dict, scope: Optional[dict]) -> bool:
    """True when `rec` was recorded on exactly the dataset named by `scope`.

    Matching is strict on all three fields: a record with no cell_line (a pre-v3
    legacy line) is NOT shown under a cell-scoped query, because attributing it to a
    cell line would be a guess. Those records remain visible through an unscoped
    query and through the all-records view.
    """
    if scope is None:
        return True
    have = record_scope(rec)
    for key in ("cell_line", "assembly", "resolution"):
        if key in scope and have[key] != scope[key]:
            return False
    return True


def load_history(boundary_id: str, log_path: str | Path,
                 scope: Optional[dict] = None) -> list[dict]:
    """Records for one boundary, optionally restricted to one dataset scope.

    `boundary_id` may carry an "@cell/assembly/resolution" suffix (see
    parse_boundary_selector); an explicit `scope` argument wins over the suffix.
    """
    parsed_id, parsed_scope = parse_boundary_selector(boundary_id)
    effective_scope = scope if scope is not None else parsed_scope
    p = Path(log_path)
    if not p.exists():
        return []
    out = []
    for line in p.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        rec = json.loads(line)
        if rec.get("boundary_id") != parsed_id:
            continue
        if not scope_matches(rec, effective_scope):
            continue
        out.append(rec)
    return sorted(out, key=lambda r: r.get("version", 0))


def next_version(boundary_id: str, log_path: str | Path) -> int:
    """Highest version already written for this boundary id, across ALL scopes.

    Deliberately unscoped: version numbers must never be reused for one
    boundary_id, otherwise two records with the same id and version but different
    content both claim to be "v1" in an export. Per-scope history therefore starts
    at whatever version the log had reached, which is expected.
    """
    hist = load_history(parse_boundary_selector(boundary_id)[0], log_path)
    return (hist[-1]["version"] if hist else 0)


def verify_snapshot(rec: dict) -> bool:
    """Recompute the per-record checksum; does not prove log completeness or identity."""
    schema = rec.get("hash_schema")
    if schema == "creditad-decision-checksum-v3":
        payload = {
            "boundary_id": rec["boundary_id"], "assembly": rec["assembly"],
            "auto_tier": rec["auto_tier"], "final_tier": rec["final_tier"],
            "evidence": rec["evidence"], "review_decision": rec["review_decision"],
            "version": rec["version"], "actor_label": rec["actor_label"],
            "created_at": rec["created_at"], "resolution": rec.get("resolution"),
            "cell_line": rec.get("cell_line"), "hash_schema": schema,
        }
        return content_hash(payload) == rec["content_hash"]

    # Backward verification only. New records never write expert/curator fields.
    # A verifier must answer true/false, never raise: a v3 record whose hash_schema has
    # been rewritten to a legacy value has no `expert_override` key, and the KeyError that
    # produced would propagate out of the integrity check instead of failing it.
    try:
        legacy_payload = {
            "boundary_id": rec["boundary_id"], "assembly": rec["assembly"],
            "auto_tier": rec["auto_tier"], "final_tier": rec["final_tier"],
            "evidence": rec["evidence"], "expert_override": rec["expert_override"],
            "version": rec["version"],
        }
    except KeyError:
        return False
    if schema != "creditad-decision-checksum-v2":
        return _legacy_content_hash(legacy_payload) == rec["content_hash"]
    payload = {
        **legacy_payload,
        "curator": rec["curator"],
        "created_at": rec["created_at"],
        "resolution": rec.get("resolution"),
        "cell_line": rec.get("cell_line"),
        "hash_schema": schema,
    }
    return content_hash(payload) == rec["content_hash"]
