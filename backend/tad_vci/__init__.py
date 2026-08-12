"""TAD-VCI: structured evidence cards for post-calling TAD-boundary review."""
from .graded_evidence import (
    Criterion, CRITERIA, annotate, combine_tier, grade_loop,
    grade_fold, grade_votes, CTCF_BANDS, RAD21_BANDS,
)
# Retired BD5/conservation experiments are outside the evidence-card package contract.
__all__ = ["Criterion", "CRITERIA", "annotate", "combine_tier",
           "grade_loop", "grade_fold", "grade_votes",
           "CTCF_BANDS", "RAD21_BANDS"]
