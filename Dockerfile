# CrediTAD headline-number reproduction image.
# Recomputes every reported figure from the shipped ~6 MB intermediates and
# exits 0 iff all checks pass. Does NOT bundle the GB-scale raw Hi-C; the full
# raw-data pipeline is documented in REPRODUCE.md.
#
#   docker build -t creditad-reproduce .
#   docker run --rm creditad-reproduce
FROM python:3.12-slim

WORKDIR /creditad
COPY reproduce/ ./reproduce/
RUN pip install --no-cache-dir -r reproduce/requirements.txt

LABEL org.opencontainers.image.title="CrediTAD headline-number reproduction"
LABEL org.opencontainers.image.description="Recomputes all reported figures from shipped intermediates; exits 0 iff all checks pass"

CMD ["python", "reproduce/recompute_headline_numbers.py"]
