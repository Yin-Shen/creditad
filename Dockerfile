# CrediTAD command-line interface and quick-start check in a clean Python environment.
#
#   docker build -t creditad:1.0.4 .
#   docker run --rm creditad:1.0.4 --help
#   docker run --rm --entrypoint python creditad:1.0.4 examples/quickstart/run_quickstart.py
#   docker run --rm -v "$PWD":/work -w /work creditad:1.0.4 annotate --boundaries ... --out out
#
# The image carries the package with its test dependencies, the tests and the data shipped
# in this repository (`docker run --rm --entrypoint pytest creditad:1.0.4 -q` runs the suite). The
# genome-wide ENCODE ChIP tracks needed by scripts/reproduce_packages.py are not included
# (see REPRODUCE.md); mount them and set CREDITAD_DATA_ROOT.
FROM python:3.12-slim

WORKDIR /opt/creditad
COPY . .
RUN pip install --no-cache-dir '.[test]'

LABEL org.opencontainers.image.title="CrediTAD" \
      org.opencontainers.image.version="1.0.4" \
      org.opencontainers.image.source="https://github.com/Yin-Shen/creditad" \
      org.opencontainers.image.licenses="MIT"

ENTRYPOINT ["creditad"]
