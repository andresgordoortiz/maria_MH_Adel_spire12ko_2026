#!/usr/bin/env python3
"""Reproduce VAST config loading with Nextflow 25.10.4 in an isolated Podman container."""
from pathlib import Path
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]

with tempfile.TemporaryDirectory(prefix="spire-vast-startup-") as tmp:
    work = Path(tmp)
    for name in ("shared", "copied", "outputs"):
        (work / name).mkdir()
    (work / "shared/main.nf").write_text("workflow { log.info 'VAST_STARTUP_OK' }\n")
    (work / "shared/nextflow.config").write_text(
        "params { config_profile_description = 'VAST startup fixture' }\n"
        "env { SINGULARITYENV_VASTDB = params.vastdb_path }\n")
    shutil.copyfile(work / "shared/main.nf", work / "copied/main.nf")
    command = [
        "podman", "run", "--rm", "--read-only", "--network", "none", "--tmpfs", "/tmp",
        "--entrypoint", "/bin/bash",
        "-v", f"{ROOT}:/repo:ro", "-v", f"{work / 'shared'}:/shared:ro",
        "-v", f"{work / 'copied'}:/copied:ro", "-v", f"{work / 'outputs'}:/checks:rw",
        "-w", "/checks", "-e", "NXF_HOME=/checks/home", "-e", "NXF_OFFLINE=true",
        "-e", "PROJECT_DIR=/repo", "-e", "VAST_PIPELINE=/shared", "-e", "VASTDB=/vastdb",
        "-e", "SCRATCH_DIR=/scratch", "-e", "READONLY_BINDS=/repo:/repo:ro",
        "docker.io/nextflow/nextflow:25.10.4", "-c", '''set -uo pipefail
mkdir -p /checks/home /checks/direct /checks/copied
cp -r /.nextflow/framework /checks/home/
cd /checks/direct
nextflow -trace nextflow.cli.CmdRun run /shared/main.nf \
    --vastdb_path /vastdb -c /repo/config/spire_vast.config \
    -c /repo/config/fastq_crg.config -ansi-log false > stdout.log 2>&1
echo "$?" > status
cd /checks/copied
nextflow -trace nextflow.cli.CmdRun run /copied/main.nf \
    --vastdb_path /vastdb -c /repo/config/spire_vast.config \
    -c /repo/config/fastq_crg.config -ansi-log false > stdout.log 2>&1
echo "$?" > status
''',
    ]
    # Trace forces the same strict config resolution as Seqera, without contacting it.
    subprocess.run(command, check=True)
    direct = work / "outputs/direct"
    copied = work / "outputs/copied"
    assert (direct / "status").read_text().strip() != "0"
    assert "Unknown config attribute `params.vastdb_path`" in (direct / "stdout.log").read_text(), (direct / "stdout.log").read_text()
    assert (copied / "status").read_text().strip() == "0", (copied / "stdout.log").read_text()
    assert "VAST_STARTUP_OK" in (copied / "stdout.log").read_text()
print("VAST startup checks passed: direct launch reproduces the failure; the Murielle layout starts successfully.")
