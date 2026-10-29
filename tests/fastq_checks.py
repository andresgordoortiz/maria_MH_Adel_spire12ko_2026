#!/usr/bin/env python3
"""Check SPIRE sheets, missing-input failures and the Slurm/Nextflow hand-off."""
import csv
import gzip
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]


with tempfile.TemporaryDirectory(prefix="spire fastq ") as tmp:
    project = Path(tmp) / "project"
    for directory in ("scripts", "config", "bin", "raw", "vast", "db/Mm2/FILES"):
        (project / directory).mkdir(parents=True)
    for name in ("prepare_fastq_samples.py", "submit_fastq.sh", "fastq_pipeline.sbatch",
                 "fastq_betas.sbatch", "run_R.sh"):
        shutil.copy(ROOT / "scripts" / name, project / "scripts" / name)
    for path in (ROOT / "config").glob("*"):
        shutil.copy(path, project / "config" / path.name)
    with (project / "config/spire_fastq.tsv").open() as handle:
        manifest = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["fastq_1"] for row in manifest] == [
        "reads_1-controle-12mars.fastq.gz", "reads_2-controle-12mars.fastq.gz",
        "reads_7-controle-14mars.fastq.gz", "reads_4-doubleKO-12mars.fastq.gz",
        "reads_5-doubleKO-12mars.fastq.gz", "reads_6-doubleKO-12mars.fastq.gz"]
    for row in manifest:
        with gzip.open(project / "raw" / row["fastq_1"], "wt") as handle:
            handle.write("@read\nACGT\n+\nIIII\n")
    for name in ("genome.fa", "genes.gtf", "vast/main.nf", "vast/nextflow.config", "db/Mm2/FILES/gDNA.1.ebwt"):
        (project / name).write_text("fixture\n")
    log = project / "calls.jsonl"
    stub = f"#!{sys.executable}\n" + '''import json, os, pathlib, sys
log = pathlib.Path(os.environ["CALL_LOG"])
previous = log.read_text().splitlines() if log.exists() else []
with log.open("a") as handle:
    handle.write(json.dumps({"tool": pathlib.Path(sys.argv[0]).name,
        "args": sys.argv[1:], "cwd": os.getcwd(),
        "binds": os.environ.get("READONLY_BINDS", ""),
        "dataset": os.environ.get("BETAS_DATASET", ""),
        "cache": os.environ.get("NXF_SINGULARITY_CACHEDIR", ""),
        "tmpdir": os.environ.get("TMPDIR", "")}) + "\\n")
if pathlib.Path(sys.argv[0]).name == "sbatch":
    print(str(101 + len(previous)) + ";crg")
if pathlib.Path(sys.argv[0]).name == "Rscript" and sys.argv[1] == "scripts/00_init_run.R":
    run = pathlib.Path(sys.argv[2])
    run.mkdir(parents=True)
    (run / "run.rds").write_text("fixture")
'''
    for name in ("sbatch", "nextflow", "Rscript", "module", "singularity"):
        path = project / "bin" / name
        path.write_text(stub if name in ("sbatch", "nextflow", "Rscript") else "#!/bin/sh\nexit 0\n")
        path.chmod(0o755)
    broken_python = project / "bin/python3"
    broken_python.write_text("#!/bin/sh\necho 'Broken login Python: missing encodings' >&2\nexit 97\n")
    broken_python.chmod(0o755)
    module_bin = project / "python_module/bin"
    module_bin.mkdir(parents=True)
    (module_bin / "python3").symlink_to(sys.executable)
    module_init = project / "modules.sh"
    module_init.write_text('''module() {
    if [[ $1 == load && $2 == Python/3.12.3-GCCcore-13.2.0 ]]; then
        export PATH="$TEST_PYTHON_MODULE_BIN:$PATH"
    fi
}
''')
    env = dict(os.environ, PATH=f"{project / 'bin'}:{os.environ['PATH']}",
               PROJECT_DIR=str(project), RAW_DIR=str(project / "raw"),
               VAST_PIPELINE=str(project / "vast"), VASTDB=str(project / "db"),
               SCRATCH_DIR=str(Path(tmp) / "scratch"),
               GENOME_FASTA=str(project / "genome.fa"), GENES_GTF=str(project / "genes.gtf"),
               BASH_ENV=str(module_init), TEST_PYTHON_MODULE_BIN=str(module_bin),
               RSCRIPT=str(project / "bin/Rscript"), SLURM_JOB_ID="12345", CALL_LOG=str(log))
    env.pop("TOWER_ACCESS_TOKEN", None)
    env.pop("SPIRE_FASTQ_INCLUSION_TABLE", None)

    def run(*args, success=True):
        result = subprocess.run(args, cwd=project, env=env, capture_output=True, text=True)
        assert (result.returncode == 0) == success, result.stdout + result.stderr
        return result

    def calls():
        return [json.loads(line) for line in log.read_text().splitlines()]

    assert "Broken login Python" in run("python3", "-I", "-c", "import encodings", success=False).stderr
    run("bash", "scripts/submit_fastq.sh")
    submitted = calls()
    assert len(submitted) == 4
    assert "--dependency=afterok:101" in submitted[1]["args"]
    assert "--dependency=afterok:101" in submitted[2]["args"]
    assert [call["args"][-2] for call in submitted[:3]] == ["rnaseq", "vast", "dge"]
    assert "--dependency=afterok:102" in submitted[3]["args"]
    assert submitted[3]["args"][-2] == "scripts/fastq_betas.sbatch"
    sheets = project / "data/processed/spire_fastq/samples"
    assert not sheets.exists()
    missing = project / "raw" / manifest[0]["fastq_1"]
    contents = missing.read_bytes()
    missing.write_bytes(b"")
    assert "Missing or empty FASTQ" in run("bash", "scripts/submit_fastq.sh", success=False).stderr
    assert len(calls()) == 4 and not sheets.exists()
    missing.write_bytes(contents)
    manifest_path = project / "config/spire_fastq.tsv"
    original_manifest = manifest_path.read_text()
    manifest_path.write_text(original_manifest.replace("SPIRE_Cont_2", "SPIRE_Cont_1"))
    assert "six unique SPIRE samples" in run("bash", "scripts/fastq_pipeline.sbatch", "rnaseq", str(project), success=False).stderr
    assert len(calls()) == 4
    manifest_path.write_text(original_manifest)
    for stage in ("rnaseq", "vast", "dge"):
        if stage == "vast":
            run("bash", "scripts/fastq_pipeline.sbatch", stage, str(project), success=False)
            assert (sheets / "rnaseq.csv").read_bytes() == original
            trimmed = project / "results/spire_fastq/rnaseq/fastp"
            trimmed.mkdir(parents=True)
            for row in manifest:
                shutil.copy(missing, trimmed / f"{row['sample']}.fastp.fastq.gz")
        if stage == "dge":
            run("bash", "scripts/fastq_pipeline.sbatch", stage, str(project), success=False)
            matrixdir = project / "results/spire_fastq/rnaseq/star_salmon"
            matrixdir.mkdir()
            for matrix in ("counts", "lengths"):
                (matrixdir / f"salmon.merged.gene_{matrix}.tsv").write_text("fixture\n")
        run("bash", "scripts/fastq_pipeline.sbatch", stage, str(project))
        if stage == "rnaseq":
            with (sheets / "rnaseq.csv").open() as handle:
                rows = list(csv.DictReader(handle))
            assert len(rows) == 6 and all(row["fastq_2"] == "" for row in rows)
            assert all(row["strandedness"] == "auto" for row in rows)
            original = (sheets / "rnaseq.csv").read_bytes()
    commands = calls()[4:]
    assert len(commands) == 3
    assert commands[0]["args"][:4] == ["run", "nf-core/rnaseq", "-r", "3.27.0"]
    assert all(flag in commands[1]["args"] for flag in ("--skip_trimming", "--skip_compare", "--skip_betas"))
    assert commands[2]["args"][:4] == ["run", "nf-core/differentialabundance", "-r", "2.0.0"]
    assert "--feature_length_matrix" in commands[2]["args"]
    assert str(project / "genes.gtf") in commands[0]["args"]
    assert str(project / "genes.gtf") in commands[2]["args"]
    assert all("-resume" in call["args"] and ":ro" in call["binds"] for call in commands)
    assert all(call["cache"] == str(Path(tmp) / "scratch/singularity") for call in commands)
    for stage, call in zip(("rnaseq", "vast", "dge"), commands):
        assert str(Path(tmp) / f"scratch/maria_MH_Adel_spire12ko_2026/work_{stage}") in call["args"]
    with (sheets / "vast.csv").open() as handle:
        rows = list(csv.DictReader(handle))
    assert [row["sample"] for row in rows] == [row["sample"] for row in manifest]
    assert all(row["type"] == "single" and row["fastq_2"] == "" for row in rows)
    assert [row["group"] for row in rows] == ["Cont"] * 3 + ["DKO"] * 3
    with (project / "config/spire_contrasts.csv").open() as handle:
        contrast, = csv.DictReader(handle)
    assert (contrast["reference"], contrast["target"], contrast["blocking"]) == ("Cont", "DKO", "")

    # Resuming betAS needs the inclusion table, not raw FASTQs or RNA-seq references.
    shutil.rmtree(project / "raw")
    (project / "genome.fa").unlink()
    (project / "genes.gtf").unlink()
    assert "Missing VAST inclusion table" in run("bash", "scripts/submit_fastq.sh", "betas", success=False).stderr
    inclusion = project / "results/spire_fastq/vast/inclusion_tables/spire_fastq_INCLUSION_LEVELS_FULL-mm10.tab"
    inclusion.parent.mkdir(parents=True)
    inclusion.write_text("fixture\n")
    run("bash", "scripts/submit_fastq.sh", "betas")
    assert not any(arg.startswith("--dependency") for arg in calls()[-1]["args"])
    start = len(calls())
    run("bash", "scripts/fastq_betas.sbatch", str(project))
    first = calls()[start:]
    assert [call["args"][0] for call in first] == [
        "scripts/01_filter.R", "scripts/00_init_run.R", "scripts/02_simulate.R", "scripts/03_collect.R"]
    assert all(call["dataset"] == "spire_fastq" for call in first)
    assert all(call["tmpdir"] == str(Path(tmp) / "scratch/maria_MH_Adel_spire12ko_2026/betas/12345") for call in first)
    start = len(calls())
    run("bash", "scripts/fastq_betas.sbatch", str(project))
    assert [call["args"][0] for call in calls()[start:]] == ["scripts/02_simulate.R", "scripts/03_collect.R"]
print("SPIRE FASTQ checks passed (broken login Python, module-loaded workers, mock Slurm/Nextflow).")
