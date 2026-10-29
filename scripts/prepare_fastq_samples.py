#!/usr/bin/env python3
"""Validate the SPIRE FASTQs and generate single-end Nextflow sample sheets."""
import argparse
import csv
from pathlib import Path
import re


def prepare(raw, outdir, trimmed=None):
    manifest = Path(__file__).resolve().parents[1] / "config/spire_fastq.tsv"
    with manifest.open() as handle:
        samples = list(csv.DictReader(handle, delimiter="\t"))
    ids = [row["sample"] for row in samples]
    filenames = [row["fastq_1"] for row in samples]
    if (len(samples) != 6 or len(set(ids)) != 6 or len(set(filenames)) != 6
            or any(not re.fullmatch(r"[A-Za-z0-9_]+", sample) for sample in ids)
            or sorted(row["condition"] for row in samples) != ["Cont"] * 3 + ["DKO"] * 3):
        raise ValueError("Expected six unique SPIRE samples: three Cont and three DKO")
    raw_files = [raw.resolve() / filename for filename in filenames]
    trimmed_files = ([trimmed.resolve() / f"{sample}.fastp.fastq.gz" for sample in ids]
                     if trimmed else [])
    for path in raw_files + trimmed_files:
        if not path.is_file() or not path.stat().st_size:
            raise ValueError(f"Missing or empty FASTQ: {path}")
    outdir.mkdir(parents=True, exist_ok=True)
    # shortcut: one FASTQ is treated as single-end; split interleaved pairs before use.
    with (outdir / "rnaseq.csv").open("w") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["sample", "fastq_1", "fastq_2", "strandedness"])
        writer.writerows((sample, path, "", "auto") for sample, path in zip(ids, raw_files))
    if trimmed:
        with (outdir / "vast.csv").open("w") as handle:
            writer = csv.writer(handle, lineterminator="\n")
            writer.writerow(["sample", "fastq_1", "fastq_2", "type", "group"])
            writer.writerows((row["sample"], path, "", "single", row["condition"])
                             for row, path in zip(samples, trimmed_files))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("raw", type=Path)
    parser.add_argument("outdir", type=Path)
    parser.add_argument("--trimmed", type=Path)
    args = parser.parse_args()
    prepare(args.raw, args.outdir, args.trimmed)
