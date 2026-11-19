#!/usr/bin/env python3
"""Combine SAIGEQTL Step 2 files into the wide table consumed by Step 3."""

import argparse
import os
import re
import sys
from pathlib import Path

os.environ.setdefault("POLARS_MAX_THREADS", "1")

def parse_args():
    parser = argparse.ArgumentParser(
        description="Concatenate Step 2 main-effect results and filter by MAF."
    )
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--file-pattern", default="*", help="Input glob [default: *]")
    parser.add_argument(
        "--input-format",
        choices=["auto", "txt", "parquet"],
        default="auto",
        help="Input format. 'auto' detects .parquet by extension and treats other files as tab-delimited text [default: auto].",
    )
    parser.add_argument(
        "--gene-regex",
        required=True,
        help="Regex applied to each filename; it must contain a named group '(?P<gene>...)'.",
    )
    parser.add_argument("--maf-min", type=float, default=0.05)
    parser.add_argument("--maf-max", type=float, default=0.95)
    parser.add_argument(
        "--threads", type=int, default=1, help="Polars worker threads [default: 1]"
    )
    return parser.parse_args()


def main():
    args = parse_args()
    os.environ["POLARS_MAX_THREADS"] = str(args.threads)
    import polars as pl

    if not 0 <= args.maf_min < args.maf_max <= 1:
        raise SystemExit("Require 0 <= --maf-min < --maf-max <= 1")

    gene_regex = re.compile(args.gene_regex)
    if "gene" not in gene_regex.groupindex:
        raise SystemExit("--gene-regex must contain a named '(?P<gene>...)' group")

    files = sorted(path for path in args.input_dir.glob(args.file_pattern) if path.is_file())
    if not files:
        raise SystemExit(f"No files matched {args.input_dir / args.file_pattern}")

    frames = []
    format_counts = {"txt": 0, "parquet": 0}
    for path in files:
        match = gene_regex.search(path.name)
        if match is None:
            print(f"[SKIP] filename does not match --gene-regex: {path.name}", file=sys.stderr)
            continue
        gene = match.group("gene")
        input_format = args.input_format
        if input_format == "auto":
            input_format = "parquet" if path.suffix.lower() in {".parquet", ".pq"} else "txt"

        if input_format == "parquet":
            source = pl.scan_parquet(path)
        else:
            source = pl.scan_csv(
                path,
                separator="\t",
                null_values=["NA", "NaN", ""],
                schema_overrides={
                    "MarkerID": pl.String,
                    "AF_Allele2": pl.Float64,
                    "p.value": pl.Float64,
                },
            )

        frame = (
            source
            .select(["MarkerID", "AF_Allele2", "p.value"])
            .filter(
                (pl.col("AF_Allele2") > args.maf_min)
                & (pl.col("AF_Allele2") < args.maf_max)
            )
            .with_columns(
                pl.lit(gene).alias("Gene"),
                pl.col("p.value").alias("pval_main"),
            )
            .select(["Gene", "MarkerID", "AF_Allele2", "pval_main"])
        )
        frames.append(frame)
        format_counts[input_format] += 1

    if not frames:
        raise SystemExit("No filenames matched --gene-regex")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    combined = pl.concat(frames, how="vertical").collect(engine="streaming")
    combined.write_csv(args.output, separator="\t")
    print(f"Loaded files: {len(frames)}")
    print(f"Input formats: txt={format_counts['txt']}, parquet={format_counts['parquet']}")
    print(f"Output rows: {combined.height}")
    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()
