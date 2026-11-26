#!/usr/bin/env python3
"""Step 2: sample, download, scan.

For each year, candidates are put in a fixed random order (seed in config.py).
Articles are processed in that order until the target number of eligible
articles (Excel supplement with a gene-symbol column) is reached. The
analysis later keeps exactly the first N eligible articles, so the sample is
a simple random sample of eligible articles.

Everything is saved to data/scan.sqlite. The run can be stopped (Ctrl-C) and
started again; finished articles are skipped. Downloaded spreadsheets are
deleted right after they are scanned.
"""
import argparse
import datetime as dt
import os
import sys

import config
import gene_scan
import store
import pipeline


def load_candidates(year, tag=""):
    """Candidate PMCIDs for a year. An article that also matched an earlier
    year's search is kept only in the earliest year (none were seen in
    checks, but this guards against double counting)."""
    p = config.CAND / ("candidates_%d%s.csv" % (year, tag))
    if not p.exists():
        raise SystemExit("missing %s. Run 01_build_candidates.py first." % p)
    ids = {int(x) for x in p.read_text().split()[1:]}
    earlier = set()
    for y in config.YEARS:
        if y >= year:
            break
        q = config.CAND / ("candidates_%d%s.csv" % (y, tag))
        if q.exists():
            earlier.update(int(x) for x in q.read_text().split()[1:])
    dup = ids & earlier
    if dup:
        print("%d: %d candidates also matched an earlier year; kept there" % (year, len(dup)))
    return sorted(ids - earlier)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--years", type=int, nargs="+", default=config.YEARS)
    ap.add_argument("--target", type=int, default=config.TARGET_ELIGIBLE_PER_YEAR,
                    help="eligible articles wanted per year (default %(default)s)")
    ap.add_argument("--max-candidates", type=int, default=None,
                    help="stop after this many candidates per year (pilots)")
    ap.add_argument("--threads", type=int, default=8, help="download threads inside each worker process")
    ap.add_argument("--procs", type=int, default=min(os.cpu_count() or 1, 6),
                    help="worker processes for parsing spreadsheets (default: CPU cores, at most 6)")
    ap.add_argument("--db", default=str(config.DB_PATH))
    ap.add_argument("--tag", default="", help="candidate file suffix (pilots)")
    ap.add_argument("--original-lists-only", action="store_true",
                    help="gene detection with the 2015 Ensembl lists only (sensitivity run; use a separate --db)")
    ap.add_argument("--source-dir", default=None, help="read from a local folder laid out like the bucket (tests)")
    args = ap.parse_args()

    hgnc_path = config.REF / "hgnc_complete_set.txt"
    if not hgnc_path.exists():
        raise SystemExit("reference data missing. Run 00_setup_reference.py first.")
    hgnc = gene_scan.load_hgnc(hgnc_path)
    sets = gene_scan.load_gene_sets(config.REF, hgnc, original_only=args.original_lists_only)
    spec = {"source_dir": args.source_dir, "original_only": args.original_lists_only,
            "threads": args.threads}

    con = store.connect(args.db)
    store.set_info(con, "seed", config.SEED)
    store.set_info(con, "gene_lists", "original_only" if args.original_lists_only else "original+HGNC")
    store.set_info(con, "gene_set_sizes", {k: len(v) for k, v in sets.items()})
    store.set_info(con, "last_run_utc", dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"))
    store.set_info(con, "target_per_year", args.target)
    store.set_info(con, "max_file_mb", config.MAX_FILE_MB)

    for year in args.years:
        ids = load_candidates(year, args.tag)
        n = pipeline.run_year(con, year, ids, spec, args.target,
                              max_candidates=args.max_candidates, procs=args.procs)
        print("%d: finished with %d eligible articles" % (year, n), flush=True)
    con.close()


if __name__ == "__main__":
    sys.exit(main())
