#!/usr/bin/env python3
"""Step 1: build the candidate list of PMC articles for each year.

Query (same idea as the 2021 study): genom* in the abstract, open access,
publication year = Y. NCBI returns at most ~10,000 ids per search, so the
PMC id range is split until each piece fits.

Writes data/candidates/candidates_<year>.csv (one pmcid per line) and
data/candidates/summary.csv. Safe to re-run; existing years are skipped
unless --force is given.
"""
import argparse
import csv
import datetime as dt
import sys

import config
import pmc_access


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--years", type=int, nargs="+", default=config.YEARS)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--uid-range", type=int, nargs=2, metavar=("LO", "HI"),
                    help="only search this PMC id range (for quick pilots)")
    ap.add_argument("--tag", default="", help="suffix for output file names (pilots)")
    args = ap.parse_args()

    config.CAND.mkdir(parents=True, exist_ok=True)
    summary = config.CAND / ("summary%s.csv" % args.tag)
    rows = {}
    if summary.exists():
        with open(summary) as fh:
            rows = {int(r["year"]): r for r in csv.DictReader(fh)}

    for year in args.years:
        out = config.CAND / ("candidates_%d%s.csv" % (year, args.tag))
        if out.exists() and not args.force:
            print("%d: already have %s (use --force to rebuild)" % (year, out.name))
            continue
        lo, hi = args.uid_range or (1, 20_000_000)
        print("%d: searching PMC ..." % year, flush=True)
        ids = sorted(set(pmc_access.candidate_ids(year, lo, hi)))
        out.write_text("pmcid\n" + "\n".join(map(str, ids)) + ("\n" if ids else ""))
        rows[year] = {"year": year, "n_candidates": len(ids),
                      "query": config.CANDIDATE_QUERY,
                      "retrieved_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")}
        print("%d: %d candidate articles -> %s" % (year, len(ids), out.name), flush=True)

    with open(summary, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["year", "n_candidates", "query", "retrieved_utc"])
        w.writeheader()
        for y in sorted(rows):
            w.writerow(rows[y])


if __name__ == "__main__":
    sys.exit(main())
