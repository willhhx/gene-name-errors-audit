#!/usr/bin/env python3
"""Step 0: download reference data and rebuild the published 2014-2020 rates.

Downloads (all public):
  * species gene lists used by the 2021 scanner (GitHub markziemann/GeneNameErrors2020)
  * HGNC complete set (current human symbols, previous symbols, aliases)
  * the 2021 study's raw result files, from which the yearly rates for
    2014-2020 are recomputed (so the calibration does not depend on typing
    numbers by hand)

Writes data/reference/ and data/reference/published_rates_2014_2020.csv.
A source_manifest.json with the download date and checksums is saved too.
"""
import argparse
import csv
import datetime as dt
import hashlib
import json
import sys
from collections import defaultdict

import config
from pmc_access import http_get


def fetch(url, dest, force=False):
    if dest.exists() and dest.stat().st_size > 0 and not force:
        return dest
    r = http_get(url, timeout=120)
    if r.status_code != 200:
        raise SystemExit("could not download %s (HTTP %d)" % (url, r.status_code))
    dest.write_bytes(r.content)
    return dest


def md5(path):
    return hashlib.md5(path.read_bytes()).hexdigest()


def published_rates(ref):
    """Yearly article-level rates from the 2021 study's raw files."""
    true_pos = set()
    year_of_tp = {}
    for line in (ref / "orig_true_positive.txt").read_text().splitlines():
        f = line.split("\t")
        if len(f) >= 4:
            true_pos.add(f[0])
            year_of_tp[f[0]] = f[3]
    rows = []
    for year in config.ORIG_YEARS:
        have_list, suspected = set(), set()
        path = ref / ("orig_results%d.txt" % year)
        for line in path.read_text(errors="replace").splitlines():
            f = line.split()
            if len(f) >= 3:                      # a gene column was found
                have_list.add(f[0])
                if len(f) >= 4:                  # ... and it held suspect cells
                    suspected.add(f[0])
        confirmed = {p for p in suspected if p in true_pos}
        rows.append({
            "year": year,
            "articles_with_gene_list": len(have_list),
            "articles_suspected": len(suspected),
            "articles_confirmed": len(confirmed),
        })
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--force", action="store_true", help="download again")
    args = ap.parse_args()

    ref = config.REF
    ref.mkdir(parents=True, exist_ok=True)
    manifest = {"downloaded_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
                "files": {}}

    for sp in config.ORIG_SPECIES_LISTS:
        p = fetch("%s/genelists/%s_genes" % (config.ORIG_RAW, sp), ref / ("%s_genes" % sp), args.force)
        manifest["files"][p.name] = {"md5": md5(p), "lines": len(p.read_text().splitlines())}
    p = fetch(config.HGNC_URL, ref / "hgnc_complete_set.txt", args.force)
    manifest["files"][p.name] = {"md5": md5(p), "url": config.HGNC_URL}
    for year in config.ORIG_YEARS:
        p = fetch("%s/results/genom/results%d.txt" % (config.ORIG_RAW, year),
                  ref / ("orig_results%d.txt" % year), args.force)
        manifest["files"][p.name] = {"md5": md5(p)}
    p = fetch("%s/results/genom/true_positive.txt" % config.ORIG_RAW,
              ref / "orig_true_positive.txt", args.force)
    manifest["files"][p.name] = {"md5": md5(p)}

    rows = published_rates(ref)
    out = ref / "published_rates_2014_2020.csv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]) + ["rate_confirmed"])
        w.writeheader()
        for r in rows:
            r["rate_confirmed"] = round(r["articles_confirmed"] / r["articles_with_gene_list"], 4)
            w.writerow(r)
    tot_l = sum(r["articles_with_gene_list"] for r in rows)
    tot_c = sum(r["articles_confirmed"] for r in rows)
    print("Published series rebuilt from the 2021 study's raw files:")
    for r in rows:
        print("  %d: %5d articles with gene lists, %4d confirmed errors (%.1f%%)" % (
            r["year"], r["articles_with_gene_list"], r["articles_confirmed"], 100 * r["rate_confirmed"]))
    print("  total: %d with gene lists, %d confirmed (%.1f%%)" % (tot_l, tot_c, 100 * tot_c / tot_l))
    print("  the paper reports 3,436 / 11,117 (30.9%). Check that these totals agree.")
    manifest["published_totals"] = {"with_gene_list": tot_l, "confirmed": tot_c}
    (ref / "source_manifest.json").write_text(json.dumps(manifest, indent=2))
    print("saved", out)


if __name__ == "__main__":
    sys.exit(main())
