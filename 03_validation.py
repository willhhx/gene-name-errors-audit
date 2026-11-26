#!/usr/bin/env python3
"""Step 3: manual validation of the scanner (plan: 200 flagged + 200 unflagged files).

  python 03_validation.py make   -> output/validation_sheet.csv (random files, labels blank)
  python 03_validation.py fetch  -> downloads those files to data/validation_files/ and
                                    re-scans them (checks the scan is repeatable)
  (open the files, fill the "label" column by hand)
  python 03_validation.py score  -> positive predictive value and miss rate

How to label
  flagged files   : open the file, go to the listed sheet/row/column. label = 1 if the
                    cell really is a damaged gene name (for example 2-Sep where the gene
                    should be SEPTIN2/SEPT2), 0 if it is a false alarm (a real date or a
                    real number that belongs there).
  unflagged files : open the file, look at the gene-symbol column(s) listed in
                    "gene_columns" (sort the column, look for dates and 5-digit numbers).
                    label = 1 if you find a damaged gene name the scanner missed, else 0.
  Write the rule or reason in "rule_used" (for example "date in gene column").
"""
import argparse
import csv
import json
import random
import sys
import tempfile
import zipfile
from pathlib import Path

import config
import gene_scan
import pmc_access
import store

SHEET = config.OUT / "validation_sheet.csv"
FILES = config.DATA / "validation_files"
COLS = ["id", "group", "pmcid", "year", "file", "gene_columns", "example_cells",
        "pmc_page", "local_copy", "label", "rule_used", "notes"]


def col_letter(n):  # 1 -> A
    s = ""
    while n:
        n, r = divmod(n - 1, 26)
        s = chr(65 + r) + s
    return s


def wilson(k, n, z=1.96):
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * ((p * (1 - p) / n + z * z / (4 * n * n)) ** 0.5) / d
    return (c - h, c + h)


def make(args):
    con = store.connect(args.db)
    keep = []
    for y in config.YEARS:
        keep += store.eligible_pmcids(con, y, args.target)
    if not keep:
        raise SystemExit("no eligible articles in the database yet; run step 2 first")
    marks = ",".join("?" * len(keep))
    flagged = con.execute(
        "SELECT f.pmcid, f.name, a.year, f.n_errors FROM files f JOIN articles a USING(pmcid) "
        "WHERE f.pmcid IN (%s) AND f.parse_status='ok' AND f.n_gene_columns>0 AND f.n_errors>0" % marks, keep).fetchall()
    clean = con.execute(
        "SELECT f.pmcid, f.name, a.year, f.n_errors FROM files f JOIN articles a USING(pmcid) "
        "WHERE f.pmcid IN (%s) AND f.parse_status='ok' AND f.n_gene_columns>0 AND f.n_errors=0" % marks, keep).fetchall()
    rng = random.Random(config.SEED)
    rng.shuffle(flagged)
    rng.shuffle(clean)
    rows = []
    for group, pool in (("flagged", flagged), ("unflagged", clean)):
        for pmcid, name, year, n_err in pool[:args.n]:
            ex = []
            if group == "flagged":
                for e in con.execute(
                        "SELECT sheet,row,col,raw,form FROM errors WHERE pmcid=? AND file=? LIMIT 6",
                        (pmcid, name)):
                    ex.append("%s!%s%d=%s (%s)" % (e[0], col_letter(e[2]), e[1], e[3], e[4]))
            rows.append({"id": len(rows) + 1, "group": group, "pmcid": pmcid, "year": year,
                         "file": name, "example_cells": " | ".join(ex),
                         "pmc_page": "https://pmc.ncbi.nlm.nih.gov/articles/PMC%d/" % pmcid})
    config.OUT.mkdir(parents=True, exist_ok=True)
    with open(SHEET, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print("wrote %s: %d flagged and %d unflagged files (available: %d and %d)" % (
        SHEET, min(args.n, len(flagged)), min(args.n, len(clean)), len(flagged), len(clean)))


def fetch(args):
    if not SHEET.exists():
        raise SystemExit("run 'make' first")
    hgnc = gene_scan.load_hgnc(config.REF / "hgnc_complete_set.txt")
    sets = gene_scan.load_gene_sets(config.REF, hgnc)
    src = pmc_access.HttpSource()
    FILES.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader(open(SHEET)))
    con = store.connect(args.db)
    n_mismatch = 0
    for r in rows:
        pmcid = int(r["pmcid"])
        meta = src.metadata(pmcid)
        version = int(meta["version"])
        name = r["file"]
        zip_name, _, member = name.partition("::")
        base = pmc_access.media_name(next(u for u in meta["media_urls"]
                                          if pmc_access.media_name(u) == zip_name))
        safe = "%s_%s" % (r["id"].zfill(3), (member or base).rsplit("/", 1)[-1].replace(" ", "_"))
        dest = FILES / safe
        if not dest.exists():
            with tempfile.TemporaryDirectory() as td:
                tmp = Path(td) / "dl"
                src.download(pmcid, version, base, tmp, None, 400 * 1024 * 1024)
                if member:
                    with zipfile.ZipFile(tmp) as zf:
                        dest.write_bytes(zf.read(member))
                else:
                    dest.write_bytes(tmp.read_bytes())
        res = gene_scan.scan_workbook(dest, sets, hgnc)
        r["local_copy"] = str(dest)
        r["gene_columns"] = "; ".join("%s!%s (%s)" % (g.sheet, col_letter(g.col + 1), g.species)
                                      for g in res.gene_columns[:8])
        stored = con.execute("SELECT n_errors FROM files WHERE pmcid=? AND name=?", (pmcid, name)).fetchone()
        if stored is not None and stored[0] != len(res.errors):
            n_mismatch += 1
            r["notes"] = "rescan found %d errors, database has %d" % (len(res.errors), stored[0])
    with open(SHEET, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS)
        w.writeheader()
        w.writerows(rows)
    print("downloaded %d files to %s; rescan disagreed with the database for %d" % (len(rows), FILES, n_mismatch))
    print("Now open each file, fill the 'label' column in %s, then run: python 03_validation.py score" % SHEET)


def score(args):
    rows = list(csv.DictReader(open(SHEET)))
    out = {}
    for group in ("flagged", "unflagged"):
        lab = [r for r in rows if r["group"] == group and r["label"].strip() in ("0", "1")]
        todo = [r for r in rows if r["group"] == group and r["label"].strip() not in ("0", "1")]
        out[group] = {"labelled": len(lab), "unlabelled": len(todo),
                      "ones": sum(r["label"].strip() == "1" for r in lab)}
    f, u = out["flagged"], out["unflagged"]
    if f["labelled"]:
        ppv = f["ones"] / f["labelled"]
        lo, hi = wilson(f["ones"], f["labelled"])
        out["ppv"] = {"value": ppv, "ci95": [lo, hi], "n": f["labelled"]}
        out["decision"] = "GO" if ppv >= 0.90 else ("BORDERLINE: tighten rules and recheck" if ppv >= 0.80 else "NO-GO")
        print("Positive predictive value: %.1f%% (%.1f to %.1f), n=%d flagged files -> %s" % (
            100 * ppv, 100 * lo, 100 * hi, f["labelled"], out["decision"]))
    if u["labelled"]:
        miss = u["ones"] / u["labelled"]
        lo, hi = wilson(u["ones"], u["labelled"])
        out["miss_rate"] = {"value": miss, "ci95": [lo, hi], "n": u["labelled"]}
        print("Miss rate among unflagged gene-list files: %.1f%% (%.1f to %.1f), n=%d" % (
            100 * miss, 100 * lo, 100 * hi, u["labelled"]))
    print("still unlabelled: %d flagged, %d unflagged" % (f["unlabelled"], u["unlabelled"]))
    (config.OUT / "validation_scores.json").write_text(json.dumps(out, indent=2))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("action", choices=["make", "fetch", "score"])
    ap.add_argument("--n", type=int, default=200, help="files per group")
    ap.add_argument("--target", type=int, default=config.TARGET_ELIGIBLE_PER_YEAR)
    ap.add_argument("--db", default=str(config.DB_PATH))
    args = ap.parse_args()
    {"make": make, "fetch": fetch, "score": score}[args.action](args)


if __name__ == "__main__":
    sys.exit(main())
