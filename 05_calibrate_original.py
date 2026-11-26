#!/usr/bin/env python3
"""Step 5: head-to-head calibration on the 2021 study's own articles.

Takes a random sample of the articles that the 2021 study found to have Excel
gene lists in 2019 and 2020, rescans them with this pipeline, and compares
article by article:
  * did we also find a gene list?
  * do we flag the same articles as having errors (agreement, Cohen's kappa)?
  * is our error rate on these same articles close to theirs?

This is the strongest answer to "is your scanner comparable with the
original?", because it removes any difference in which articles were sampled.

Writes data/calibration.sqlite, output/calibration_original_articles.csv and
output/calibration_original_summary.json.
"""
import argparse
import json
import random
import sys

import pandas as pd

import config
import pipeline
import store


def original_articles(year):
    """pmcid -> (suspected, confirmed) for articles with a gene list in the 2021 study."""
    tp = set()
    for line in (config.REF / "orig_true_positive.txt").read_text().splitlines():
        f = line.split("\t")
        if f and f[0].startswith("PMC"):
            tp.add(f[0])
    has, sus = set(), set()
    for line in (config.REF / ("orig_results%d.txt" % year)).read_text(errors="replace").splitlines():
        f = line.split()
        if len(f) >= 3 and f[0].startswith("PMC"):
            has.add(f[0])
            if len(f) >= 4:
                sus.add(f[0])
    return {int(p[3:]): (int(p in sus), int(p in sus and p in tp)) for p in has}


def kappa(a, b):
    a, b = pd.Series(a), pd.Series(b)
    po = (a == b).mean()
    pe = a.mean() * b.mean() + (1 - a.mean()) * (1 - b.mean())
    return float("nan") if pe == 1 else (po - pe) / (1 - pe)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--years", type=int, nargs="+", default=[2019, 2020])
    ap.add_argument("--n", type=int, default=250, help="articles per year")
    ap.add_argument("--threads", type=int, default=12)
    ap.add_argument("--db", default=str(config.DATA / "calibration.sqlite"))
    args = ap.parse_args()

    if not (config.REF / "orig_true_positive.txt").exists():
        raise SystemExit("run 00_setup_reference.py first")
    config.EXCLUDED_ARTICLE_TYPES = set()   # the original scanned every article type
    pipeline.init_worker({"threads": args.threads})
    con = store.connect(args.db)

    rows = []
    for year in args.years:
        orig = original_articles(year)
        ids = sorted(orig)
        random.Random(config.SEED + year).shuffle(ids)
        ids = ids[:args.n]
        done = {r[0] for r in con.execute(
            "SELECT pmcid FROM articles WHERE year=? AND status!='fetch_error'", (year,))}
        todo = [(i, p, year) for i, p in enumerate(ids) if p not in done]
        print("%d: %d sampled from %d original gene-list articles, %d to scan" % (
            year, len(ids), len(orig), len(todo)), flush=True)
        for k in range(0, len(todo), 48):
            for art, files, errs in pipeline.work_chunk(todo[k:k + 48]):
                store.save_result(con, art, files, errs)
            con.commit()
        for p in ids:
            r = con.execute("SELECT status, has_gene_list, any_error_broad, any_error_strict "
                            "FROM articles WHERE pmcid=?", (p,)).fetchone()
            sus, conf = orig[p]
            rows.append({"year": year, "pmcid": p, "status": r[0] if r else "missing",
                         "ours_gene_list": r[1] if r else None, "ours_error": r[2] if r else None,
                         "ours_error_strict": r[3] if r else None,
                         "orig_suspected": sus, "orig_confirmed": conf})

    df = pd.DataFrame(rows)
    config.OUT.mkdir(parents=True, exist_ok=True)
    df.to_csv(config.OUT / "calibration_original_articles.csv", index=False)

    out = {"n_sampled": len(df), "status": df.status.value_counts().to_dict()}
    ok = df[df.status == "scanned"].copy()
    out["n_scanned"] = len(ok)
    out["gene_list_found_share"] = float(ok.ours_gene_list.mean())
    both = ok[ok.ours_gene_list == 1]
    out["n_both_found_gene_list"] = len(both)
    res = {}
    for ref in ("orig_suspected", "orig_confirmed"):
        t = pd.crosstab(both[ref], both.ours_error)
        res[ref] = {
            "orig_rate": float(both[ref].mean()), "ours_rate": float(both.ours_error.mean()),
            "agreement": float((both[ref] == both.ours_error).mean()),
            "kappa": kappa(both[ref], both.ours_error),
            "ours_flags_orig_positive": float(both[both[ref] == 1].ours_error.mean()),
            "ours_flags_orig_negative": float(both[both[ref] == 0].ours_error.mean()),
            "table": {"%s=%d,ours=%d" % (ref, i, j): int(t.loc[i, j]) if i in t.index and j in t.columns else 0
                      for i in (0, 1) for j in (0, 1)}}
    out["vs"] = res
    out["by_year"] = {int(y): {"n": len(g), "orig_confirmed_rate": float(g.orig_confirmed.mean()),
                               "orig_suspected_rate": float(g.orig_suspected.mean()),
                               "ours_rate": float(g.ours_error.mean())}
                      for y, g in both.groupby("year")}
    (config.OUT / "calibration_original_summary.json").write_text(json.dumps(out, indent=2))

    print("\nScanned %d of %d sampled articles; we found a gene list in %.1f%%" % (
        len(ok), len(df), 100 * out["gene_list_found_share"]))
    for ref, v in res.items():
        print("vs %s: theirs %.1f%%, ours %.1f%%, agreement %.1f%%, kappa %.2f "
              "(we flag %.1f%% of their positives, %.1f%% of their negatives)" % (
                  ref, 100 * v["orig_rate"], 100 * v["ours_rate"], 100 * v["agreement"], v["kappa"],
                  100 * v["ours_flags_orig_positive"], 100 * v["ours_flags_orig_negative"]))


if __name__ == "__main__":
    sys.exit(main())
