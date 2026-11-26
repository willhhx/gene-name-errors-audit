#!/usr/bin/env python3
"""Step 4: yearly rates, trend tests, before/after comparison, sensitivity
analyses, calibration against the published 2014-2020 series, figures, tables.

Reads data/scan.sqlite. Writes everything to output/.
Keeps exactly the first N eligible articles per year (N = --target) in the
seeded random order, so the sample is a simple random sample.
"""
import argparse
import json
import re
import sys

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from statsmodels.stats.proportion import (confint_proportions_2indep,
                                          proportion_confint,
                                          test_proportions_2indep)

import config
import store

CATS = ["renamed_human", "non_renamed_human", "mouse", "other"]
CAT_LABEL = {"renamed_human": "Renamed human genes (SEPT, MARCH, DEC1 ...)",
             "non_renamed_human": "Other human symbols",
             "mouse": "Mouse symbols", "other": "Other (numbers, unresolved)"}
# validated categorical slots 1-4 (blue, orange, aqua, yellow), light surface
COLORS = {"renamed_human": "#2a78d6", "non_renamed_human": "#eb6834",
          "mouse": "#1baf7a", "other": "#eda100"}
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"


def wilson(k, n):
    if n == 0:
        return (np.nan, np.nan)
    lo, hi = proportion_confint(k, n, alpha=0.05, method="wilson")
    return (lo, hi)


def norm_journal(s):
    import html
    s = html.unescape(s or "").lower().replace("&", "and")
    s = re.sub(r"^the\s+", "", s)
    s = re.sub(r"\s*[:(].*$", "", s)   # "Science (New York, N.Y.)" -> "science"
    return re.sub(r"[^a-z0-9]", "", s)


# ---- loading ---------------------------------------------------------------
def load_sample(con, target):
    frames = []
    for y in config.YEARS:
        ids = store.eligible_pmcids(con, y, target)
        if ids:
            frames.append(pd.read_sql_query(
                "SELECT * FROM articles WHERE pmcid IN (%s)" % ",".join(map(str, ids)), con))
    if not frames:
        raise SystemExit("no eligible articles yet. Run 02_scan_articles.py first.")
    df = pd.concat(frames, ignore_index=True)
    df["other_types"] = ((df.any_non_renamed_human + df.any_mouse + df.any_other) > 0).astype(int)
    return df


def yearly_table(df, col="any_error_broad"):
    rows = []
    for y, g in df.groupby("year"):
        n, k = len(g), int(g[col].sum())
        lo, hi = wilson(k, n)
        rows.append({"year": y, "n_eligible": n, "n_error": k,
                     "rate": k / n, "ci_low": lo, "ci_high": hi})
    return pd.DataFrame(rows)


def before_after(df, col, before=config.BEFORE_YEARS, after=config.AFTER_YEARS):
    b = df[df.year.isin(before)][col]
    a = df[df.year.isin(after)][col]
    if len(a) == 0 or len(b) == 0:
        return None
    kb, nb, ka, na = int(b.sum()), len(b), int(a.sum()), len(a)
    d = ka / na - kb / nb
    out = {"before_k": kb, "before_n": nb, "before_rate": kb / nb,
           "after_k": ka, "after_n": na, "after_rate": ka / na, "diff_points": 100 * d}
    try:
        lo, hi = confint_proportions_2indep(ka, na, kb, nb, method="newcomb", compare="diff")
        out.update(diff_ci_low=100 * lo, diff_ci_high=100 * hi)
        out["p_value"] = float(test_proportions_2indep(ka, na, kb, nb, method="score",
                                                      compare="diff").pvalue)
    except Exception:
        out.update(diff_ci_low=np.nan, diff_ci_high=np.nan, p_value=np.nan)
    return out


def trend_tests(df, col="any_error_broad"):
    """Logistic regression of error on year; check for a non-linear pattern."""
    d = df[["year", col]].copy()
    years = sorted(d.year.unique())
    out = {"years": [int(y) for y in years]}
    if len(years) < 3 or d[col].nunique() < 2:
        out["note"] = "not enough variation to fit a trend"
        return out
    x = (d.year - years[0]).astype(float)
    lin = sm.GLM(d[col], sm.add_constant(x), family=sm.families.Binomial()).fit()
    b, se = lin.params.iloc[1], lin.bse.iloc[1]
    out.update(odds_ratio_per_year=float(np.exp(b)),
               or_ci_low=float(np.exp(b - 1.96 * se)), or_ci_high=float(np.exp(b + 1.96 * se)),
               p_trend=float(lin.pvalues.iloc[1]))
    dummies = pd.get_dummies(d.year, drop_first=True).astype(float)
    cat = sm.GLM(d[col], sm.add_constant(dummies), family=sm.families.Binomial()).fit()
    lr = 2 * (cat.llf - lin.llf)
    df_lr = len(years) - 2
    out.update(nonlinearity_lr=float(lr), nonlinearity_df=df_lr,
               nonlinearity_p=float(stats.chi2.sf(lr, df_lr)))
    return out


def bootstrap_diff_of_changes(df, col_a="any_renamed_human", col_b="other_types",
                              reps=4000):
    """Is the change in renamed-gene errors different from the change in the
    other error types? Resamples articles (the two columns are paired)."""
    b = df[df.year.isin(config.BEFORE_YEARS)]
    a = df[df.year.isin(config.AFTER_YEARS)]
    if len(a) == 0 or len(b) == 0:
        return None
    rng = np.random.default_rng(config.SEED)
    stat = lambda B, A: (A[col_a].mean() - B[col_a].mean()) - (A[col_b].mean() - B[col_b].mean())
    obs = stat(b, a)
    vals = np.empty(reps)
    ba, bb = b[[col_a, col_b]].to_numpy(), a[[col_a, col_b]].to_numpy()
    for i in range(reps):
        sb = ba[rng.integers(0, len(ba), len(ba))].mean(axis=0)
        sa = bb[rng.integers(0, len(bb), len(bb))].mean(axis=0)
        vals[i] = (sa[0] - sb[0]) - (sa[1] - sb[1])
    lo, hi = np.percentile(vals, [2.5, 97.5])
    return {"diff_of_changes_points": 100 * obs, "ci_low": 100 * lo, "ci_high": 100 * hi}


def publisher_adjusted(df, col="any_error_broad", top=12):
    """Before/after odds ratio adjusted for publisher mix (the mix of journals
    in PMC changed a lot between 2019 and 2025). Logistic regression of error
    on period plus publisher (the `top` largest, the rest pooled)."""
    d = df[df.year.isin(config.BEFORE_YEARS + config.AFTER_YEARS)].copy()
    if len(d) == 0 or d[col].nunique() < 2:
        return None
    d["after"] = d.year.isin(config.AFTER_YEARS).astype(float)
    pub = d.publisher.fillna("").str.strip().str.lower().replace("", "unknown")
    keep = pub.value_counts().index[:top]
    d["pub"] = pub.where(pub.isin(keep), "other")
    X = pd.concat([d[["after"]], pd.get_dummies(d.pub, prefix="pub", drop_first=True).astype(float)], axis=1)
    try:
        crude = sm.GLM(d[col], sm.add_constant(d[["after"]]), family=sm.families.Binomial()).fit()
        adj = sm.GLM(d[col], sm.add_constant(X), family=sm.families.Binomial()).fit()
    except Exception as exc:
        return {"note": "model failed: %s" % exc}
    def orci(m):
        b, se = m.params["after"], m.bse["after"]
        return float(np.exp(b)), float(np.exp(b - 1.96 * se)), float(np.exp(b + 1.96 * se)), float(m.pvalues["after"])
    c, a = orci(crude), orci(adj)
    mix = (d.groupby(["after", "pub"]).size() / d.groupby("after").size()).unstack(0).fillna(0)
    mix.columns = ["share_before", "share_after"]
    return {"crude_or": c[0], "crude_ci": [c[1], c[2]], "crude_p": c[3],
            "adjusted_or": a[0], "adjusted_ci": [a[1], a[2]], "adjusted_p": a[3],
            "publishers_in_model": int(d.pub.nunique()),
            "publisher_mix": mix.round(4).reset_index().to_dict("records")}


def direct_files_only(con, df):
    """Sensitivity: count only Excel files attached directly to the article
    (not inside zip archives), closer to how most readers meet the files."""
    ids = ",".join(map(str, df.pmcid)) or "0"
    f = pd.read_sql_query("SELECT pmcid, source, n_gene_columns FROM files WHERE parse_status='ok' "
                          "AND pmcid IN (%s)" % ids, con)
    e = pd.read_sql_query("SELECT DISTINCT pmcid, file FROM errors WHERE pmcid IN (%s)" % ids, con)
    has = set(f[(f.source == "direct") & (f.n_gene_columns > 0)].pmcid)
    err = set(e[~e.file.str.contains("::", regex=False)].pmcid)
    d = df[df.pmcid.isin(has)].copy()
    d["direct_error"] = d.pmcid.isin(err).astype(int)
    return d


# ---- exploratory -----------------------------------------------------------
def exploratory(con, target, df):
    out = {}
    rows, fmt_rows = [], []
    for y in config.YEARS:
        cut = store.cutoff_order_idx(con, y, target)
        if cut < 0:
            continue
        a = pd.read_sql_query(
            "SELECT * FROM articles WHERE year=? AND order_idx<=? AND status IN "
            "('no_excel','scanned','excluded_type')", con, params=(y, cut))
        n = len(a)
        has_x = int(((a.n_excel_files > 0)).sum())
        lo, hi = wilson(has_x, n)
        rows.append({"year": y, "candidates_checked": n, "with_excel_supplement": has_x,
                     "share": has_x / n if n else np.nan, "ci_low": lo, "ci_high": hi})
        counts = {}
        for s in a.supp_formats.fillna(""):
            for kv in filter(None, s.split(";")):
                k, v = kv.split("=")
                counts[k] = counts.get(k, 0) + 1
        for k, v in sorted(counts.items(), key=lambda kv: -kv[1]):
            fmt_rows.append({"year": y, "format": k, "articles_with_format": v,
                             "share_of_candidates": v / n})
    out["excel_share"] = pd.DataFrame(rows)
    out["format_mix"] = pd.DataFrame(fmt_rows)

    ids = ",".join(map(str, df.pmcid)) or "0"
    files = pd.read_sql_query(
        "SELECT a.year, f.parse_status, f.ext FROM files f JOIN articles a USING(pmcid) "
        "WHERE f.pmcid IN (%s)" % ids, con)
    files["problem"] = (files.parse_status != "ok").astype(int)
    out["file_problems"] = files.groupby("year").agg(
        files=("problem", "size"), not_scanned=("problem", "sum")).reset_index()
    out["file_problems"]["share_not_scanned"] = (
        out["file_problems"].not_scanned / out["file_problems"].files)
    reasons = files[files.problem == 1].groupby(["year", "parse_status"]).size().reset_index(name="files")
    out["file_problem_reasons"] = reasons

    err = pd.read_sql_query(
        "SELECT a.year, e.pmcid, e.form, e.category, e.candidate, e.raw FROM errors e "
        "JOIN articles a USING(pmcid) WHERE e.pmcid IN (%s)" % ids, con)
    if len(err):
        err["symbol"] = np.where(err.candidate.fillna("") != "", err.candidate,
                                 err.form + ":" + err.raw.fillna(""))
        top = (err.drop_duplicates(["year", "pmcid", "symbol"])
               .groupby(["year", "symbol"]).size().reset_index(name="articles")
               .sort_values(["year", "articles"], ascending=[True, False]))
        out["top_symbols"] = top.groupby("year").head(15)
    else:
        out["top_symbols"] = pd.DataFrame(columns=["year", "symbol", "articles"])
    out["dup_check"] = int(df.pmcid.duplicated().sum())
    return out


# ---- figures ---------------------------------------------------------------
def style_axes(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=9)
    ax.grid(axis="y", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def fig1(ours, pub, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7.2, 4.2), dpi=200)
    style_axes(ax)
    if pub is not None and len(pub):
        lo, hi = zip(*[wilson(k, n) for k, n in zip(pub.articles_confirmed, pub.articles_with_gene_list)])
        r = pub.articles_confirmed / pub.articles_with_gene_list
        ax.errorbar(pub.year, 100 * r, yerr=[100 * (r - np.array(lo)), 100 * (np.array(hi) - r)],
                    fmt="o-", color="#8a8984", ecolor="#8a8984", lw=1.5, ms=6, capsize=3,
                    label="2014-2020, published (Abeysooriya et al. 2021)")
    ax.errorbar(ours.year, 100 * ours.rate,
                yerr=[100 * (ours.rate - ours.ci_low), 100 * (ours.ci_high - ours.rate)],
                fmt="o-", color=COLORS["renamed_human"], ecolor=COLORS["renamed_human"],
                lw=2, ms=7, capsize=3, mfc=COLORS["renamed_human"], mec="#fcfcfb",
                label="2019-2025, this study")
    ax.axvline(config.HGNC_RENAMING_YEAR + 0.5, color=INK2, lw=1, ls=(0, (4, 3)))
    ax.axvline(config.EXCEL_UPDATE, color=INK2, lw=1, ls=(0, (4, 3)))
    top = max(ax.get_ylim()[1], 52)
    ax.set_ylim(0, top)
    ax.text(config.HGNC_RENAMING_YEAR + 0.58, 1.0, "HGNC renames\nerror-prone genes (2020)", ha="left",
            fontsize=8, color=INK2, va="bottom")
    ax.text(config.EXCEL_UPDATE + 0.08, 1.0, "Excel opt-out\nsetting (Oct 2023)", ha="left", fontsize=8,
            color=INK2, va="bottom")
    ax.set_xlabel("Publication year", color=INK2)
    ax.set_ylabel("Articles with gene-name errors (%)", color=INK2)
    ax.set_xticks(range(2014, 2026))
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    ax.set_title("Share of articles with Excel gene lists that contain damaged gene names",
                 fontsize=10, color=INK, loc="left")
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"))
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def fig2(df, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    years = sorted(df.year.unique())
    prio = {c: np.zeros(len(years)) for c in CATS}
    for i, y in enumerate(years):
        g = df[df.year == y]
        left = pd.Series(True, index=g.index)
        for c, col in zip(CATS, ["any_renamed_human", "any_non_renamed_human", "any_mouse", "any_other"]):
            hit = left & (g[col] == 1)
            prio[c][i] = 100 * hit.sum() / len(g)
            left &= ~hit
    fig, ax = plt.subplots(figsize=(7.2, 4.2), dpi=200)
    style_axes(ax)
    bottom = np.zeros(len(years))
    for c in CATS:
        ax.bar(years, prio[c], bottom=bottom, color=COLORS[c], width=0.7,
               edgecolor="#fcfcfb", linewidth=1.5, label=CAT_LABEL[c])
        bottom += prio[c]
    ax.set_xlabel("Publication year", color=INK2)
    ax.set_ylabel("Articles with errors (% of eligible)", color=INK2)
    ax.set_xticks(years)
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    ax.set_title("Error types by year (each article counted once, in the order listed)",
                 fontsize=10, color=INK, loc="left")
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"))
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)
    return pd.DataFrame(prio, index=years).rename_axis("year").reset_index()


# ---- main ------------------------------------------------------------------
def fmt_ba(r):
    if r is None:
        return "n/a"
    return "%.1f%% -> %.1f%% (change %+.1f points, 95%% CI %.1f to %.1f, p=%.3g)" % (
        100 * r["before_rate"], 100 * r["after_rate"], r["diff_points"],
        r["diff_ci_low"], r["diff_ci_high"], r["p_value"])


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--db", default=str(config.DB_PATH))
    ap.add_argument("--db-orig", default=str(config.DATA / "scan_original_lists.sqlite"),
                    help="optional second run made with --original-lists-only")
    ap.add_argument("--target", type=int, default=config.TARGET_ELIGIBLE_PER_YEAR)
    ap.add_argument("--out", default=str(config.OUT))
    args = ap.parse_args()
    out = __import__("pathlib").Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    con = store.connect(args.db)
    df = load_sample(con, args.target)
    summary = {"target_per_year": args.target}

    # primary table
    broad = yearly_table(df, "any_error_broad")
    strict = yearly_table(df, "any_error_strict")
    t1 = broad.merge(strict[["year", "n_error", "rate", "ci_low", "ci_high"]], on="year",
                     suffixes=("", "_strict"))
    for c in CATS:
        col = {"renamed_human": "any_renamed_human", "non_renamed_human": "any_non_renamed_human",
               "mouse": "any_mouse", "other": "any_other"}[c]
        t1["n_" + c] = df.groupby("year")[col].sum().reindex(t1.year).to_numpy()
    t1.to_csv(out / "table1_yearly.csv", index=False)

    # published series and calibration
    pub_path = config.REF / "published_rates_2014_2020.csv"
    pub = pd.read_csv(pub_path) if pub_path.exists() else None
    cal = {}
    if pub is not None:
        for y in config.BEFORE_YEARS:
            mine = broad[broad.year == y]
            p = pub[pub.year == y]
            if len(mine) and len(p):
                pr = (p.articles_confirmed / p.articles_with_gene_list).iloc[0]
                ps = (p.articles_suspected / p.articles_with_gene_list).iloc[0]
                cal[y] = {"ours": float(mine.rate.iloc[0]), "published_confirmed": float(pr),
                          "published_suspected": float(ps),
                          "gap_points_vs_confirmed": 100 * float(mine.rate.iloc[0] - pr),
                          "gap_points_vs_suspected": 100 * float(mine.rate.iloc[0] - ps)}
        if cal:
            gap = max(abs(v["gap_points_vs_confirmed"]) for v in cal.values())
            gap_s = max(abs(v["gap_points_vs_suspected"]) for v in cal.values())
            best = min(gap, gap_s)
            cal["decision"] = ("GO" if best <= 5 else "BORDERLINE: find the cause, report it"
                               if best <= 10 else "NO-GO unless explained")
            cal["largest_gap_points"] = best
    # same comparison restricted to our sampled articles that the 2021 study also
    # counted as having a gene list (removes differences in coverage)
    for y in config.BEFORE_YEARS:
        f = config.REF / ("orig_results%d.txt" % y)
        if y in cal and f.exists():
            has = set()
            for line in f.read_text(errors="replace").splitlines():
                parts = line.split()
                if len(parts) >= 3 and parts[0].startswith("PMC"):
                    has.add(int(parts[0][3:]))
            g = df[(df.year == y) & df.pmcid.isin(has)]
            cal[y]["overlap_n"] = int(len(g))
            cal[y]["overlap_share_of_our_sample"] = float(len(g) / max((df.year == y).sum(), 1))
            cal[y]["ours_rate_on_overlap"] = float(g.any_error_broad.mean()) if len(g) else None
    summary["calibration"] = cal

    # figures
    fig1(broad, pub, out / "fig1_yearly_rate")
    prio = fig2(df, out / "fig2_error_types")
    prio.to_csv(out / "fig2_data.csv", index=False)

    # primary tests
    summary["trend_broad"] = trend_tests(df, "any_error_broad")
    summary["trend_strict"] = trend_tests(df, "any_error_strict")
    ba_rows = []
    defs = [("any error (broad)", df, "any_error_broad"),
            ("dates only (strict)", df, "any_error_strict"),
            ("renamed human genes", df, "any_renamed_human"),
            ("other types merged", df, "other_types"),
            ("other human symbols", df, "any_non_renamed_human"),
            ("mouse symbols", df, "any_mouse"),
            ("other forms", df, "any_other")]
    for label, d, col in defs:
        r = before_after(d, col)
        if r:
            r["comparison"] = label
            r["small_count_note"] = "fewer than 20 cases; read with care" if (r["before_k"] + r["after_k"]) < 20 else ""
            ba_rows.append(r)
    # sensitivity
    sens = []
    def add_sens(name, d, col="any_error_broad", before=config.BEFORE_YEARS, after=config.AFTER_YEARS):
        if len(d) == 0:
            return
        r = before_after(d, col, before, after)
        if r:
            r["comparison"] = name
            sens.append(r)
    add_sens("all eligible (primary)", df)
    add_sens("human-only files", df[df.human_only == 1])
    j16 = {norm_journal(j) for j in config.JOURNALS_2016}
    add_sens("journals of the 2016 study (all 18)", df[df.journal.map(norm_journal).isin(j16)])
    add_sens("dates only (strict)", df, "any_error_strict")
    add_sens("2025 only vs 2019-2020", df[df.year.isin([2019, 2020, 2025])], before=config.BEFORE_YEARS, after=(2025,))
    dd = direct_files_only(con, df)
    add_sens("Excel files attached directly (zips ignored)", dd, "direct_error")
    adj = publisher_adjusted(df)
    summary["publisher_adjusted"] = adj
    if adj and "adjusted_or" in adj:
        pd.DataFrame(adj["publisher_mix"]).to_csv(out / "table_publisher_mix.csv", index=False)
    pd.DataFrame(ba_rows).to_csv(out / "table_before_after.csv", index=False)
    pd.DataFrame(sens).to_csv(out / "table_sensitivity.csv", index=False)
    summary["diff_of_changes_renamed_vs_other"] = bootstrap_diff_of_changes(df)

    # original-list rerun, if present
    try:
        import os
        if os.path.exists(args.db_orig):
            c2 = store.connect(args.db_orig)
            d2 = load_sample(c2, args.target)
            yearly_table(d2).to_csv(out / "table_sensitivity_original_lists_yearly.csv", index=False)
            r = before_after(d2, "any_error_broad")
            if r:
                r["comparison"] = "original 2015 gene lists only"
                sens.append(r)
                pd.DataFrame(sens).to_csv(out / "table_sensitivity.csv", index=False)
    except SystemExit:
        pass

    # exploratory
    ex = exploratory(con, args.target, df)
    for k in ("excel_share", "format_mix", "file_problems", "file_problem_reasons", "top_symbols"):
        ex[k].to_csv(out / ("exploratory_%s.csv" % k), index=False)
    summary["duplicate_articles_in_sample"] = ex["dup_check"]

    # record-level table (for the Zenodo deposit)
    cols = ["pmcid", "year", "order_idx", "journal", "publisher", "article_type", "license", "supp_formats",
            "n_excel_files", "n_files_scanned", "n_files_unparseable", "n_files_too_large", "species",
            "human_only", "n_errors", "any_error_broad", "any_error_strict", "any_renamed_human",
            "any_non_renamed_human", "any_mouse", "any_other"]
    df[cols].to_csv(out / "record_level_articles.csv", index=False)

    # validation scores, if present
    vs = out / "validation_scores.json"
    summary["validation"] = json.loads(vs.read_text()) if vs.exists() else "not done yet (03_validation.py)"

    co = out / "calibration_original_summary.json"
    summary["calibration_original_articles"] = json.loads(co.read_text()) if co.exists() else "not run (05_calibrate_original.py)"
    n_min = int(df.groupby("year").size().min())
    p0 = float(df.any_error_broad.mean())
    summary["precision"] = {"n_per_year_min": n_min, "overall_rate": p0,
                            "ci_half_width_points_per_year": 100 * 1.96 * (p0 * (1 - p0) / n_min) ** 0.5}

    # download decision
    per_year = df.groupby("year").size()
    worst = int(per_year.min())
    summary["eligible_per_year"] = {int(k): int(v) for k, v in per_year.items()}
    summary["years_missing"] = [y for y in config.YEARS if y not in per_year.index]
    summary["download_decision"] = ("not applicable: this run targets %d per year (GO line is 800)" % args.target
                                    if args.target < 800 else "GO" if worst >= 800 and not summary["years_missing"]
                                    else "BORDERLINE: widen the sample" if worst >= 400 and not summary["years_missing"]
                                    else "NO-GO after two tries, or a year is missing")
    summary["before_after_broad"] = ba_rows[0] if ba_rows else None
    summary["before_after_all"] = ba_rows
    summary["sensitivity"] = sens
    (out / "results_summary.json").write_text(json.dumps(summary, indent=2, default=float))

    # readable summary
    L = ["# Results summary (generated by 04_analyze.py)", "",
         "Numbers are expectations of nothing: they are whatever the data gave.", "",
         "## Sample", ""]
    L += ["- %d: %d eligible articles" % (y, n) for y, n in per_year.items()]
    L += ["- Download decision: **%s**" % summary["download_decision"], "", "## Yearly rates (broad definition)", ""]
    for _, r in broad.iterrows():
        L.append("- %d: %d of %d (%.1f%%, 95%% CI %.1f to %.1f)" % (
            r.year, r.n_error, r.n_eligible, 100 * r.rate, 100 * r.ci_low, 100 * r.ci_high))
    if cal:
        L += ["", "## Calibration against the published series", ""]
        for y in config.BEFORE_YEARS:
            if y in cal:
                v = cal[y]
                L.append("- %d: ours %.1f%%, published confirmed %.1f%% (suspected %.1f%%)" % (
                    y, 100 * v["ours"], 100 * v["published_confirmed"], 100 * v["published_suspected"]))
                if v.get("overlap_n"):
                    L.append("  - on the %d sampled articles the 2021 study also had as gene-list articles (%.0f%% of ours): %.1f%%" % (
                        v["overlap_n"], 100 * v["overlap_share_of_our_sample"], 100 * v["ours_rate_on_overlap"]))
        L.append("- Decision: **%s** (largest gap %.1f points)" % (cal["decision"], cal["largest_gap_points"]))
    L += ["", "## Before and after (2019-20 vs 2024-25)", ""]
    for r in ba_rows:
        L.append("- %s: %s %s" % (r["comparison"], fmt_ba(r), r["small_count_note"]))
    t = summary["trend_broad"]
    if "odds_ratio_per_year" in t:
        L += ["", "## Trend", "",
              "- odds ratio per year %.3f (95%% CI %.3f to %.3f), p=%.3g; non-linearity test p=%.3g" % (
                  t["odds_ratio_per_year"], t["or_ci_low"], t["or_ci_high"], t["p_trend"], t["nonlinearity_p"])]
    d = summary["diff_of_changes_renamed_vs_other"]
    if d:
        L.append("- change in renamed-gene errors minus change in other types: %+.1f points (95%% CI %.1f to %.1f)" % (
            d["diff_of_changes_points"], d["ci_low"], d["ci_high"]))
    if adj and "adjusted_or" in adj:
        L.append("- before/after odds ratio, crude %.2f (%.2f to %.2f); adjusted for publisher mix %.2f (%.2f to %.2f)" % (
            adj["crude_or"], *adj["crude_ci"], adj["adjusted_or"], *adj["adjusted_ci"]))
    co = summary["calibration_original_articles"]
    if isinstance(co, dict):
        v = co["vs"]["orig_confirmed"]
        L += ["", "## Head-to-head on the 2021 study's own articles", "",
              "- %d articles rescanned; gene list found in %.1f%%" % (co["n_scanned"], 100 * co["gene_list_found_share"]),
              "- their confirmed rate %.1f%%, ours %.1f%%, agreement %.1f%%, kappa %.2f" % (
                  100 * v["orig_rate"], 100 * v["ours_rate"], 100 * v["agreement"], v["kappa"])]
    L += ["", "- 95%% CI half-width per year at n=%d: about %.1f points" % (
        summary["precision"]["n_per_year_min"], summary["precision"]["ci_half_width_points_per_year"])]
    L += ["", "## Sensitivity", ""]
    for r in sens:
        L.append("- %s: %s" % (r["comparison"], fmt_ba(r)))
    (out / "results_summary.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    sys.exit(main())
