"""Synthetic multi-year database -> run 04_analyze.py and check outputs."""
import importlib.util
import json
import random
import sys
from pathlib import Path

import pytest

import config
import store

ROOT = Path(__file__).resolve().parent.parent


def load_script(name):
    spec = importlib.util.spec_from_file_location(name.replace(".py", "").replace("0", "s", 1),
                                                  ROOT / name)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def synth_db(path, n_per_year=300, rates=None):
    rates = rates or {2019: .32, 2020: .31, 2021: .30, 2022: .28, 2023: .25, 2024: .18, 2025: .15}
    con = store.connect(path)
    rng = random.Random(1)
    pid = 1000
    for y in config.YEARS:
        for i in range(n_per_year):
            pid += 1
            err = rng.random() < rates[y]
            ren = err and rng.random() < (0.6 if y <= 2020 else 0.25)
            oth = err and (not ren or rng.random() < 0.3)
            art = {"pmcid": pid, "year": y, "order_idx": i, "status": "scanned",
                   "supp_formats": "xlsx=1;pdf=1", "n_excel_files": 1,
                   "journal": "Nucleic Acids Research" if i % 5 == 0 else "Sci Rep",
                   "article_type": "research-article", "has_gene_list": 1,
                   "publisher": "Pub A" if i % 3 else "Pub B",
                   "n_gene_list_files": 1, "species": "Hsapiens", "human_only": 1,
                   "n_errors": 1 if err else 0, "any_error_broad": int(err),
                   "any_error_strict": int(err), "any_renamed_human": int(ren),
                   "any_non_renamed_human": int(oth), "any_mouse": 0, "any_other": 0,
                   "n_files_scanned": 1}
            errs = []
            if err:
                errs.append({"pmcid": pid, "file": "a.xlsx", "sheet": "S", "row": 5, "col": 1,
                             "species": "Hsapiens", "raw": "2023-09-02", "form": "date", "month": 9,
                             "day": 2, "category": "renamed_human" if ren else "non_renamed_human",
                             "candidate": "SEPT2" if ren else "OCT4"})
            files = [{"pmcid": pid, "name": "a.xlsx", "ext": "xlsx", "source": "direct", "bytes": 1,
                      "parse_status": "ok", "n_sheets": 1, "n_gene_columns": 1, "species": "Hsapiens",
                      "n_cells_checked": 10, "n_errors": len(errs)}]
            store.save_result(con, art, files, errs)
        # some candidates without Excel
        for i in range(n_per_year, n_per_year + 200):
            pid += 1
            store.save_result(con, {"pmcid": pid, "year": y, "order_idx": i - n_per_year, "status": "no_excel",
                                    "supp_formats": "pdf=1"}, [], [])
    con.commit()
    con.close()


def test_analysis_runs_and_finds_the_planted_decline(tmp_path, monkeypatch):
    db = tmp_path / "s.sqlite"
    synth_db(db)
    out = tmp_path / "out"
    mod = load_script("04_analyze.py")
    monkeypatch.setattr(sys, "argv", ["04_analyze.py", "--db", str(db), "--target", "300",
                                      "--out", str(out), "--db-orig", str(tmp_path / "none.sqlite")])
    mod.main()
    s = json.loads((out / "results_summary.json").read_text())
    ba = s["before_after_broad"]
    assert ba["before_n"] == 600 and ba["after_n"] == 600
    assert ba["diff_points"] < -8 and ba["p_value"] < 0.001      # planted: ~31% -> ~16%
    assert s["trend_broad"]["odds_ratio_per_year"] < 1 and s["trend_broad"]["p_trend"] < 0.001
    assert s["download_decision"].startswith("not applicable")   # target 300 is a pilot size
    assert s["duplicate_articles_in_sample"] == 0
    d = s["diff_of_changes_renamed_vs_other"]
    assert d["diff_of_changes_points"] < 0
    for f in ["table1_yearly.csv", "table_before_after.csv", "table_sensitivity.csv",
              "fig1_yearly_rate.png", "fig1_yearly_rate.pdf", "fig2_error_types.png",
              "exploratory_excel_share.csv", "exploratory_format_mix.csv",
              "exploratory_top_symbols.csv", "record_level_articles.csv", "results_summary.md"]:
        assert (out / f).exists() and (out / f).stat().st_size > 0, f
    import pandas as pd
    t1 = pd.read_csv(out / "table1_yearly.csv")
    assert list(t1.year) == config.YEARS and (t1.n_eligible == 300).all()
    ex = pd.read_csv(out / "exploratory_excel_share.csv")
    assert abs(ex.share.iloc[0] - 300 / 500) < 1e-9
    adj = s["publisher_adjusted"]
    assert adj["adjusted_or"] < 1 and adj["adjusted_ci"][1] < 1 and adj["publishers_in_model"] == 2
    assert s["precision"]["n_per_year_min"] == 300
    sens = pd.read_csv(out / "table_sensitivity.csv")
    assert any("directly" in c for c in sens.comparison)
    assert any("2016" in c for c in sens.comparison)


def test_target_keeps_first_n_eligible_only(tmp_path, monkeypatch):
    db = tmp_path / "s.sqlite"
    synth_db(db, n_per_year=100)
    mod = load_script("04_analyze.py")
    con = store.connect(db)
    df = mod.load_sample(con, 40)
    assert (df.groupby("year").size() == 40).all()
    assert df.order_idx.max() == 39


def test_validation_score(tmp_path, monkeypatch, capsys):
    import csv
    mod = load_script("03_validation.py")
    sheet = tmp_path / "v.csv"
    rows = [{"id": i, "group": "flagged", "label": "1" if i < 46 else "0"} for i in range(50)]
    rows += [{"id": 100 + i, "group": "unflagged", "label": "1" if i == 0 else "0"} for i in range(40)]
    with open(sheet, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["id", "group", "label"])
        w.writeheader()
        w.writerows(rows)
    monkeypatch.setattr(mod, "SHEET", sheet)
    monkeypatch.setattr(config, "OUT", tmp_path)
    mod.score(None)
    s = json.loads((tmp_path / "validation_scores.json").read_text())
    assert abs(s["ppv"]["value"] - 0.92) < 1e-9 and s["decision"] == "GO"
    assert abs(s["miss_rate"]["value"] - 0.025) < 1e-9
