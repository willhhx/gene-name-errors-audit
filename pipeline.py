"""Per-article processing and the sampling loop (used by 02_scan_articles.py)."""
from __future__ import annotations

import datetime as dt
import random
import shutil
import tempfile
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from pathlib import Path

import config
import gene_scan
import pmc_access
import store

IMAGE_EXT = {"jpg", "jpeg", "gif", "png", "tif", "tiff", "svg", "bmp", "eps"}
EXCEL_EXT = ("xls", "xlsx")


def ext_of(name: str) -> str:
    return name.rsplit(".", 1)[-1].lower() if "." in name else ""


def eligible_orders(con, year):
    """order_idx of eligible articles (scanned, with a gene list) in order."""
    return [r[0] for r in con.execute(
        "SELECT order_idx FROM articles WHERE year=? AND status='scanned' AND has_gene_list=1 "
        "ORDER BY order_idx", (year,))]


def _summarise(art, file_rows, errors):
    scanned = [f for f in file_rows if f["parse_status"] == "ok"]
    gl = [f for f in scanned if f["n_gene_columns"] > 0]
    art["n_files_scanned"] = len(scanned)
    art["n_files_unparseable"] = sum(f["parse_status"].startswith("unparseable")
                                     or f["parse_status"] == "timeout" for f in file_rows)
    art["n_files_too_large"] = sum(f["parse_status"] == "too_large" for f in file_rows)
    art["n_gene_list_files"] = len(gl)
    art["has_gene_list"] = int(bool(gl))
    species = set()
    for f in gl:
        species.update(s for s in f["species"].split(",") if s)
    art["species"] = ",".join(sorted(species))
    art["human_only"] = int(species == {"Hsapiens"})
    art["n_errors"] = len(errors)
    art["any_error_broad"] = int(bool(errors))
    art["any_error_strict"] = int(any(e["form"] in ("date", "date_string") for e in errors))
    cats = {e["category"] for e in errors}
    art["any_renamed_human"] = int("renamed_human" in cats)
    art["any_non_renamed_human"] = int("non_renamed_human" in cats)
    art["any_mouse"] = int("mouse" in cats)
    art["any_other"] = int("other" in cats)


def _scan_one(path, pmcid, display_name, source, ctx, file_rows, errors):
    res = gene_scan.scan_workbook(path, ctx["sets"], ctx["hgnc"])
    cols = res.gene_columns
    file_rows.append({
        "pmcid": pmcid, "name": display_name[:200], "ext": ext_of(display_name),
        "source": source, "bytes": Path(path).stat().st_size if Path(path).exists() else 0,
        "parse_status": res.parse_status, "n_sheets": res.n_sheets,
        "n_gene_columns": len(cols),
        "species": ",".join(sorted({c.species for c in cols})),
        "n_cells_checked": res.n_cells_checked, "n_errors": len(res.errors)})
    for e in res.errors:
        errors.append({"pmcid": pmcid, "file": display_name[:200], "sheet": e.sheet[:60],
                       "row": e.row, "col": e.col, "species": e.species, "raw": e.raw,
                       "form": e.form, "month": e.month, "day": e.day,
                       "category": e.category, "candidate": e.candidate})


def process_article(pmcid: int, year: int, order_idx: int, ctx: dict):
    """Fetch, scan and summarise one article. Returns (article, files, errors)."""
    src = ctx["source"]
    art = {"pmcid": pmcid, "year": year, "order_idx": order_idx, "status": "",
           "scanned_at": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")}
    files, errors = [], []
    try:
        meta = src.metadata(pmcid)
    except Exception as exc:
        art.update(status="fetch_error", note=str(exc)[:200])
        return art, files, errors
    if meta is None:
        art["status"] = "not_found"
        return art, files, errors

    version = int(meta.get("version") or 1)
    art.update(license=meta.get("license_code"), doi=meta.get("doi"),
               is_retracted=int(bool(meta.get("is_retracted"))))
    if not meta.get("is_pmc_openaccess", True):
        art["status"] = "not_open_access"
        return art, files, errors

    media = [(pmc_access.media_name(u), u) for u in meta.get("media_urls") or []]
    counts = Counter(ext_of(n) for n, _ in media if ext_of(n) not in IMAGE_EXT)
    art["supp_formats"] = ";".join("%s=%d" % kv for kv in sorted(counts.items()))
    excel = [(n, u) for n, u in media if ext_of(n) in EXCEL_EXT]
    zips = [(n, u) for n, u in media if ext_of(n) == "zip"]
    art["n_excel_files"] = len(excel)
    art["n_zip_files"] = len(zips)
    if not excel and not zips:
        art["status"] = "no_excel"
        return art, files, errors

    try:
        front = pmc_access.parse_xml_head(src.xml_head(pmcid, version))
    except Exception as exc:
        art.update(status="fetch_error", note="xml: " + str(exc)[:190])
        return art, files, errors
    art.update(article_type=front["article_type"], journal=front["journal"][:120],
               publisher=front["publisher"][:120])
    if front["article_type"] in config.EXCLUDED_ARTICLE_TYPES:
        art["status"] = "excluded_type"
        return art, files, errors

    tmp = Path(tempfile.mkdtemp(prefix="gne_"))
    max_bytes = config.MAX_FILE_MB * 1024 * 1024
    try:
        for i, (name, url) in enumerate(excel):
            dest = tmp / ("f%d.%s" % (i, ext_of(name)))
            n = src.download(pmcid, version, name, dest, pmc_access.media_md5(url), max_bytes)
            if n < 0:
                files.append({"pmcid": pmcid, "name": name[:200], "ext": ext_of(name),
                              "source": "direct", "bytes": 0, "parse_status": "too_large",
                              "n_sheets": 0, "n_gene_columns": 0, "species": "",
                              "n_cells_checked": 0, "n_errors": 0})
                continue
            _scan_one(dest, pmcid, name, "direct", ctx, files, errors)
            dest.unlink(missing_ok=True)
        for j, (name, url) in enumerate(zips):
            zdest = tmp / ("z%d.zip" % j)
            n = src.download(pmcid, version, name, zdest, pmc_access.media_md5(url), 2 * max_bytes)
            if n < 0:
                files.append({"pmcid": pmcid, "name": name[:200], "ext": "zip", "source": "direct",
                              "bytes": 0, "parse_status": "too_large", "n_sheets": 0,
                              "n_gene_columns": 0, "species": "", "n_cells_checked": 0, "n_errors": 0})
                continue
            zdir = tmp / ("zx%d" % j)
            zdir.mkdir()
            try:
                members = list(gene_scan.iter_excel_in_zip(zdest, zdir))
            except gene_scan.ParseError as exc:
                files.append({"pmcid": pmcid, "name": name[:200], "ext": "zip", "source": "direct",
                              "bytes": n, "parse_status": "unparseable:" + str(exc), "n_sheets": 0,
                              "n_gene_columns": 0, "species": "", "n_cells_checked": 0, "n_errors": 0})
                continue
            for mname, mpath in members:
                disp = "%s::%s" % (name, mname)
                if mpath is None:
                    files.append({"pmcid": pmcid, "name": disp[:200], "ext": ext_of(mname),
                                  "source": "zip", "bytes": 0, "parse_status": "too_large",
                                  "n_sheets": 0, "n_gene_columns": 0, "species": "",
                                  "n_cells_checked": 0, "n_errors": 0})
                    continue
                _scan_one(mpath, pmcid, disp, "zip", ctx, files, errors)
                mpath.unlink(missing_ok=True)
            zdest.unlink(missing_ok=True)
    except Exception as exc:
        art.update(status="fetch_error", note=str(exc)[:200])
        return art, [], []
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    art["status"] = "scanned"
    _summarise(art, files, errors)
    return art, files, errors


def candidate_order(year: int, ids: list[int]) -> list[int]:
    """The fixed random order of candidates for a year (seeded)."""
    ids = sorted(set(ids))
    random.Random(config.SEED + year).shuffle(ids)
    return ids


_WORKER = {}


def init_worker(spec: dict):
    """Build the scanner context once per worker process."""
    hgnc = gene_scan.load_hgnc(config.REF / "hgnc_complete_set.txt")
    sets = gene_scan.load_gene_sets(config.REF, hgnc, original_only=spec.get("original_only", False))
    source = (pmc_access.LocalSource(spec["source_dir"]) if spec.get("source_dir")
              else pmc_access.HttpSource())
    _WORKER["ctx"] = {"source": source, "sets": sets, "hgnc": hgnc}
    _WORKER["threads"] = spec.get("threads", 8)


def work_chunk(items):
    """items: list of (order_idx, pmcid, year). Threads hide network waiting;
    several worker processes spread the spreadsheet parsing over CPU cores."""
    ctx = _WORKER["ctx"]
    with ThreadPoolExecutor(max_workers=_WORKER["threads"]) as ex:
        return list(ex.map(lambda it: process_article(it[1], it[2], it[0], ctx), items))


def run_year(con, year: int, ids: list[int], spec: dict, target: int,
             max_candidates: int | None = None, procs: int = 1,
             chunk: int = 16, log=print):
    """Process candidates in seeded random order until `target` eligible
    articles are found (or candidates / max_candidates run out)."""
    order = candidate_order(year, ids)
    if max_candidates:
        order = order[:max_candidates]
    done = {r[0] for r in con.execute(
        "SELECT pmcid FROM articles WHERE year=? AND status!='fetch_error'", (year,))}
    todo = [(i, p, year) for i, p in enumerate(order) if p not in done]
    n_elig = len(eligible_orders(con, year))
    log("%d: %d candidates, %d already done, %d eligible so far (target %d)" % (
        year, len(order), len(done), n_elig, target))
    t0 = time.time()
    n_new = 0
    batch = chunk * max(procs, 1) * 2
    if procs <= 1:
        init_worker(spec)
        run_chunks = lambda chunks: [work_chunk(c) for c in chunks]
        pool = None
    else:
        import multiprocessing as mp
        pool = ProcessPoolExecutor(max_workers=procs, mp_context=mp.get_context("spawn"),
                                   initializer=init_worker, initargs=(spec,))
        run_chunks = lambda chunks: list(pool.map(work_chunk, chunks))
    try:
        for k in range(0, len(todo), batch):
            if n_elig >= target:
                break
            part = todo[k:k + batch]
            chunks = [part[j:j + chunk] for j in range(0, len(part), chunk)]
            for results in run_chunks(chunks):
                for art, files, errs in results:
                    store.save_result(con, art, files, errs)
            con.commit()
            n_new += len(part)
            n_elig = len(eligible_orders(con, year))
            rate = n_new / max(time.time() - t0, 1e-9)
            log("%d: %d/%d candidates done, %d eligible, %.1f candidates/s" % (
                year, len(order) - len(todo) + k + len(part), len(order), n_elig, rate))
    finally:
        if pool:
            pool.shutdown()
    return n_elig
