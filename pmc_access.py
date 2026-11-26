"""Access to PubMed Central: candidate search (E-utilities) and per-article
files (PMC Cloud Service on AWS, https://pmc-oa-opendata.s3.amazonaws.com).

Two interchangeable sources:
  HttpSource   - the real thing
  LocalSource  - a folder laid out like the bucket, used by the tests
"""
from __future__ import annotations

import hashlib
import json
import re
import threading
import time
from pathlib import Path
from urllib.parse import quote, urlparse, parse_qs

import requests

import config

_local = threading.local()
_esearch_lock = threading.Lock()
_last_esearch = [0.0]

UA = "excel_gene_errors/1.0 (research audit; %s)" % config.NCBI_EMAIL


def _session() -> requests.Session:
    s = getattr(_local, "session", None)
    if s is None:
        s = requests.Session()
        s.headers["User-Agent"] = UA
        adapter = requests.adapters.HTTPAdapter(pool_connections=4, pool_maxsize=4)
        s.mount("https://", adapter)
        _local.session = s
    return s


def http_get(url, params=None, headers=None, retries=5, timeout=60, stream=False):
    """GET with retry and back-off on network errors, 429 and 5xx."""
    last = None
    for attempt in range(retries):
        try:
            r = _session().get(url, params=params, headers=headers,
                               timeout=timeout, stream=stream)
            if r.status_code in (429, 500, 502, 503, 504):
                last = RuntimeError("HTTP %d for %s" % (r.status_code, url))
                r.close()
            else:
                return r
        except requests.RequestException as exc:
            last = exc
        time.sleep(min(2 ** attempt, 20))
    raise RuntimeError("giving up on %s: %s" % (url, last))


# ---- E-utilities: candidate lists ------------------------------------------
def _esearch(term, retmax):
    min_gap = 0.11 if config.NCBI_API_KEY else 0.35   # NCBI: 10/s with key, 3/s without
    with _esearch_lock:
        wait = min_gap - (time.monotonic() - _last_esearch[0])
        if wait > 0:
            time.sleep(wait)
        _last_esearch[0] = time.monotonic()
    params = {"db": "pmc", "term": term, "retmax": retmax, "retmode": "json",
              "tool": config.NCBI_TOOL, "email": config.NCBI_EMAIL}
    if config.NCBI_API_KEY:
        params["api_key"] = config.NCBI_API_KEY
    r = http_get(config.EUTILS + "/esearch.fcgi", params=params)
    r.raise_for_status()
    res = r.json()["esearchresult"]
    return int(res["count"]), [int(x) for x in res.get("idlist", [])]


def candidate_ids(year: int, lo: int = 1, hi: int = 20_000_000, _depth=0,
                  query: str | None = None, log=print) -> list[int]:
    """All PMC UIDs matching the candidate query for one publication year.

    esearch returns at most about 10,000 ids per query, so the UID range is
    split in half until each piece is small enough.
    """
    query = query or config.CANDIDATE_QUERY
    term = "(%s) AND %d[pdat] AND %d:%d[uid]" % (query, year, lo, hi)
    count, ids = _esearch(term, 9999)
    if count == 0:
        return []
    if count <= 9999:
        return ids
    if lo >= hi:
        log("warning: %d ids at uid %d, keeping the first 9,999" % (count, lo))
        return ids
    mid = (lo + hi) // 2
    return (candidate_ids(year, lo, mid, _depth + 1, query, log)
            + candidate_ids(year, mid + 1, hi, _depth + 1, query, log))


# ---- article metadata and files --------------------------------------------
def media_name(url: str) -> str:
    """File name from an s3:// url in the metadata JSON (md5 query dropped)."""
    return urlparse(url).path.rsplit("/", 1)[-1]


def media_md5(url: str) -> str | None:
    q = parse_qs(urlparse(url).query)
    return q.get("md5", [None])[0]


class HttpSource:
    def __init__(self, base: str = config.S3_BASE):
        self.base = base.rstrip("/")

    def _get_json(self, key):
        r = http_get("%s/%s" % (self.base, key))
        if r.status_code in (403, 404):
            return None
        r.raise_for_status()
        return r.json()

    def metadata(self, pmcid: int):
        """Metadata JSON of the newest version of an article, or None."""
        meta = self._get_json("metadata/PMC%d.1.json" % pmcid)
        if meta is not None:
            return meta
        # version 1 missing: look for other versions
        r = http_get("%s/" % self.base, params={
            "list-type": "2", "prefix": "PMC%d." % pmcid, "delimiter": "/"})
        r.raise_for_status()
        versions = re.findall(r"<Prefix>PMC%d\.(\d+)/</Prefix>" % pmcid, r.text)
        if not versions:
            return None
        return self._get_json("metadata/PMC%d.%d.json" % (pmcid, max(map(int, versions))))

    def xml_head(self, pmcid: int, version: int, nbytes: int = 60000) -> str:
        r = http_get("%s/PMC%d.%d/PMC%d.%d.xml" % (self.base, pmcid, version, pmcid, version),
                     headers={"Range": "bytes=0-%d" % (nbytes - 1)})
        if r.status_code not in (200, 206):
            return ""
        return r.text[:nbytes]

    def download(self, pmcid: int, version: int, name: str, dest: Path,
                 md5: str | None = None, max_bytes: int | None = None) -> int:
        """Stream a file to dest. Returns bytes written, or -1 if over max_bytes."""
        url = "%s/PMC%d.%d/%s" % (self.base, pmcid, version, quote(name))
        r = http_get(url, stream=True, timeout=120)
        if r.status_code != 200:
            r.close()
            raise RuntimeError("HTTP %d for %s" % (r.status_code, url))
        size = int(r.headers.get("Content-Length") or 0)
        if max_bytes and size > max_bytes:
            r.close()
            return -1
        h = hashlib.md5()
        n = 0
        with open(dest, "wb") as fh:
            for chunk in r.iter_content(1 << 16):
                fh.write(chunk)
                h.update(chunk)
                n += len(chunk)
                if max_bytes and n > max_bytes:
                    r.close()
                    return -1
        if md5 and h.hexdigest() != md5:
            raise RuntimeError("md5 mismatch for %s" % url)
        return n


class LocalSource:
    """Reads a folder laid out like the bucket (used by the tests)."""

    def __init__(self, root: Path):
        self.root = Path(root)

    def metadata(self, pmcid: int):
        p = self.root / "metadata" / ("PMC%d.1.json" % pmcid)
        return json.loads(p.read_text()) if p.exists() else None

    def xml_head(self, pmcid: int, version: int, nbytes: int = 60000) -> str:
        p = self.root / ("PMC%d.%d" % (pmcid, version)) / ("PMC%d.%d.xml" % (pmcid, version))
        return p.read_text()[:nbytes] if p.exists() else ""

    def download(self, pmcid, version, name, dest, md5=None, max_bytes=None) -> int:
        src = self.root / ("PMC%d.%d" % (pmcid, version)) / name
        data = src.read_bytes()
        if max_bytes and len(data) > max_bytes:
            return -1
        Path(dest).write_bytes(data)
        return len(data)


# ---- article XML front matter ----------------------------------------------
_RE_ARTICLE_TYPE = re.compile(r"<article\b[^>]*?\barticle-type=\"([^\"]+)\"")
_RE_JOURNAL = re.compile(r"<journal-title>(.*?)</journal-title>", re.S)
_RE_NLM_TA = re.compile(r"<journal-id[^>]*nlm-ta[^>]*>(.*?)</journal-id>", re.S)
_RE_PUBLISHER = re.compile(r"<publisher-name>(.*?)</publisher-name>", re.S)
_TAG = re.compile(r"<[^>]+>")


def parse_xml_head(text: str) -> dict:
    def first(rx):
        m = rx.search(text)
        return _TAG.sub("", m.group(1)).strip() if m else ""
    m = _RE_ARTICLE_TYPE.search(text)
    return {
        "article_type": m.group(1).strip().lower() if m else "",
        "journal": first(_RE_JOURNAL) or first(_RE_NLM_TA),
        "publisher": first(_RE_PUBLISHER),
    }
