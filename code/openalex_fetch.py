"""Builds the co-authorship dataset (institution coordinates as
positions, productivity marks) from the OpenAlex API.

The date range is split into quarterly slices fetched concurrently,
and the institution geo lookups run concurrently as well.  The work
is network-bound, not CPU-bound: a global limiter keeps the request
rate at --rps (default 5 per second, safely under the API's
polite-pool ceiling of about 10), so throughput saturates at a
handful of effective streams no matter how many workers are
configured.

    python3 openalex_fetch.py --mailto you@example.edu \
        [--country NL] [--concept C121332964] \
        [--from-year 2018] [--to-year 2022] [--workers 30]

Outputs, written to the current directory:
    oa_homes.csv.gz    user_id,lat,lon,checkins
    oa_edges.txt.gz    tab-separated co-authorship edges (id1 < id2)
"""
import argparse, gzip, json, threading, time
import urllib.request, urllib.parse, urllib.error
from concurrent.futures import ThreadPoolExecutor, as_completed

API = "https://api.openalex.org"

class RateLimiter:
    def __init__(self, rps):
        self.min_int = 1.0 / max(rps, 0.1)
        self.lock = threading.Lock()
        self.next_t = 0.0
    def wait(self):
        with self.lock:
            now = time.monotonic()
            t = max(self.next_t, now)
            self.next_t = t + self.min_int
        d = t - now
        if d > 0:
            time.sleep(d)

class ApiError(Exception):
    def __init__(self, code, body, url):
        super().__init__(f"HTTP {code} for {url[:160]} :: {body[:300]}")
        self.code = code

STATS = {"ok": 0, "retries": 0, "last": ""}
STATS_LOCK = threading.Lock()

def _note(ok=False, retry=None):
    with STATS_LOCK:
        if ok:
            STATS["ok"] += 1
        if retry:
            STATS["retries"] += 1
            STATS["last"] = retry

def get(url, mailto, rl, tries=5):
    q = ("&" if "?" in url else "?") + urllib.parse.urlencode(
        {"mailto": mailto})
    req = urllib.request.Request(
        url + q, headers={"User-Agent":
                          f"colas-fetch/1.0 (mailto:{mailto})"})
    for t in range(tries):
        rl.wait()
        try:
            with urllib.request.urlopen(req, timeout=25) as r:
                out = json.load(r)
                _note(ok=True)
                return out
        except urllib.error.HTTPError as e:
            if e.code == 429 or e.code >= 500:
                ra = e.headers.get("Retry-After")
                try:
                    wait_s = float(ra)
                except (TypeError, ValueError):
                    wait_s = 2.0 ** (t + 1)
                _note(retry=f"HTTP {e.code}, waiting {wait_s:.0f}s")
                time.sleep(wait_s)
                continue
            try:
                body = e.read().decode("utf-8", "replace")
            except Exception:
                body = ""
            raise ApiError(e.code, body, url)   # 4xx: fail fast, say why
        except Exception as e:
            if t == tries - 1:
                raise
            wait_s = 2.0 ** (t + 1)
            _note(retry=f"{type(e).__name__}: {e}; waiting {wait_s:.0f}s")
            time.sleep(wait_s)
    raise RuntimeError("retries exhausted: " + url)

def heartbeat(t0):
    def loop():
        while True:
            time.sleep(20)
            with STATS_LOCK:
                msg = (f"[oa] alive {time.time()-t0:.0f}s: "
                       f"{STATS['ok']} requests ok, "
                       f"{STATS['retries']} retried"
                       + (f" (last: {STATS['last']})"
                          if STATS["last"] else ""))
            print(msg, flush=True)
    th = threading.Thread(target=loop, daemon=True)
    th.start()

def probe(a, rl):
    """One request before anything else; explain failure plainly."""
    try:
        d = get(f"{API}/works?per-page=1", a.mailto, rl, tries=2)
        n = (d.get("meta") or {}).get("count", "?")
        print(f"[oa] probe ok: API reachable ({n} works total)",
              flush=True)
    except ApiError as e:
        raise SystemExit(
            f"[oa] PROBE FAILED: the API refused this machine "
            f"(HTTP {e.code}).  This is usually a temporary rate-limit "
            f"penalty from earlier failed runs; wait 15-30 minutes and "
            f"rerun (checkpointed slices are not refetched).  "
            f"Detail: {e}")
    except Exception as e:
        raise SystemExit(
            f"[oa] PROBE FAILED: cannot reach {API} "
            f"({type(e).__name__}: {e}).  Check connectivity/proxy, "
            f"or wait and rerun; completed slices are checkpointed.")

class Counter:
    def __init__(self):
        self.n = 0
        self.lock = threading.Lock()
    def add(self, k):
        with self.lock:
            self.n += k
            return self.n

def quarters(y0, y1):
    out = []
    for y in range(y0, y1 + 1):
        out += [(f"{y}-01-01", f"{y}-03-31"), (f"{y}-04-01", f"{y}-06-30"),
                (f"{y}-07-01", f"{y}-09-30"), (f"{y}-10-01", f"{y}-12-31")]
    return out

CKPT_DIR = ".oa_ckpt"

def fetch_slice(d0, d1, a, rl, counter, stop, t0, plock):
    import os, pickle
    ck = os.path.join(CKPT_DIR, f"slice_{d0}.pkl")
    if os.path.exists(ck):
        with open(ck, "rb") as fh:
            aw, ai, edges, works = pickle.load(fh)
        counter.add(works)
        with plock:
            print(f"[oa] slice {d0} loaded from checkpoint "
                  f"({works} works)", flush=True)
        return aw, ai, edges, works
    filt = (f"institutions.country_code:{a.country},"
            f"concepts.id:{a.concept},"
            f"from_publication_date:{d0},"
            f"to_publication_date:{d1},"
            f"type:article")
    cursor = "*"
    aw, ai, edges, works = {}, {}, set(), 0
    while cursor and not stop.is_set():
        d = get(f"{API}/works?filter={filt}&select=id,authorships"
                f"&per-page=200&cursor={cursor}", a.mailto, rl)
        res = d.get("results", [])
        for wrec in res:
            auths = []
            for au in wrec.get("authorships", []):
                aid = (au.get("author") or {}).get("id")
                if not aid:
                    continue
                aid = aid.rsplit("/", 1)[-1]
                auths.append(aid)
                aw[aid] = aw.get(aid, 0) + 1
                for inst in au.get("institutions", []):
                    iid = (inst.get("id") or "").rsplit("/", 1)[-1]
                    if iid:
                        ai.setdefault(aid, {})
                        ai[aid][iid] = ai[aid].get(iid, 0) + 1
            auths = sorted(set(auths))
            if len(auths) > 25:      # drop hyper-authored consortia papers
                continue
            for i in range(len(auths)):
                for j in range(i + 1, len(auths)):
                    edges.add((auths[i], auths[j]))
            works += 1
        cursor = d.get("meta", {}).get("next_cursor")
        total = counter.add(len(res))
        if total >= a.max_works:
            stop.set()
        if res and (total // 5000) != ((total - len(res)) // 5000):
            with plock:
                print(f"[oa] ~{total} works fetched, "
                      f"{time.time()-t0:.0f}s", flush=True)
    import os, pickle
    os.makedirs(CKPT_DIR, exist_ok=True)
    with open(ck + ".tmp", "wb") as fh:
        pickle.dump((aw, ai, edges, works), fh)
    os.replace(ck + ".tmp", ck)
    return aw, ai, edges, works

def _geo_of(rec):
    gg = rec.get("geo") or {}
    if gg.get("latitude") is not None:
        return (gg["latitude"], gg["longitude"])
    return None

def fetch_geo(inst_ids, a, rl, workers):
    chunks = [inst_ids[k:k + 50] for k in range(0, len(inst_ids), 50)]
    geo = {}
    t0 = time.time()

    def one(ch, key):
        d = get(f"{API}/institutions?filter={key}:{'|'.join(ch)}"
                f"&select=id,geo&per-page=50", a.mailto, rl)
        out = {}
        for rec in d.get("results", []):
            g = _geo_of(rec)
            if g:
                out[rec["id"].rsplit("/", 1)[-1]] = g
        return out

    # probe which batch-filter key this API version accepts
    key_ok = None
    for key in ("openalex_id", "ids.openalex", "openalex"):
        try:
            geo.update(one(chunks[0], key))
            key_ok = key
            print(f"[oa] geo batch filter '{key}' accepted", flush=True)
            break
        except ApiError as e:
            print(f"[oa] geo batch filter '{key}' rejected "
                  f"(HTTP {e.code}); trying next", flush=True)

    if key_ok is not None:
        done = 1
        ex = ThreadPoolExecutor(max_workers=workers)
        try:
            futs = [ex.submit(one, c, key_ok) for c in chunks[1:]]
            for fut in as_completed(futs):
                geo.update(fut.result())
                done += 1
                if done % 50 == 0:
                    print(f"[oa] geo {done}/{len(chunks)} chunks, "
                          f"{time.time()-t0:.0f}s", flush=True)
        except BaseException:
            ex.shutdown(wait=False, cancel_futures=True)
            raise
        ex.shutdown(wait=True)
        return geo

    # fallback: one request per institution (slower, always valid)
    print("[oa] falling back to per-institution geo lookups "
          f"({len(inst_ids)} requests)", flush=True)
    def single(iid):
        d = get(f"{API}/institutions/{iid}?select=id,geo", a.mailto, rl)
        g = _geo_of(d)
        return (iid, g)
    done = 0
    ex = ThreadPoolExecutor(max_workers=workers)
    try:
        futs = [ex.submit(single, i) for i in inst_ids]
        for fut in as_completed(futs):
            iid, g = fut.result()
            if g:
                geo[iid] = g
            done += 1
            if done % 2000 == 0:
                print(f"[oa] geo {done}/{len(inst_ids)} institutions, "
                      f"{time.time()-t0:.0f}s", flush=True)
    except BaseException:
        ex.shutdown(wait=False, cancel_futures=True)
        raise
    ex.shutdown(wait=True)
    return geo

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--country", default="NL")
    ap.add_argument("--concept", default="C121332964")
    ap.add_argument("--from-year", type=int, default=2018)
    ap.add_argument("--to-year", type=int, default=2022)
    ap.add_argument("--mailto", required=True)
    ap.add_argument("--max-works", type=int, default=120000)
    ap.add_argument("--min-works-per-author", type=int, default=2)
    ap.add_argument("--workers", type=int, default=12,
                    help="concurrent fetch threads (throughput is capped "
                         "by --rps, not by thread count)")
    ap.add_argument("--rps", type=float, default=5.0,
                    help="global request rate limit per second")
    a = ap.parse_args()

    rl = RateLimiter(a.rps)
    counter, stop, plock = Counter(), threading.Event(), threading.Lock()
    slices = quarters(a.from_year, a.to_year)
    t0 = time.time()
    print(f"[oa] {len(slices)} quarterly slices, {a.workers} workers, "
          f"{a.rps} requests/s", flush=True)
    probe(a, rl)
    heartbeat(t0)

    author_works, author_inst, edges = {}, {}, set()
    total_works = 0
    with ThreadPoolExecutor(max_workers=a.workers) as ex:
        futs = [ex.submit(fetch_slice, d0, d1, a, rl, counter, stop,
                          t0, plock) for d0, d1 in slices]
        for fut in as_completed(futs):
            aw, ai, ed, w = fut.result()
            total_works += w
            edges |= ed
            for k, v in aw.items():
                author_works[k] = author_works.get(k, 0) + v
            for k, m in ai.items():
                dst = author_inst.setdefault(k, {})
                for iid, c in m.items():
                    dst[iid] = dst.get(iid, 0) + c
    print(f"[oa] works done: {total_works} kept, "
          f"{len(author_works)} authors, {len(edges)} raw edges, "
          f"{time.time()-t0:.0f}s", flush=True)

    inst_ids = sorted({i for m in author_inst.values() for i in m})
    print(f"[oa] fetching geo for {len(inst_ids)} institutions", flush=True)
    geo = fetch_geo(inst_ids, a, rl, a.workers)
    print(f"[oa] geo found for {len(geo)} of {len(inst_ids)} "
          f"institutions", flush=True)

    idmap = {}
    kept = 0
    with gzip.open("oa_homes.csv.gz", "wt") as fh:
        fh.write("user_id,lat,lon,checkins\n")
        for aid, wcount in author_works.items():
            if wcount < a.min_works_per_author:
                continue
            insts = {i: c for i, c in author_inst.get(aid, {}).items()
                     if i in geo}
            if not insts:
                continue
            top = max(insts, key=insts.get)
            la, lo = geo[top]
            idmap[aid] = len(idmap) + 1
            fh.write(f"{idmap[aid]},{la},{lo},{wcount}\n")
            kept += 1
    ne = 0
    with gzip.open("oa_edges.txt.gz", "wt") as fh:
        for u, v in sorted(edges):
            if u in idmap and v in idmap:
                x, y = idmap[u], idmap[v]
                if x > y:
                    x, y = y, x
                fh.write(f"{x}\t{y}\n")
                ne += 1
    print(f"[oa] wrote oa_homes.csv.gz ({kept} authors) and "
          f"oa_edges.txt.gz ({ne} edges), {time.time()-t0:.0f}s total",
          flush=True)

if __name__ == "__main__":
    main()
