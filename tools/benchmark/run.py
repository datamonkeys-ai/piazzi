#!/usr/bin/env python3
"""Reproducible astrometric accuracy benchmark for Piazzi using public data (ZTF/IRSA + Gaia DR3 + SkyBoT + JPL Horizons).

Standard library only. Run with `python3 -I tools/benchmark/run.py ...` from the repository root.

  select  [--n 30] [--seed 1]          generates tools/benchmark/fields.json (fixed, versioned list)
  select-sdss [--n 35]                 generates tools/benchmark/fields_sdss.json (independent TEST set: SDSS DR17, MPC 645)
  run     --label base [--piazzi BIN] [--only ID,...] [--max N]
                                       runs piazzi on each sequence and writes results/<label>.json and .md
  inject  --scenario trail [--fields F] [--n 6] [--seeds 1,2] --label L
                                       synthetic object injection (never on fields_test.json): completeness and O−C
  vs-mpc  --label L [--fields F]      Piazzi vs the station's official astrometry (I41/645) at the MPC -> results/<L>-vsmpc.*
  compare a.json b.json                compares two results
  summary a.json                       reprints the summary

Network cache in tools/benchmark/cache/ (gitignored): FITS cutouts, Gaia catalogs, SkyBoT/Horizons responses.
"""
import argparse
import csv
import hashlib
import io
import json
import math
import os
import random
import re
import statistics
import subprocess
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
CACHE = os.path.join(HERE, "cache")
RESULTS = os.path.join(HERE, "results")
FIELDS = os.path.join(HERE, "fields.json")


def rel_path(p):
    """Path to record in results/reports: relative to the repository root when inside it; otherwise just the
    file name (never an absolute path on the local machine)."""
    if not p:
        return p
    ap = os.path.abspath(p)
    if ap == ROOT or ap.startswith(ROOT + os.sep):
        return os.path.relpath(ap, ROOT)
    return os.path.basename(ap)

IRSA_SEARCH = "https://irsa.ipac.caltech.edu/ibe/search/ztf/products/sci"
IRSA_DATA = "https://irsa.ipac.caltech.edu/ibe/data/ztf/products/sci"
SKYBOT = "https://ssp.imcce.fr/webservices/skybot/api/conesearch.php"
HORIZONS = "https://ssd.jpl.nasa.gov/api/horizons.api"
UA = "piazzi-benchmark/1.0 (+https://github.com/; astrometry benchmark)"

# ZTF configuration (1.01"/px camera, 48" Samuel Oschin, Palomar I41). OBSJD/EXPTIME in the cutout header;
# OBSJD = exposure start (UTC-OBS) and piazzi, with time_reference=0 (start), adds EXPTIME/2 -> mid-exposure.
ZTF_SETTINGS = {
    "focal_length_mm": "3063", "pixel_width_mm": "0.015", "pixel_height_mm": "0.015", "mpc_code": "I41", "time_reference": "0",
    "catalog_faint_mag": "19.0", "saturation": "60000", "observer": "ZTF", "measurer": "Piazzi benchmark",
    "submitter": "Piazzi benchmark", "telescope": "1.2-m f/2.4 Schmidt + CCD",
}
CUTOUT_ARCMIN = 12.0
MATCH_ARCSEC = 10.0

_last_call = {}


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def http_get(url, cache_name=None, binary=False, pause=1.0, tries=6, timeout=120):
    """GET with on-disk cache, per-host pacing and retries with backoff."""
    path = None
    if cache_name:
        path = os.path.join(CACHE, cache_name)
        if os.path.exists(path) and os.path.getsize(path) > 0:
            with open(path, "rb") as f:
                d = f.read()
            return d if binary else d.decode("utf-8", "replace")
    host = urllib.parse.urlparse(url).netloc
    err = None
    for k in range(tries):
        wait = _last_call.get(host, 0) + pause - time.time()
        if wait > 0:
            time.sleep(wait)
        _last_call[host] = time.time()
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                d = r.read()
            if path:
                os.makedirs(os.path.dirname(path), exist_ok=True)
                with open(path + ".tmp", "wb") as f:
                    f.write(d)
                os.replace(path + ".tmp", path)
            return d if binary else d.decode("utf-8", "replace")
        except (urllib.error.URLError, OSError) as e:
            err = e
            code = getattr(e, "code", None)
            if code is not None and 400 <= code < 500 and code != 429:
                break
            time.sleep(min(30, 5 * 2 ** k))
    raise RuntimeError("failed on %s: %s" % (url[:120], err))


def h(s):
    return hashlib.sha1(s.encode()).hexdigest()[:16]


# ------------------------------------------------------------------ simple astrometry
def sep_arcsec(ra1, de1, ra2, de2):
    r1, d1, r2, d2 = map(math.radians, (ra1, de1, ra2, de2))
    c = math.sin(d1) * math.sin(d2) + math.cos(d1) * math.cos(d2) * math.cos(r1 - r2)
    return math.degrees(math.acos(max(-1.0, min(1.0, c)))) * 3600


def hms2deg(s):
    p = [float(x) for x in s.split()]
    return 15 * (p[0] + p[1] / 60 + p[2] / 3600)


def dms2deg(s):
    s = s.strip()
    sg = -1 if s.startswith("-") else 1
    p = [abs(float(x)) for x in s.replace("+", "").replace("-", "").split()]
    return sg * (p[0] + p[1] / 60 + p[2] / 3600)


def ecl2eq(lon, lat):
    e = math.radians(23.4393)
    l, b = math.radians(lon), math.radians(lat)
    sd = math.sin(b) * math.cos(e) + math.cos(b) * math.sin(e) * math.sin(l)
    d = math.asin(sd)
    ra = math.atan2(math.sin(l) * math.cos(e) - math.tan(b) * math.sin(e), math.cos(l))
    return math.degrees(ra) % 360, math.degrees(d)


def jd_to_iso(jd):
    z = jd + 0.5
    Z = int(z)
    F = z - Z
    a = int((Z - 1867216.25) / 36524.25)
    A = Z + 1 + a - a // 4
    B = A + 1524
    C = int((B - 122.1) / 365.25)
    D = int(365.25 * C)
    E = int((B - D) / 30.6001)
    day = B - D - int(30.6001 * E)
    mon = E - 1 if E < 14 else E - 13
    yr = C - 4716 if mon > 2 else C - 4715
    return yr, mon, day, F * 24


def iso_to_jd(s):
    m = re.match(r"(\d+)-(\d+)-(\d+)T(\d+):(\d+):([\d.]+)", s)
    y, mo, d, hh, mi, ss = int(m[1]), int(m[2]), int(m[3]), int(m[4]), int(m[5]), float(m[6])
    if mo <= 2:
        y, mo = y - 1, mo + 12
    a = y // 100
    b = 2 - a + a // 4
    return int(365.25 * (y + 4716)) + int(30.6001 * (mo + 1)) + d + b - 1524.5 + (hh + mi / 60 + ss / 3600) / 24


# ------------------------------------------------------------------ services
def irsa_search(ra, dec, jd0, jd1):
    q = {"POS": "%.5f,%.5f" % (ra, dec), "ct": "csv",
         "COLUMNS": "obsjd,field,ccdid,qid,filtercode,filefracday,exptime,seeing,maglimit,infobits",
         "WHERE": "obsjd>%.1f AND obsjd<%.1f" % (jd0, jd1)}
    url = IRSA_SEARCH + "?" + urllib.parse.urlencode(q)
    t = http_get(url, "irsa/%s.csv" % h(url), pause=1.0)
    return list(csv.DictReader(io.StringIO(t)))


def cutout_url(e, ra, dec, size=CUTOUT_ARCMIN):
    ff = str(e["filefracday"])
    return "%s/%s/%s/%s/ztf_%s_%06d_%s_c%02d_o_q%d_sciimg.fits?center=%.6f,%.6f&size=%garcmin&gzip=false" % (
        IRSA_DATA, ff[:4], ff[4:8], ff[8:], ff, int(e["field"]), e["filtercode"], int(e["ccdid"]), int(e["qid"]),
        ra, dec, size)


def skybot(ra, dec, jd, rad_deg, loc="I41"):
    q = {"-ep": "%.6f" % jd, "-ra": "%.6f" % ra, "-dec": "%.6f" % dec, "-rd": "%.4f" % rad_deg, "-mime": "text",
         "-output": "basic", "-loc": loc, "-filter": "0", "-objFilter": "100"}
    url = SKYBOT + "?" + urllib.parse.urlencode(q)
    t = http_get(url, "skybot/%s.txt" % h(url), pause=2.0)
    out = []
    for line in t.splitlines():
        if line.startswith("#") or "|" not in line:
            continue
        p = [x.strip() for x in line.split("|")]
        if len(p) < 10:
            continue
        try:
            out.append({"num": p[0], "name": p[1], "ra": hms2deg(p[2]), "dec": dms2deg(p[3]), "cls": p[4],
                        "V": float(p[5]), "err": float(p[6]), "dra_h": float(p[8]), "ddec_h": float(p[9])})
        except ValueError:
            pass
    return out


def horizons_command(obj):
    if obj.get("num") and obj["num"] not in ("-", ""):
        return "%s;" % obj["num"]
    return "DES=%s;" % obj["name"]


def horizons(obj, jds, site="I41"):
    """Astrometric RA/Dec (ICRF) from Horizons for the given MPC site (I41 = ZTF, 645 = SDSS) at the given UTC times, with the
    formal ephemeris uncertainty per axis (Horizons quantity 36, RA_3sigma/DEC_3sigma in ″ on the sky plane; here
    divided by 3 -> 1σ). Returns [(ra, dec, sig_ra, sig_dec)]; sig = None when Horizons does not provide it ('n.a.')."""
    q = {"format": "text", "COMMAND": "'%s'" % horizons_command(obj), "CENTER": "'%s@399'" % site, "MAKE_EPHEM": "YES",
         "EPHEM_TYPE": "OBSERVER", "TLIST": "'" + " ".join("%.7f" % j for j in jds) + "'", "TLIST_TYPE": "JD",
         "TIME_TYPE": "UT", "QUANTITIES": "'1,36'", "ANG_FORMAT": "DEG", "EXTRA_PREC": "YES", "CSV_FORMAT": "YES",
         "OBJ_DATA": "NO"}
    url = HORIZONS + "?" + urllib.parse.urlencode(q)
    t = http_get(url, "horizons/%s.txt" % h(url), pause=1.5)
    if "$$SOE" not in t:
        os.remove(os.path.join(CACHE, "horizons/%s.txt" % h(url)))
        raise RuntimeError("Horizons returned no ephemeris for %s: %s" % (obj["name"], t[-300:].strip()))
    head = [x.strip() for x in t.split("$$SOE")[0].strip().splitlines()[-2].split(",")]
    rows = t.split("$$SOE")[1].split("$$EOE")[0].strip().splitlines()
    ira = next((i for i, c in enumerate(head) if c.startswith("R.A.")), 3)
    isr = next((i for i, c in enumerate(head) if c == "RA_3sigma"), None)
    out = []
    for r in rows:
        p = [x.strip() for x in r.split(",")]
        sra = sde = None
        if isr is not None and len(p) > isr + 1 and fnum(p[isr]) is not None and fnum(p[isr + 1]) is not None:
            sra, sde = float(p[isr]) / 3, float(p[isr + 1]) / 3
        out.append((float(p[ira]), float(p[ira + 1]), sra, sde))
    return out


# ------------------------------------------------------------------ select
def group_sequences(rows, nmin=3, span_min_h=0.5, span_max_h=3.0, keyf=None):
    """Groups exposures of the same quadrant (field/ccd/qid/filter) within a night, in windows of ≤3 h, with ≥nmin."""
    g = {}
    for r in rows:
        try:
            if int(r.get("infobits") or 0) != 0:
                continue
        except ValueError:
            pass
        k = (r["field"], r["ccdid"], r["qid"])
        if keyf and not keyf(r):
            continue
        g.setdefault(k, []).append(r)
    seqs = []
    for k, v in g.items():
        v.sort(key=lambda r: float(r["obsjd"]))
        used = set()
        for i in range(len(v)):
            if i in used:
                continue
            j0 = float(v[i]["obsjd"])
            win = [x for x in range(i, len(v)) if float(v[x]["obsjd"]) - j0 <= span_max_h / 24]
            span = (float(v[win[-1]]["obsjd"]) - j0) * 24
            if len(win) < nmin or span < span_min_h:
                continue
            # pick nmin (up to 5) exposures spread along the arc
            m = min(5, len(win))
            pick = sorted({win[round(t * (len(win) - 1) / (m - 1))] for t in range(m)})
            if len(pick) < nmin:
                continue
            used.update(win)
            seqs.append({"key": k + ("".join(sorted({v[x]["filtercode"][1:] for x in pick})),),
                         "exps": [v[x] for x in pick], "span_h": span})
    return seqs


def load_exclude(paths):
    """Set of (field, night) and ids already used in other field files (e.g. the development set)."""
    ex = set()
    for p in paths or []:
        if not os.path.exists(p):
            continue
        for q in json.load(open(p))["sequences"]:
            ex.add(q["id"])
            ex.add(("fn", int(q["field"]), q["exposures"][0]["filefracday"][:8]))
    return ex


def excluded(s, excl):
    return ("fn", int(s["key"][0]), s["exps"][0]["filefracday"][:8]) in excl


def cmd_select(a):
    rnd = random.Random(a.seed)
    excl = load_exclude(a.exclude)
    # points on the ecliptic (|β| ≤ 8°) visible from Palomar (δ > −25°), 2019–2023
    cand = []
    bins = {}
    tries = 0
    while tries < a.probes and sum(min(len(b), a.per_bin) for b in bins.values()) < a.n * 1.0:
        tries += 1
        lon = rnd.uniform(0, 360)
        lat = rnd.uniform(-8, 8)
        ra, dec = ecl2eq(lon, lat)
        if dec < -25:
            continue
        jd0 = rnd.uniform(2458550, 2460000)
        try:
            rows = irsa_search_fp(ra, dec, jd0, jd0 + 365)  # with the quadrant corners (whole cutout inside it)
        except RuntimeError as e:
            log("IRSA:", e)
            continue
        seqs = [x for x in group_sequences(rows, span_min_h=a.span_min_h) if not excluded(x, excl)]
        log("probe %d: (%.2f,%.2f) %d exposures, %d sequences" % (tries, ra, dec, len(rows), len(seqs)))
        rnd.shuffle(seqs)
        for s in seqs[:4]:
            jds = [float(e["obsjd"]) + float(e.get("exptime") or 30) / 2 / 86400 for e in s["exps"]]
            jm = (jds[0] + jds[-1]) / 2
            try:
                objs = skybot(ra, dec, jm, 0.25)
            except RuntimeError as e:
                log("SkyBoT:", e)
                continue
            # target: asteroid with V 15–20 within 8' of the searched point (cutout center = mean target position)
            good = [o for o in objs if a.vmin <= o["V"] <= a.vmax and sep_arcsec(ra, dec, o["ra"], o["dec"]) < 480
                    and o["err"] < 5]
            if not good:
                continue
            rnd.shuffle(good)
            o = good[0]
            vbin = int(o["V"])
            rate = math.hypot(o["dra_h"], o["ddec_h"]) / 60  # "/min
            if len(bins.get(vbin, [])) >= a.per_bin:
                continue
            # cutout sized by the motion over the arc (as in select-targets) and fully inside the quadrant in all
            # exposures: a cutout touching the edge comes back truncated and the sequence becomes data_error in run
            size = max(CUTOUT_ARCMIN, math.ceil(rate * s["span_h"] + 8))
            if size > 30 or not all(inside_quadrant(e, o["ra"], o["dec"], size / 2 / 60) for e in s["exps"]):
                continue
            others = [x for x in skybot(o["ra"], o["dec"], jm, size / 60 * 0.75)]
            seq = {
                "id": "ztf_%s_f%s_c%s_q%s_%s" % (s["exps"][0]["filefracday"][:8], s["key"][0], s["key"][1],
                                                  s["key"][2], s["key"][3]),
                "ra": round(o["ra"], 6), "dec": round(o["dec"], 6), "jd_mid": round(jm, 6),
                "night": "%04d-%02d-%02d" % jd_to_iso(float(s["exps"][0]["obsjd"]) - 0.3)[:3],
                "field": int(s["key"][0]), "ccdid": int(s["key"][1]), "qid": int(s["key"][2]),
                "filters": s["key"][3], "span_h": round(s["span_h"], 3), "cutout_arcmin": size,
                "exposures": [{"obsjd": float(e["obsjd"]), "filefracday": e["filefracday"],
                               "exptime": float(e.get("exptime") or 30), "filtercode": e["filtercode"], "seeing": e.get("seeing"),
                               "maglimit": e.get("maglimit")} for e in s["exps"]],
                "target": {"num": o["num"], "name": o["name"], "V": o["V"], "rate_arcsec_min": round(rate, 3),
                           "dra_h": o["dra_h"], "ddec_h": o["ddec_h"], "cls": o["cls"]},
                "known": [{"num": x["num"], "name": x["name"], "V": x["V"], "ra": x["ra"], "dec": x["dec"],
                           "dra_h": x["dra_h"], "ddec_h": x["ddec_h"]} for x in others],
            }
            seq["gal_b"] = round(gal_b(seq["ra"], seq["dec"]), 2)
            if any(c["id"] == seq["id"] for c in cand) or seq["id"] in excl:
                continue
            bins.setdefault(vbin, []).append(seq)
            cand.append(seq)
            log("  + %s target %s V=%.1f %.2f\"/min, %d exp over %.2f h" % (seq["id"], o["name"], o["V"], rate,
                                                                       len(seq["exposures"]), s["span_h"]))
    cand.sort(key=lambda s: (s["target"]["V"], s["id"]))
    cand = cand[:a.n]
    doc = {"description": "ZTF (IRSA) sequences for the Piazzi benchmark; generated by `run.py select`.",
           "seed": a.seed, "vrange": [a.vmin, a.vmax], "cutout_arcmin": CUTOUT_ARCMIN, "sequences": cand}
    write_fields(a.out, doc, a.append)


def write_fields(out, doc, append):
    if append and os.path.exists(out):
        old = json.load(open(out))
        ids = {q["id"] for q in old["sequences"]}
        old["sequences"] += [q for q in doc["sequences"] if q["id"] not in ids]
        old.setdefault("parts", []).append({k: v for k, v in doc.items() if k != "sequences"})
        doc = old
    with open(out, "w") as f:
        json.dump(doc, f, indent=1, ensure_ascii=False)
    log("%d sequences in %s" % (len(doc["sequences"]), out))


# ------------------------------------------------------------------ select-targets (known targets: NEAs, bright objects)
CAD_API = "https://ssd-api.jpl.nasa.gov/cad.api"


def gal_b(ra, dec):
    """Galactic latitude (degrees) — only to describe field density."""
    r, d = math.radians(ra), math.radians(dec)
    rp, dp = math.radians(192.85948), math.radians(27.12825)
    return math.degrees(math.asin(math.sin(d) * math.sin(dp) + math.cos(d) * math.cos(dp) * math.cos(r - rp)))


def horizons_ephem(cmd, jd0, jd1, step="20m"):
    """I41 ephemeris: list of (jd, sun_flag, ra, dec, dra*cosd "/h, ddec "/h, V)."""
    q = {"format": "text", "COMMAND": "'%s'" % cmd, "CENTER": "'I41@399'", "MAKE_EPHEM": "YES",
         "EPHEM_TYPE": "OBSERVER", "START_TIME": "'JD%.4f'" % jd0, "STOP_TIME": "'JD%.4f'" % jd1,
         "STEP_SIZE": "'%s'" % step, "QUANTITIES": "'1,3,9'", "ANG_FORMAT": "DEG", "CSV_FORMAT": "YES",
         "OBJ_DATA": "NO", "CAL_FORMAT": "JD", "TIME_TYPE": "UT"}
    url = HORIZONS + "?" + urllib.parse.urlencode(q)
    t = http_get(url, "horizons/%s.txt" % h(url), pause=1.5)
    if "$$SOE" not in t:
        os.remove(os.path.join(CACHE, "horizons/%s.txt" % h(url)))
        raise RuntimeError("Horizons: %s" % t[-200:].strip())
    out = []
    for r in t.split("$$SOE")[1].split("$$EOE")[0].strip().splitlines():
        p = [x.strip() for x in r.split(",")]
        try:
            out.append((float(p[0]), p[1], float(p[3]), float(p[4]), float(p[5]), float(p[6]),
                        float(p[7]) if fnum(p[7]) is not None else 99.0))
        except (ValueError, IndexError):
            pass
    return out


def interp_eph(eph, jd):
    for a_, b_ in zip(eph, eph[1:]):
        if a_[0] <= jd <= b_[0]:
            f = (jd - a_[0]) / (b_[0] - a_[0])
            dra = b_[2] - a_[2]
            if dra > 180:
                dra -= 360
            elif dra < -180:
                dra += 360
            return (a_[2] + f * dra) % 360, a_[3] + f * (b_[3] - a_[3])
    return None


def irsa_search_fp(ra, dec, jd0, jd1):
    q = {"POS": "%.5f,%.5f" % (ra, dec), "ct": "csv",
         "COLUMNS": "obsjd,field,ccdid,qid,filtercode,filefracday,exptime,seeing,maglimit,infobits,"
                    "ra,dec,ra1,dec1,ra2,dec2,ra3,dec3,ra4,dec4",
         "WHERE": "obsjd>%.4f AND obsjd<%.4f" % (jd0, jd1)}
    url = IRSA_SEARCH + "?" + urllib.parse.urlencode(q)
    t = http_get(url, "irsa/%s.csv" % h(url), pause=1.0)
    return list(csv.DictReader(io.StringIO(t)))


def tangent(ra0, de0, ra, de):
    r0, d0, r, d = map(math.radians, (ra0, de0, ra, de))
    c = math.sin(d0) * math.sin(d) + math.cos(d0) * math.cos(d) * math.cos(r - r0)
    x = math.cos(d) * math.sin(r - r0) / c
    y = (math.cos(d0) * math.sin(d) - math.sin(d0) * math.cos(d) * math.cos(r - r0)) / c
    return math.degrees(x), math.degrees(y)


def inside_quadrant(e, ra, dec, margin_deg):
    """Point inside the quadrant (polygon of the 4 corners) with margin_deg clearance from the edges."""
    try:
        c0 = (float(e["ra"]), float(e["dec"]))
        cs = [tangent(c0[0], c0[1], float(e["ra%d" % i]), float(e["dec%d" % i])) for i in range(1, 5)]
    except (KeyError, ValueError, TypeError):
        return False
    px, py = tangent(c0[0], c0[1], ra, dec)
    # shrink the polygon toward the center by the margin
    sh = []
    for x, y in cs:
        r = math.hypot(x, y)
        k = max(0.0, (r - margin_deg * 1.42) / r) if r else 0
        sh.append((x * k, y * k))
    ins = False
    for i in range(4):
        x1, y1 = sh[i]
        x2, y2 = sh[(i + 1) % 4]
        if (y1 > py) != (y2 > py) and px < (x2 - x1) * (py - y1) / (y2 - y1) + x1:
            ins = not ins
    return ins


def cad_objects(rnd, nmax, date_min="2019-03-01", date_max="2023-10-01"):
    q = {"date-min": date_min, "date-max": date_max, "dist-max": "0.08", "h-max": "23", "body": "Earth"}
    url = CAD_API + "?" + urllib.parse.urlencode(q)
    d = json.loads(http_get(url, "cad/%s.json" % h(url), pause=1.0))
    f = d["fields"]
    rows = [dict(zip(f, r)) for r in d["data"]]
    rnd.shuffle(rows)
    seen, out = set(), []
    for r in rows:
        if r["des"] in seen:
            continue
        seen.add(r["des"])
        out.append((r["des"], float(r["jd"])))
    return out[:nmax]


def cmd_select_targets(a):
    rnd = random.Random(a.seed)
    excl = load_exclude((a.exclude or []) + ([a.out] if a.append else []))
    if a.objects == "cad":
        objs = cad_objects(rnd, 400, a.date_min, a.date_max)
    else:
        objs = [(x.strip(), None) for x in a.objects.split(",") if x.strip()]
    # excluded objects: --exclude-objects + targets of the files in --exclude
    norm = lambda x: re.sub(r"[\s()]", "", str(x)).upper()
    bad = {norm(x) for x in (a.exclude_objects or "").split(",") if x.strip()}
    for p in a.exclude or []:
        if os.path.exists(p):
            for q in json.load(open(p))["sequences"]:
                t = q.get("target") or {}
                for k in ("name", "num", "des"):
                    if t.get(k):
                        bad.add(norm(t[k]))
    objs = [o for o in objs if norm(o[0]) not in bad]
    cand = []
    for des, jd_ca in objs:
        if len(cand) >= a.n:
            break
        cmd = ("%s;" % des) if re.match(r"^\d+$", des) else ("DES=%s;" % des)
        got = 0
        for w in range(a.windows):
            if got >= a.per_object or len(cand) >= a.n:
                break
            # 6-day window: around the close approach (NEA) or random within 2019-03..2023-09
            if jd_ca:
                j0 = jd_ca + rnd.uniform(-12, 6) - 3
            else:
                j0 = rnd.uniform(2458550, 2460200)
            try:
                eph = horizons_ephem(cmd, j0, j0 + 6, "20m")
            except RuntimeError as e:
                log("Horizons %s: %s" % (des, str(e)[-100:]))
                break
            ok = [r for r in eph if r[1] in ("", "A", "N") and r[3] > -25 and a.vmin <= r[6] <= a.vmax
                  and a.rate_min <= math.hypot(r[4], r[5]) / 60 <= a.rate_max]
            log("%s window %d: %d/%d usable points" % (des, w, len(ok), len(eph)))
            nights = {}
            for r in ok:
                nights.setdefault(int(r[0] - 0.25), []).append(r)
            for nk in rnd.sample(sorted(nights), len(nights)):
                if got >= a.per_object:
                    break
                rows_n = nights[nk]
                for r in [rows_n[len(rows_n) // 2]] + ([rows_n[0], rows_n[-1]] if len(rows_n) > 6 else []):
                    rate = math.hypot(r[4], r[5]) / 60
                    try:
                        rows = irsa_search_fp(r[2], r[3], r[0] - 0.08, r[0] + 0.08) if rate > 4 else \
                            irsa_search_fp(r[2], r[3], r[0] - 0.15, r[0] + 0.15)
                    except RuntimeError as e:
                        log("IRSA:", str(e)[-100:])
                        continue

                    def keyf(e):
                        p = interp_eph(eph, float(e["obsjd"]))
                        return p is not None and inside_quadrant(e, p[0], p[1], 0.06)
                    span_max = 3.0 if rate < 5 else max(0.4, 15.0 / rate)  # arc ≤ ~15'
                    seqs = [x for x in group_sequences(rows, span_min_h=a.span_min_h, span_max_h=min(3.0, span_max),
                                                       keyf=keyf) if not excluded(x, excl)]
                    log("  %s night %d: %.2f\"/min V=%.1f, %d exposures, %d sequences" % (
                        des, nk, rate, r[6], len(rows), len(seqs)))
                    if not seqs:
                        continue
                    s = max(seqs, key=lambda x: (len(x["exps"]), x["span_h"]))
                    seq = make_target_seq(s, eph, des, a)
                    if seq is None or seq["id"] in excl or any(c["id"] == seq["id"] for c in cand):
                        continue
                    cand.append(seq)
                    excl.add(("fn", seq["field"], seq["exposures"][0]["filefracday"][:8]))
                    got += 1
                    log("  + %s target %s (%s) V=%.1f %.2f\"/min, %d exp over %.2f h, cutout %g'" % (
                        seq["id"], seq["target"]["name"], seq["target"]["cls"], seq["target"]["V"],
                        seq["target"]["rate_arcsec_min"], len(seq["exposures"]), seq["span_h"], seq["cutout_arcmin"]))
                    break
    doc = {"description": "ZTF sequences from known targets (Horizons); generated by `run.py select-targets`.",
           "label": a.label, "seed": a.seed, "sequences": cand}
    write_fields(a.out, doc, a.append)


def make_target_seq(s, eph, des, a):
    exps = s["exps"]
    jds = [float(e["obsjd"]) + float(e.get("exptime") or 30) / 2 / 86400 for e in exps]
    jm = (jds[0] + jds[-1]) / 2
    p0, p1, pm = interp_eph(eph, jds[0]), interp_eph(eph, jds[-1]), interp_eph(eph, jm)
    if not (p0 and p1 and pm):
        return None
    motion = sep_arcsec(p0[0], p0[1], p1[0], p1[1]) / 60
    size = max(CUTOUT_ARCMIN, math.ceil(motion + 8))
    if size > 30:
        return None
    # the whole cutout (not just the target) inside the quadrant in all exposures: otherwise it comes back truncated (data_error)
    if not all(inside_quadrant(e, pm[0], pm[1], size / 2 / 60) for e in exps):
        return None
    # center = arc midpoint (not the position at jm, which for curved paths is nearly the same)
    cra, cde = pm
    try:
        objs = skybot(cra, cde, jm, size / 60 * 0.75 + motion / 120)
    except RuntimeError as e:
        log("SkyBoT:", e)
        return None
    num = des if re.match(r"^\d+$", des) else None
    tgt = None
    for o in objs:
        if (num and o["num"] == num) or o["name"].replace(" ", "") == des.replace(" ", ""):
            tgt = o
    i = min(range(len(eph)), key=lambda k: abs(eph[k][0] - jm))
    if tgt is None:
        log("  target %s not in SkyBoT; using Horizons" % des)
        tgt = {"num": num or "-", "name": des, "V": eph[i][6], "cls": "?", "ra": cra, "dec": cde,
               "dra_h": eph[i][4], "ddec_h": eph[i][5]}
        objs.append(tgt)
    # the target's reference position/rate comes from Horizons at time jm (SkyBoT is linear and less accurate
    # for fast NEAs); association in `run` uses ±10"
    tgt = dict(tgt, ra=cra, dec=cde, dra_h=eph[i][4], ddec_h=eph[i][5])
    objs = [tgt if (o["name"] == tgt["name"]) else o for o in objs]
    rate = math.hypot(tgt["dra_h"], tgt["ddec_h"]) / 60
    k = s["key"]
    return {
        "id": "ztf_%s_f%s_c%s_q%s_%s" % (exps[0]["filefracday"][:8], k[0], k[1], k[2], k[3]),
        "ra": round(cra, 6), "dec": round(cde, 6), "jd_mid": round(jm, 6),
        "night": "%04d-%02d-%02d" % jd_to_iso(float(exps[0]["obsjd"]) - 0.3)[:3],
        "field": int(k[0]), "ccdid": int(k[1]), "qid": int(k[2]), "filters": k[3], "span_h": round(s["span_h"], 3),
        "cutout_arcmin": size, "gal_b": round(gal_b(cra, cde), 2), "motion_arcmin": round(motion, 2),
        "exposures": [{"obsjd": float(e["obsjd"]), "filefracday": e["filefracday"],
                       "exptime": float(e.get("exptime") or 30), "filtercode": e["filtercode"], "seeing": e.get("seeing"),
                       "maglimit": e.get("maglimit")} for e in exps],
        "target": {"num": tgt["num"], "name": tgt["name"], "V": tgt["V"], "rate_arcsec_min": round(rate, 3),
                   "dra_h": tgt["dra_h"], "ddec_h": tgt["ddec_h"], "cls": tgt["cls"], "V_horizons": eph[i][6],
                   "source": "horizons"},
        "known": [{"num": x["num"], "name": x["name"], "V": x["V"], "ra": x["ra"], "dec": x["dec"],
                   "dra_h": x["dra_h"], "ddec_h": x["ddec_h"]} for x in objs],
    }


# ------------------------------------------------------------------ SDSS (independent test set)
# Drift-scan camera of the 2.5 m Sloan telescope, Apache Point (MPC 645): each field (run/camcol/field)
# is imaged through the filters in the order r, i, u, z, g, ~71.7 s between adjacent filters. Calibrated DR17 frames (nanomaggies)
# are converted to electrons (img/calib + sky)·gain, unsigned integers (BITPIX 16, BZERO 32768), GAIN = 1.
# Scale 0.396″/px: pixel 0.024 mm and focal length 0.024/0.396·206265 = 12501 mm. DATE-OBS written = UTC of MID-exposure
# at the target's row in that filter -> time_reference=1 (mid). See tools/benchmark/README.md.
SDSS_FRAMES = "https://data.sdss.org/sas/dr17/eboss/photoObj/frames/301"
SKYSERVER_SQL = "https://skyserver.sdss.org/dr17/SkyServerWS/SearchTools/SqlSearch"
FIELDS_SDSS = os.path.join(HERE, "fields_sdss.json")
SDSS_SETTINGS = {
    "focal_length_mm": "12501", "pixel_width_mm": "0.024", "pixel_height_mm": "0.024", "mpc_code": "645", "time_reference": "1",
    "catalog_faint_mag": "20.0", "saturation": "60000", "observer": "SDSS", "measurer": "Piazzi benchmark",
    "submitter": "Piazzi benchmark", "telescope": "2.5-m f/5 Ritchey-Chretien + CCD (drift scan)",
}
SDSS_ORDER = "riuzg"           # temporal order of the filters in the camera
SDSS_FILTERS = "rizg"          # filters used by default (u is shallow and much fainter for asteroids)
SDSS_NROW = 2048               # rows a sky point takes to cross the CCD (integration time = EXPTIME)
# approximate gain (e-/DN) per filter and camcol (DR7 table; affects only σ/SNR, not positions)
SDSS_GAIN = {"u": [1.62, 1.71, 1.59, 1.60, 1.47, 2.17], "g": [3.32, 3.855, 3.845, 3.995, 4.05, 4.035],
             "r": [4.71, 4.60, 4.72, 4.76, 4.725, 4.895], "i": [5.165, 6.565, 4.86, 4.885, 4.64, 4.76],
             "z": [4.745, 5.155, 4.885, 4.775, 3.48, 4.69]}


def tai_minus_utc(mjd):
    return 32.0 if mjd < 53736 else (33.0 if mjd < 54832 else 34.0)  # 2006-01-01, 2009-01-01


def tai_to_jd_utc(tai_s):
    mjd = tai_s / 86400.0
    return mjd - tai_minus_utc(mjd) / 86400.0 + 2400000.5


def skyserver(sql):
    url = SKYSERVER_SQL + "?" + urllib.parse.urlencode({"cmd": sql, "format": "csv"})
    t = http_get(url, "skyserver/%s.csv" % h(url), pause=1.0, timeout=300)
    t = "\n".join(l for l in t.splitlines() if not l.startswith("#"))
    return list(csv.DictReader(io.StringIO(t)))


def fits_hdus(d):
    """List of (header, data_offset, size) for all HDUs of an in-memory FITS."""
    out, pos = [], 0
    while pos + 2880 <= len(d):
        hd, k = {}, pos
        while True:
            c = d[k:k + 80].decode("ascii", "replace")
            k += 80
            if c.startswith("END") and c[3:].strip() == "":
                break
            if c[8:10] == "= ":
                v = c[10:]
                v = v[1:].split("'")[0].strip() if v.strip().startswith("'") else v.split("/")[0].strip()
                hd[c[:8].strip()] = v
        data = (k + 2879) // 2880 * 2880
        n = abs(int(hd.get("BITPIX", 8))) // 8
        na = int(hd.get("NAXIS", 0))
        size = 0
        if na:
            size = n
            for i in range(1, na + 1):
                size *= int(hd["NAXIS%d" % i])
        size += int(hd.get("PCOUNT", 0))
        out.append((hd, data, size))
        pos = data + (size + 2879) // 2880 * 2880
    return out


def be_floats(b):
    import array
    a = array.array("f")
    a.frombytes(b)
    if sys.byteorder == "little":
        a.byteswap()
    return a


class SdssWcs:
    def __init__(self, hd):
        self.ra0, self.de0 = float(hd["CRVAL1"]), float(hd["CRVAL2"])
        self.x0, self.y0 = float(hd["CRPIX1"]), float(hd["CRPIX2"])
        self.cd = [[float(hd["CD1_1"]), float(hd["CD1_2"])], [float(hd["CD2_1"]), float(hd["CD2_2"])]]

    def sky2pix(self, ra, de):
        """(x, y) 0-based (column, row) via the header TAN projection."""
        a, d, a0, d0 = map(math.radians, (ra, de, self.ra0, self.de0))
        cc = math.sin(d) * math.sin(d0) + math.cos(d) * math.cos(d0) * math.cos(a - a0)
        xi = math.degrees(math.cos(d) * math.sin(a - a0) / cc)
        eta = math.degrees((math.sin(d) * math.cos(d0) - math.cos(d) * math.sin(d0) * math.cos(a - a0)) / cc)
        (p, q), (r, s) = self.cd
        det = p * s - q * r
        dx, dy = (s * xi - q * eta) / det, (-r * xi + p * eta) / det
        return self.x0 + dx - 1, self.y0 + dy - 1


def sdss_frame_url(f, run, camcol, field):
    return "%s/%d/%d/frame-%s-%06d-%d-%04d.fits.bz2" % (SDSS_FRAMES, run, camcol, f, run, camcol, field)


def sdss_frame(f, run, camcol, field):
    """Decompressed DR17 frame (bytes), with the .bz2 cached in cache/sdss/raw."""
    import bz2
    name = "sdss/raw/frame-%s-%06d-%d-%04d.fits.bz2" % (f, run, camcol, field)
    d = http_get(sdss_frame_url(f, run, camcol, field), name, binary=True, pause=0.5, timeout=300)
    return bz2.decompress(d)


def sdss_frame_info(d):
    hdus = fits_hdus(d)
    hd = hdus[0][0]
    return {"tai_row0": float(hd["TAI"]), "exptime": float(hd["EXPTIME"]), "nx": int(hd["NAXIS1"]),
            "ny": int(hd["NAXIS2"]), "wcs": SdssWcs(hd), "hdr": hd, "hdus": hdus}


def sdss_row_tai(info, row):
    """TAI (s) of MID-exposure for a point at row `row` (0-based). The header TAI is that of the 1st row at the end
    of integration; the per-object photoObj times (TAI_<f>) follow TAI_hdr − EXPTIME/2 + row·EXPTIME/2048
    (checked against SkyServer: residual < 0.1 s)."""
    e = info["exptime"]
    return info["tai_row0"] - e / 2 + row * e / SDSS_NROW


def sdss_convert(d, info, f, camcol, jd_mid, out):
    """Frame (nanomaggies) -> electrons, uint16 (BITPIX 16 + BZERO), with DATE-OBS = UTC of mid-exposure."""
    import array
    hdus = info["hdus"]
    hd0, off0, _ = hdus[0]
    nx, ny = info["nx"], info["ny"]
    img = be_floats(d[off0:off0 + 4 * nx * ny])
    h1, off1, _ = hdus[1]
    calib = be_floats(d[off1:off1 + 4 * nx])
    # sky: HDU2 (binary table), ALLSKY column (coarse grid in DN); we use the median (constant)
    h2, off2, _ = hdus[2]
    sky = 0.0
    try:
        tf = h2["TFORM1"]
        nsky = int(re.match(r"(\d+)E", tf)[1])
        v = sorted(x for x in be_floats(d[off2:off2 + 4 * nsky]) if x > 0)
        sky = v[len(v) // 2] if v else 0.0
    except (KeyError, TypeError):
        pass
    gain = SDSS_GAIN[f][camcol - 1]
    inv = [gain / c if c > 0 else 0.0 for c in calib]
    skye = sky * gain
    o = array.array("h", bytes(2 * nx * ny))
    for y in range(ny):
        b = y * nx
        for x in range(nx):
            v = img[b + x] * inv[x] + skye
            v = 0 if v < 0 else (65535 if v > 65535 else int(v + 0.5))
            o[b + x] = v - 32768
    if sys.byteorder == "little":
        o.byteswap()
    yy, mo, dd, hh = jd_to_iso(jd_mid)
    sec = hh * 3600
    ts = "%04d-%02d-%02dT%02d:%02d:%06.3f" % (yy, mo, dd, int(sec // 3600), int(sec % 3600 // 60), sec % 60)
    cards = [("SIMPLE", "T", ""), ("BITPIX", "16", ""), ("NAXIS", "2", ""), ("NAXIS1", str(nx), ""),
             ("NAXIS2", str(ny), ""), ("BZERO", "32768", ""), ("BSCALE", "1", ""),
             ("DATE-OBS", "'%s'" % ts, "UTC mid-exposure at target row"),
             ("EXPTIME", "%.6f" % info["exptime"], "s"), ("GAIN", "1.0", "data in electrons"),
             ("FILTER", "'%s'" % f, ""), ("SKYE", "%.2f" % skye, "sky added (e-)"),
             ("TAIROW0", "%.3f" % info["tai_row0"], "TAI from original header")]
    for k in ("RUN", "CAMCOL", "CTYPE1", "CTYPE2", "CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2", "CD1_1", "CD1_2",
              "CD2_1", "CD2_2"):
        if k in hd0:
            v = hd0[k]
            cards.append((k, "'%s'" % v if k.startswith("CTYPE") else v, ""))
    hdr = "".join(("%-8s= %20s / %s" % (k, v, c) if c else "%-8s= %20s" % (k, v)).ljust(80)[:80] for k, v, c in cards)
    hdr += "END".ljust(80)
    hdr = hdr.ljust((len(hdr) + 2879) // 2880 * 2880)
    data = o.tobytes()
    data += b"\0" * ((-len(data)) % 2880)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out + ".tmp", "wb") as fh:
        fh.write(hdr.encode("ascii") + data)
    os.replace(out + ".tmp", out)


def sdss_paths(seq):
    return [os.path.join(CACHE, "fits", seq["id"], "%d_%s.fits" % (k, e["filter"])) for k, e in enumerate(seq["exposures"])]


def fetch_sdss(seq):
    """Converts (or finds in cache) the SDSS sequence frames; returns the paths in temporal order."""
    out = []
    for p, e in zip(sdss_paths(seq), seq["exposures"]):
        if not (os.path.exists(p) and os.path.getsize(p) > 0):
            d = sdss_frame(e["filter"], seq["run"], seq["camcol"], seq["field"])
            sdss_convert(d, sdss_frame_info(d), e["filter"], seq["camcol"], e["jd_mid"], p)
        out.append(p)
    return out


def sdss_target_times(tgt, infos, filters):
    """Target times (JD UTC, mid) and pixels in each filter: iterates time <-> row (the time depends on the row)."""
    rows = {f: infos[f]["ny"] / 2 for f in filters}
    for _ in range(3):
        jds = [tai_to_jd_utc(sdss_row_tai(infos[f], rows[f])) for f in filters]
        eph = horizons(tgt, jds, site="645")
        pix = {f: infos[f]["wcs"].sky2pix(eph[k][0], eph[k][1]) for k, f in enumerate(filters)}
        rows = {f: pix[f][1] for f in filters}
    return jds, eph, pix


def cmd_select_sdss(a):
    rnd = random.Random(a.seed)
    filters = "".join(f for f in SDSS_ORDER if f in a.filters)
    fl = skyserver("select run,camcol,field,ra,dec,mjd_r,quality from Field where field %% %d = %d and quality >= 2"
                   % (a.modulus, a.residue))
    log("SkyServer: %d fields" % len(fl))
    cand = []
    for r in fl:
        ra, de = float(r["ra"]), float(r["dec"])
        # ecliptic latitude
        e = math.radians(23.4393)
        rr, dd = math.radians(ra), math.radians(de)
        beta = math.degrees(math.asin(math.sin(dd) * math.cos(e) - math.cos(dd) * math.sin(e) * math.sin(rr)))
        if abs(beta) <= a.beta_max:
            cand.append(dict(r, beta=beta))
    rnd.shuffle(cand)
    log("%d fields with |beta| <= %g" % (len(cand), a.beta_max))
    vbins = {}
    rbins = {}
    rate_edges = [0.3, 0.6, 1.0]
    rbin = lambda x: sum(1 for t in rate_edges if x >= t)
    seqs = []
    used_runs = {}
    for fr in cand[:a.probes]:
        if len(seqs) >= a.n:
            break
        run, camcol, field = int(fr["run"]), int(fr["camcol"]), int(fr["field"])
        if used_runs.get(run, 0) >= a.per_run:
            continue
        jd_r = float(fr["mjd_r"]) + 2400000.5 + 26.96 / 86400 - tai_minus_utc(float(fr["mjd_r"])) / 86400
        try:
            objs = skybot(float(fr["ra"]), float(fr["dec"]), jd_r, 0.11, loc="645")
        except RuntimeError as e:
            log("SkyBoT:", str(e)[-100:])
            continue
        objs = [o for o in objs if o["V"] < 30]
        pool = [o for o in objs if a.vmin <= o["V"] < a.vmax and vbins.get(int(o["V"]), 0) < a.per_bin]
        if not pool:
            continue
        # prioritize the least-filled V bin and rate bin; then the fastest
        rate = lambda o: math.hypot(o["dra_h"], o["ddec_h"]) / 60
        pool.sort(key=lambda o: (vbins.get(int(o["V"]), 0) + rbins.get(rbin(rate(o)), 0), -rate(o)))
        try:
            infos = {}
            for f in filters:
                infos[f] = sdss_frame_info(sdss_frame(f, run, camcol, field))
        except (RuntimeError, OSError, KeyError, ValueError) as e:
            log("frames %d/%d/%d: %s" % (run, camcol, field, str(e)[-100:]))
            continue
        got = None
        for tgt in pool[:3]:
            try:
                jds, eph, pix = sdss_target_times(tgt, infos, filters)
            except RuntimeError as e:
                log("  Horizons %s: %s" % (tgt["name"], str(e)[-100:]))
                continue
            m = a.margin
            if all(m <= pix[f][0] <= infos[f]["nx"] - 1 - m and m <= pix[f][1] <= infos[f]["ny"] - 1 - m
                   for f in filters):
                got = (tgt, jds, eph, pix)
                break
        if not got:
            continue
        tgt, jds, eph, pix = got
        r_ = rate(tgt)
        motion = sep_arcsec(eph[0][0], eph[0][1], eph[-1][0], eph[-1][1])
        jm = (jds[0] + jds[-1]) / 2
        cra = statistics.mean(x[0] for x in eph)
        cde = statistics.mean(x[1] for x in eph)
        jref = jd_r
        s = {
            "id": "sdss_%06d_%d_%04d" % (run, camcol, field), "source": "sdss", "mpc": "645",
            "run": run, "camcol": camcol, "field": field, "filters": filters,
            "ra": round(cra, 6), "dec": round(cde, 6), "jd_mid": round(jref, 7), "beta": round(fr["beta"], 2),
            "night": "%04d-%02d-%02d" % jd_to_iso(jd_r - 0.3)[:3], "gal_b": round(gal_b(cra, cde), 2),
            "motion_arcsec": round(motion, 2), "span_s": round((jds[-1] - jds[0]) * 86400, 1),
            "exposures": [{"filter": f, "jd_mid": round(jds[k], 8), "row": round(pix[f][1], 1),
                           "col": round(pix[f][0], 1), "tai_row0": infos[f]["tai_row0"],
                           "exptime": infos[f]["exptime"]} for k, f in enumerate(filters)],
            "target": {"num": tgt["num"], "name": tgt["name"], "V": tgt["V"], "rate_arcsec_min": round(r_, 3),
                       "dra_h": tgt["dra_h"], "ddec_h": tgt["ddec_h"], "cls": tgt["cls"], "source": "skybot"},
            "known": [{"num": x["num"], "name": x["name"], "V": x["V"], "ra": x["ra"], "dec": x["dec"],
                       "dra_h": x["dra_h"], "ddec_h": x["ddec_h"]} for x in objs],
        }
        try:
            fetch_sdss(s)
        except (RuntimeError, OSError, ValueError) as e:
            log("conversion %s: %s" % (s["id"], e))
            continue
        seqs.append(s)
        used_runs[run] = used_runs.get(run, 0) + 1
        vbins[int(tgt["V"])] = vbins.get(int(tgt["V"]), 0) + 1
        rbins[rbin(r_)] = rbins.get(rbin(r_), 0) + 1
        log("[%d] %s %s V=%.1f %.2f\"/min (%.1f\" in %.0f s) b=%.0f" % (
            len(seqs), s["id"], tgt["name"], tgt["V"], r_, motion, s["span_s"], s["gal_b"]))
        write_fields(a.out, {"source": "sdss", "filters": filters, "seed": a.seed,
                             "created": time.strftime("%Y-%m-%d"), "settings": SDSS_SETTINGS,
                             "sequences": seqs}, False)
    log("V:", dict(sorted(vbins.items())), "rate (<0.3/0.3–0.6/0.6–1/>1 \"/min):", dict(sorted(rbins.items())))


# ------------------------------------------------------------------ run
def fits_header(path):
    d = open(path, "rb").read()
    hd = {}
    for i in range(0, min(len(d), 2880 * 20), 80):
        c = d[i:i + 80].decode("ascii", "replace")
        if c.startswith("END "):
            break
        if c[8:10] == "= ":
            v = c[10:].split("/")[0].strip().strip("'").strip()
            hd[c[:8].strip()] = v
    return hd


def parse_ades(path):
    if not os.path.exists(path):
        return []
    lines = [l for l in open(path, encoding="utf-8", errors="replace") if l.strip() and not l.startswith("#")]
    hdr = None
    out = []
    for l in lines:
        p = [x.strip() for x in l.rstrip("\n").split("|")]
        if hdr is None:
            if "obsTime" in p:
                hdr = p
            continue
        out.append(dict(zip(hdr, p)))
    return out


def parse_log(text):
    plates = []
    for m in re.finditer(r"image (\d+): (\d+) detections \| (\d+) of (\d+) refs, dRA ([-\d.]+)\" dDe ([-\d.]+)\"", text):
        plates.append({"img": int(m[1]), "ndet": int(m[2]), "nref": int(m[3]), "nmatch": int(m[4]),
                       "dra": float(m[5]), "dde": float(m[6])})
    movers = [(float(m[1]), float(m[2])) for m in re.finditer(r"moving object: ([\d.]+)\"/min PA ([-\d.]+)", text)]
    return plates, movers


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def predicted(obj, jd, jd_ref):
    dt = (jd - jd_ref) * 24
    dec = obj["dec"] + obj["ddec_h"] * dt / 3600
    ra = obj["ra"] + obj["dra_h"] * dt / 3600 / math.cos(math.radians(obj["dec"]))
    return ra, dec


EXTRA_SETTINGS = {}  # --set K=V from `run` (experiments), applied on top of the set's settings


def seq_settings(seq):
    return dict(SDSS_SETTINGS if seq.get("source") == "sdss" else ZTF_SETTINGS, **EXTRA_SETTINGS)


def fetch_cutouts(seq, tries=12):
    """Downloads (or finds in cache) the sequence cutouts; returns the paths in exposure order."""
    if seq.get("source") == "sdss":
        return fetch_sdss(seq)
    out = []
    for k, e in enumerate(seq["exposures"]):
        ex = dict(e, field=seq["field"], ccdid=seq["ccdid"], qid=seq["qid"])
        url = cutout_url(ex, seq["ra"], seq["dec"], seq.get("cutout_arcmin", CUTOUT_ARCMIN))
        p = os.path.join(CACHE, "fits", seq["id"], "%d_%s.fits" % (k, e["filefracday"]))
        if not (os.path.exists(p) and os.path.getsize(p) > 0):
            d = http_get(url, pause=1.0, binary=True, timeout=180, tries=tries)
            if not d.startswith(b"SIMPLE"):
                raise RuntimeError("invalid cutout: %r" % d[:200])
            os.makedirs(os.path.dirname(p), exist_ok=True)
            open(p, "wb").write(d)
        out.append(p)
    return out


def cmd_fetch(a):
    """Pre-downloads all cutouts, in several passes (the IRSA cutout service sometimes returns 503)."""
    seqs = json.load(open(a.fields))["sequences"]
    done = set()
    for ps in range(a.passes):
        miss = 0
        for s in seqs:
            if s["id"] in done:
                continue
            try:
                fetch_cutouts(s, tries=3)
                if a.piazzi:
                    _vizier_down.clear()
                    prepare_catalog(os.path.abspath(a.piazzi), *seq_inputs(s))
                done.add(s["id"])
            except (RuntimeError, subprocess.TimeoutExpired) as e:
                miss += 1
                log("passada %d: %s: %s" % (ps + 1, s["id"], str(e)[-100:]))
        log("pass %d: %d incomplete sequences" % (ps + 1, miss))
        if miss:
            time.sleep(a.sleep)
        if not miss:
            break


_vizier_down = []
NOIRLAB_TAP = "https://datalab.noirlab.edu/tap/sync"
ESA_TAP = "https://gea.esac.esa.int/tap-server/tap/sync"
# TAP services with the Gaia DR3 table (same columns), tried in this order when VizieR fails
GAIA_TAPS = [(NOIRLAB_TAP, "gaia_dr3.gaia_source"), (ESA_TAP, "gaiadr3.gaia_source")]


def gaia_catalog(caturl):
    """Gaia DR3 catalog for the box: first the VizieR URL printed by piazzi; if VizieR fails, the same box
    via the NOIRLab Astro Data Lab TAP (gaia_dr3.gaia_source), converted to VizieR TSV (RAJ2000/DEJ2000 = position
    propagated linearly from 2016.0 to 2000.0 with proper motion, as VizieR does)."""
    catp = os.path.join(CACHE, "gaia", h(caturl) + ".tsv")
    if os.path.exists(catp) and os.path.getsize(catp) > 0:
        return catp, open(catp + ".src").read().strip() if os.path.exists(catp + ".src") else "cache"
    os.makedirs(os.path.dirname(catp), exist_ok=True)
    try:
        if _vizier_down:
            raise RuntimeError("VizieR down during this run")
        t = http_get(caturl, pause=2.0, tries=1, timeout=60)
        if t.count("\n") > 50 and "RAJ2000" in t:
            open(catp, "w").write(t)
            open(catp + ".src", "w").write("vizier")
            return catp, "vizier"
    except RuntimeError as e:
        _vizier_down.append(1)
        log("   VizieR failed (%s); using NOIRLab TAP" % str(e)[-80:])
    q = urllib.parse.parse_qs(urllib.parse.urlparse(caturl).query)
    ra0, de0 = [float(x) for x in q["-c"][0].split()]
    r = float(q["-c.r"][0]) / 60  # VizieR -c.geom=b: box; we use ±r (superset)
    rp = float(q["RPmag"][0].lstrip("<"))
    dra = r / max(0.05, math.cos(math.radians(de0)))
    rows, src, err = [], "", None
    for tap, table in GAIA_TAPS:
        adql = ("SELECT ra,dec,pmra,pmdec,phot_g_mean_mag,phot_bp_mean_mag,phot_rp_mean_mag,parallax,"
                "ra_error,dec_error,pmra_error,pmdec_error,ruwe "
                "FROM %s WHERE ra BETWEEN %.6f AND %.6f AND dec BETWEEN %.6f AND %.6f "
                "AND phot_rp_mean_mag < %.2f" % (table, ra0 - dra, ra0 + dra, de0 - r, de0 + r, rp))
        url = tap + "?" + urllib.parse.urlencode({"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv",
                                                  "MAXREC": "500000", "QUERY": adql})
        try:
            t = http_get(url, pause=2.0, timeout=300, tries=2)
            rows = list(csv.DictReader(io.StringIO(t)))
            if len(rows) < 20:
                raise RuntimeError("Gaia catalog empty/short (%d stars): %s" % (len(rows), t[:200]))
            src = "noirlab" if tap == NOIRLAB_TAP else "esa"
            break
        except RuntimeError as e:
            err = e
            log("   TAP %s failed (%s)" % (tap.split("/")[2], str(e)[-80:]))
    if not src:
        raise RuntimeError("no Gaia service responded: %s" % err)
    L = ["# Gaia DR3 via TAP (%s): %s" % (src, adql),
         "RAJ2000\tDEJ2000\tpmRA\tpmDE\tGmag\tBPmag\tRPmag\tPlx\te_RA_ICRS\te_DE_ICRS\te_pmRA\te_pmDE\tRUWE",
         "\t".join(["------"] * 13), "deg\tdeg\tmas/yr\tmas/yr\tmag\tmag\tmag\tmas\tmas\tmas\tmas/yr\tmas/yr\t"]
    for x in rows:
        ra, de = float(x["ra"]), float(x["dec"])
        pa, pd = fnum(x["pmra"]) or 0.0, fnum(x["pmdec"]) or 0.0
        de2 = de - pd * 16 / 3.6e6
        ra2 = ra - pa * 16 / 3.6e6 / math.cos(math.radians(de))
        f = lambda k, n=4: ("%.*f" % (n, float(x[k]))) if fnum(x[k]) is not None else ""
        L.append("%.8f\t%+.8f\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (
            ra2, de2, f("pmra", 3), f("pmdec", 3), f("phot_g_mean_mag"), f("phot_bp_mean_mag"), f("phot_rp_mean_mag"),
            f("parallax"), f("ra_error"), f("dec_error"), f("pmra_error"), f("pmdec_error"), f("ruwe", 3)))
    open(catp, "w").write("\n".join(L) + "\n")
    open(catp + ".src", "w").write(src)
    return catp, src


def prepare_catalog(piazzi, sets, center, imgs):
    """Runs piazzi without a catalog (it only prints the box VizieR URL) and downloads/finds in cache the Gaia catalog."""
    pr = subprocess.run([piazzi] + sets + center + imgs, capture_output=True, text=True, timeout=300)
    m = re.search(r"curl -o \S+ '([^']+)'", pr.stdout)
    if not m:
        raise RuntimeError("piazzi did not print the catalog URL: %s" % (pr.stdout + pr.stderr)[-300:])
    return gaia_catalog(m[1])


def seq_inputs(seq):
    imgs = fetch_cutouts(seq, tries=3)
    sets = []
    for k, v in seq_settings(seq).items():
        sets += ["--set", "%s=%s" % (k, v)]
    return sets, ["--center-deg", "%.6f" % seq["ra"], "%.6f" % seq["dec"]], imgs


def ades_tracks(path):
    """Tracks from piazzi's ADESReport.psv, in order. piazzi writes tracks one after another in time order:
    a new track starts when obsTime does not advance or when trkSub/provID changes."""
    tracks = {}
    prev_t, prev_k, n = None, None, 0
    for r in parse_ades(path):
        k = (r.get("trkSub") or "") + "/" + (r.get("provID") or r.get("permID") or "")
        if prev_t is None or r["obsTime"] <= prev_t or k != prev_k:
            n += 1
        prev_t, prev_k = r["obsTime"], k
        tracks.setdefault("%s#%d" % (k, n), []).append(r)
    return tracks


# Radius of the "whole field" SkyBoT query (degrees): the entire cutout/frame, wherever the target is in it.
# SDSS: 2048×1489 px frame at 0.396″/px (13.5′×9.8′, half-diagonal 8.4′) with the target anywhere -> 2 half-diagonals;
# ZTF: CUTOUT_ARCMIN cutout centered on the target -> half-diagonal.
FIELD_RADIUS_DEG = {"sdss": 0.30, "ztf": CUTOUT_ARCMIN / 60 * 0.75}
# Association radius to a SkyBoT object: the SkyBoT ephemeris carries the position uncertainty (ErrPos column, ″);
# up to 3σ of it is accepted, with floor MATCH_ARCSEC (measurement errors + linear ephemeris interpolation) and cap
# FIELD_MATCH_MAX (beyond that, association becomes chance in asteroid-dense fields).
FIELD_MATCH_MAX = 20.0


def field_objects(seq):
    """All SkyBoT objects in the sequence's whole field, at the central time (position + linear rate)."""
    src = "sdss" if seq.get("source") == "sdss" else "ztf"
    try:
        return skybot(seq["ra"], seq["dec"], seq["jd_mid"], FIELD_RADIUS_DEG[src], loc=seq.get("mpc", "I41"))
    except RuntimeError as e:
        log("   field SkyBoT query failed: %s" % str(e)[:200])
        return []


def field_match_radius(o):
    return max(MATCH_ARCSEC, min(3 * (o.get("err") or 0), FIELD_MATCH_MAX))


def match_field(rr, objs, jref):
    """(name, median separation ″) of the field object closest to track rr, if within its radius; otherwise None."""
    best = None
    for o in objs:
        ds = []
        for r in rr:
            pra, pde = predicted(o, iso_to_jd(r["obsTime"]), jref)
            ds.append(sep_arcsec(float(r["ra"]), float(r["dec"]), pra, pde))
        d = statistics.median(ds)
        if d < field_match_radius(o) and (best is None or d < best[1]):
            best = (o["name"], round(d, 2))
    return best


# Ephemeris uncertainty limit (1σ per axis, ″) below which the Horizons truth is treated as exact for the "clean"
# metrics (numbered objects/long arcs): Piazzi's error dominates the O−C.
SIG_EPH_CLEAN = 0.05


def score_track(rr, o, site):
    """O−C of each ADES measurement of track rr against Horizons (ephemeris position and uncertainty): dra/ddec/abs (″), ADES
    σ, ratio |O−C|/σ, ephemeris σ per axis (1σ) and the corrected ratio |O−C|/√(σ² + σ_eph²)."""
    jds = [iso_to_jd(r["obsTime"]) for r in rr]
    eph = horizons(o, jds, site=site)
    meas = []
    for r, jd, ep in zip(rr, jds, eph):
        ra, de = float(r["ra"]), float(r["dec"])
        dra = (ra - ep[0]) * 3600 * math.cos(math.radians(de))
        if dra > 648000:
            dra -= 1296000 * math.cos(math.radians(de))
        dde = (de - ep[1]) * 3600
        sra, sde = fnum(r.get("rmsRA")), fnum(r.get("rmsDec"))
        tot = math.hypot(dra, dde)
        sig = math.hypot(sra, sde) if sra and sde else None
        era, ede = ep[2], ep[3]
        seph = math.hypot(era, ede) if era is not None and ede is not None else None
        m = {"obsTime": r["obsTime"], "dra": round(dra, 3), "ddec": round(dde, 3), "abs": round(tot, 3),
             "rmsRA": sra, "rmsDec": sde, "ratio": round(tot / sig, 3) if sig else None,
             "mag": fnum(r.get("mag")), "logSNR": fnum(r.get("logSNR")),
             "sig_eph_ra": round(era, 4) if era is not None else None,
             "sig_eph_dec": round(ede, 4) if ede is not None else None,
             "ratio_c": round(tot / math.hypot(sig, seph), 3) if sig and seph is not None else None}
        meas.append(m)
    return meas


def eph_clean(m):
    """Measurement whose ephemeris σ is known and < SIG_EPH_CLEAN on both axes."""
    a, d = m.get("sig_eph_ra"), m.get("sig_eph_dec")
    return a is not None and d is not None and max(a, d) < SIG_EPH_CLEAN


def run_sequence(seq, piazzi, workdir):
    t0 = time.time()
    res = {"id": seq["id"], "target": seq["target"], "detected": False, "status": "ok", "tracks": [],
           "false_positives": 0}
    imgs = []
    for p in fetch_cutouts(seq, tries=12):
        hd = fits_header(p)
        nx, ny = int(hd.get("NAXIS1", 0)), int(hd.get("NAXIS2", 0))
        if min(nx, ny) < 300:
            raise RuntimeError("small cutout (%dx%d, quadrant edge)" % (nx, ny))
        if ("OBSJD" not in hd and "DATE-OBS" not in hd) or "EXPTIME" not in hd:
            raise RuntimeError("FITS without OBSJD/DATE-OBS/EXPTIME")
        imgs.append(p)
    sets = []
    for k, v in seq_settings(seq).items():
        sets += ["--set", "%s=%s" % (k, v)]
    center = ["--center-deg", "%.6f" % seq["ra"], "%.6f" % seq["dec"]]
    catp, res["catalog_source"] = prepare_catalog(piazzi, sets, center, imgs)
    out = os.path.join(workdir, seq["id"])
    os.makedirs(out, exist_ok=True)
    for f in os.listdir(out):
        os.remove(os.path.join(out, f))
    t1 = time.time()
    pr = subprocess.run([piazzi] + sets + ["--catalog", catp] + center + ["--out", out] + imgs,
                        capture_output=True, text=True, timeout=600)
    res["piazzi_seconds"] = round(time.time() - t1, 2)
    logtxt = pr.stdout + pr.stderr
    open(os.path.join(out, "log.txt"), "w").write(logtxt)
    plates, movers = parse_log(logtxt)
    res["plates"] = plates
    res["returncode"] = pr.returncode
    if pr.returncode != 0:
        res["status"] = "piazzi_rc_%d" % pr.returncode
        res["error"] = logtxt[-400:]
    if not plates or any(p["nref"] < 6 for p in plates):
        res["plate_fail"] = True
    tracks = ades_tracks(os.path.join(out, "ADESReport.psv"))
    fobj = field_objects(seq)
    known = seq["known"]
    jref = seq["jd_mid"]
    for tk, rr in tracks.items():
        best = None
        for o in known:
            ds = []
            for r in rr:
                jd = iso_to_jd(r["obsTime"])
                pra, pde = predicted(o, jd, jref)
                ds.append(sep_arcsec(float(r["ra"]), float(r["dec"]), pra, pde))
            d = statistics.median(ds)
            if best is None or d < best[0]:
                best = (d, o)
        t = {"trk": tk, "n": len(rr), "assoc": None, "assoc_sep": round(best[0], 2) if best else None}
        if best and best[0] < MATCH_ARCSEC:
            o = best[1]
            t["assoc"] = o["name"]
            try:
                t["meas"] = score_track(rr, o, seq.get("mpc", "I41"))
            except RuntimeError as e:
                t["horizons_error"] = str(e)[:300]
                t["meas"] = []
            if o["name"] == seq["target"]["name"]:
                res["detected"] = True
        else:
            res["false_positives_target_only"] = res.get("false_positives_target_only", 0) + 1
            fm = match_field(rr, fobj, jref)
            if fm:
                t["assoc_field"], t["assoc_field_sep"] = fm
            else:
                res["false_positives"] += 1
        res["tracks"].append(t)
    res["movers"] = movers
    res["nrefs"] = [p["nref"] for p in plates]
    res["seconds"] = round(time.time() - t0, 2)
    return res


def pct(v, q):
    v = sorted(v)
    if not v:
        return None
    k = (len(v) - 1) * q
    f = int(k)
    c = min(f + 1, len(v) - 1)
    return v[f] + (v[c] - v[f]) * (k - f)


def metrics(doc):
    seqs = doc["results"]
    ok = [s for s in seqs if s["status"] in ("ok",) or s["status"].startswith("piazzi")]
    det = [s for s in ok if s["detected"]]
    allm = [m for s in ok for t in s["tracks"] for m in t.get("meas", []) if t["assoc"] == s["target"]["name"]]
    km = [m for s in ok for t in s["tracks"] for m in t.get("meas", [])]
    absv = [m["abs"] for m in allm]
    rat = [m["ratio"] for m in allm if m["ratio"] is not None]
    bins = {}
    for s in ok:
        b = "%d–%d" % (int(s["target"]["V"]), int(s["target"]["V"]) + 1)
        n, d = bins.get(b, (0, 0))
        bins[b] = (n + 1, d + (1 if s["detected"] else 0))
    rms = lambda v: math.sqrt(sum(x * x for x in v) / len(v)) if v else None
    mean = lambda v: statistics.mean(v) if v else None
    clean = [m for m in allm if eph_clean(m)]
    ratc = [m["ratio_c"] for m in allm if m.get("ratio_c") is not None]
    sig_eph = [math.hypot(m["sig_eph_ra"], m["sig_eph_dec"]) for m in allm if m.get("sig_eph_ra") is not None]
    # per-axis ratio (|O−C_axis|/σ_axis; 0.674 for calibrated Gaussian errors), clean measurements only
    zra = [abs(m["dra"]) / m["rmsRA"] for m in clean if m["rmsRA"]]
    zde = [abs(m["ddec"]) / m["rmsDec"] for m in clean if m["rmsDec"]]
    within, wra, wde = [], [], []
    for s in ok:
        for t in s["tracks"]:
            ms = [m for m in t.get("meas", []) if t["assoc"] == s["target"]["name"] and m["rmsRA"] and m["rmsDec"]]
            if len(ms) < 2:
                continue
            k = math.sqrt(len(ms) / (len(ms) - 1))
            ma, md = statistics.mean(m["dra"] for m in ms), statistics.mean(m["ddec"] for m in ms)
            for m in ms:
                a, d = k * (m["dra"] - ma), k * (m["ddec"] - md)
                within.append(math.hypot(a, d) / math.hypot(m["rmsRA"], m["rmsDec"]))
                wra.append(abs(a) / m["rmsRA"])
                wde.append(abs(d) / m["rmsDec"])
    return {
        "n_sequences": len(seqs), "n_run": len(ok), "n_failed": len(seqs) - len(ok),
        "n_plate_fail": sum(1 for s in ok if s.get("plate_fail")),
        "n_detected": len(det), "detection_rate": round(len(det) / len(ok), 3) if ok else None,
        "detection_by_V": {k: {"n": v[0], "det": v[1], "rate": round(v[1] / v[0], 3)} for k, v in sorted(bins.items())},
        "n_meas_target": len(allm), "n_meas_known": len(km),
        "oc_median": pct(absv, 0.5), "oc_rms": rms(absv), "oc_p95": pct(absv, 0.95),
        "oc_rms_ra": rms([m["dra"] for m in allm]), "oc_rms_dec": rms([m["ddec"] for m in allm]),
        "oc_mean_ra": mean([m["dra"] for m in allm]),
        "oc_mean_dec": mean([m["ddec"] for m in allm]),
        "oc_median_all_known": pct([m["abs"] for m in km], 0.5),
        "ratio_median": pct(rat, 0.5), "sigma_median": pct([math.hypot(m["rmsRA"], m["rmsDec"]) for m in allm
                                                           if m["rmsRA"] and m["rmsDec"]], 0.5),
        # ephemeris error (Horizons, 1σ, axes in quadrature) and metrics that separate it from Piazzi's error:
        # (b) only measurements with σ_eph < SIG_EPH_CLEAN per axis; (c) corrected ratio |O−C|/√(σ² + σ_eph²)
        "sig_eph_median": pct(sig_eph, 0.5), "n_meas_eph": len(sig_eph),
        "n_meas_clean": len(clean), "n_seq_clean": len({id(s) for s in ok for t in s["tracks"]
                                                        for m in t.get("meas", []) if t["assoc"] == s["target"]["name"] and eph_clean(m)}),
        "oc_median_clean": pct([m["abs"] for m in clean], 0.5), "oc_rms_clean": rms([m["abs"] for m in clean]),
        "oc_p95_clean": pct([m["abs"] for m in clean], 0.95),
        "oc_mean_ra_clean": mean([m["dra"] for m in clean]), "oc_mean_dec_clean": mean([m["ddec"] for m in clean]),
        "oc_rms_ra_clean": rms([m["dra"] for m in clean]), "oc_rms_dec_clean": rms([m["ddec"] for m in clean]),
        "ratio_median_clean": pct([m["ratio"] for m in clean if m["ratio"] is not None], 0.5),
        "ratio_c_median": pct(ratc, 0.5),
        "ratio_c_median_clean": pct([m["ratio_c"] for m in clean if m.get("ratio_c") is not None], 0.5),
        "z_ra_median_clean": pct(zra, 0.5), "z_dec_median_clean": pct(zde, 0.5),
        # internal consistency: residual minus the sequence mean (the ephemeris error is constant over an arc of minutes;
        # the factor √(n/(n−1)) restores the variance removed by the mean), in units of σ — measures only Piazzi's error
        "ratio_within_median": pct(within, 0.5), "z_ra_within_median": pct(wra, 0.5), "z_dec_within_median": pct(wde, 0.5),
        "fp_per_seq": round(sum(s["false_positives"] for s in ok) / len(ok), 3) if ok else None,
        "fp_total": sum(s["false_positives"] for s in ok),
        "fp_per_seq_target_only": round(sum(s.get("false_positives_target_only", s["false_positives"]) for s in ok)
                                        / len(ok), 3) if ok else None,
        "n_known_field": sum(1 for s in ok for t in s["tracks"] if t.get("assoc_field")),
        "seconds_total": round(sum(s.get("seconds", 0) for s in seqs), 1),
        "piazzi_seconds_total": round(sum(s.get("piazzi_seconds", 0) for s in seqs), 1),
    }


def fmt(x, nd=3):
    return "—" if x is None else ("%.*f" % (nd, x) if isinstance(x, float) else str(x))


def summary_md(doc):
    m = doc["metrics"]
    L = ["# Benchmark Piazzi — %s" % doc["label"], "",
         "piazzi: `%s` (%s) · %s · fields.json: %d sequences" % (doc["piazzi"], doc.get("git", "?"), doc["date"],
                                                                 m["n_sequences"]), "",
         "| metric | value |", "|---|---|"]
    rows = [("sequences (run / data failures)", "%d (%d / %d)" % (m["n_sequences"], m["n_run"], m["n_failed"])),
            ("plate not solved (some image < 6 refs)", m["n_plate_fail"]),
            ("target detected", "%d (%.0f%%)" % (m["n_detected"], 100 * (m["detection_rate"] or 0))),
            ("target measurements", m["n_meas_target"]),
            ("|O−C| median (\")", fmt(m["oc_median"])), ("|O−C| rms (\")", fmt(m["oc_rms"])),
            ("|O−C| p95 (\")", fmt(m["oc_p95"])),
            ("O−C mean RA·cosδ / Dec (\")", "%s / %s" % (fmt(m["oc_mean_ra"]), fmt(m["oc_mean_dec"]))),
            ("O−C rms RA·cosδ / Dec (\")", "%s / %s" % (fmt(m["oc_rms_ra"]), fmt(m["oc_rms_dec"]))),
            ("median ADES sigma (\", RA/Dec in quadrature)", fmt(m["sigma_median"])),
            ("median |O−C|/sigma ratio", fmt(m["ratio_median"])),
            ("median ephemeris sigma (Horizons, 1σ, quadrature) (\")", "%s (%s measurements)" % (fmt(m.get("sig_eph_median")), m.get("n_meas_eph", "—"))),
            ("median corrected ratio |O−C|/√(σ²+σ_eph²)", fmt(m.get("ratio_c_median"))),
            ("clean measurements (σ_eph < %g\" per axis)" % SIG_EPH_CLEAN, "%s in %s sequences" % (m.get("n_meas_clean", "—"), m.get("n_seq_clean", "—"))),
            ("  |O−C| median / rms / p95 (\")", "%s / %s / %s" % (fmt(m.get("oc_median_clean")), fmt(m.get("oc_rms_clean")), fmt(m.get("oc_p95_clean")))),
            ("  O−C mean RA·cosδ / Dec (\")", "%s / %s" % (fmt(m.get("oc_mean_ra_clean")), fmt(m.get("oc_mean_dec_clean")))),
            ("  O−C rms RA·cosδ / Dec (\")", "%s / %s" % (fmt(m.get("oc_rms_ra_clean")), fmt(m.get("oc_rms_dec_clean")))),
            ("  median |O−C|/σ ratio (0.83 calibrated) / corrected", "%s / %s" % (fmt(m.get("ratio_median_clean")), fmt(m.get("ratio_c_median_clean")))),
            ("  median per-axis |O−C|/σ RA / Dec (0.674 calibrated)", "%s / %s" % (fmt(m.get("z_ra_median_clean")), fmt(m.get("z_dec_median_clean")))),
            ("internal consistency (residual − sequence mean) median |O−C|/σ / per axis RA / Dec", "%s / %s / %s" % (fmt(m.get("ratio_within_median")), fmt(m.get("z_ra_within_median")), fmt(m.get("z_dec_within_median")))),
            ("tracks with no known object (false positives) per sequence", fmt(m["fp_per_seq"])),
            ("  same, associating only target/known (old metric)", fmt(m.get("fp_per_seq_target_only"))),
            ("tracks of other SkyBoT objects in the field", m.get("n_known_field", "—")),
            ("total time (s) / piazzi only (s)", "%s / %s" % (m["seconds_total"], m["piazzi_seconds_total"]))]
    L += ["| %s | %s |" % r for r in rows]
    L += ["", "Detection by V bin:", "", "| V | n | detected | rate |", "|---|---|---|---|"]
    for k, v in m["detection_by_V"].items():
        L.append("| %s | %d | %d | %.0f%% |" % (k, v["n"], v["det"], 100 * v["rate"]))
    L += ["", "| sequence | target | V | \"/min | refs | det | |O−C| med (\") | σ_eph (\") | FP | status |", "|---|---|---|---|---|---|---|---|---|---|"]
    for s in doc["results"]:
        oc = [mm["abs"] for t in s.get("tracks", []) if t["assoc"] == s["target"]["name"] for mm in t.get("meas", [])]
        se = [math.hypot(mm["sig_eph_ra"], mm["sig_eph_dec"]) for t in s.get("tracks", []) if t["assoc"] == s["target"]["name"]
              for mm in t.get("meas", []) if mm.get("sig_eph_ra") is not None]
        L.append("| %s | %s | %.1f | %.2f | %s | %s | %s | %s | %d | %s |" % (
            s["id"], s["target"]["name"], s["target"]["V"], s["target"]["rate_arcsec_min"],
            ",".join(map(str, s.get("nrefs", []))) or "—", "yes" if s["detected"] else "no",
            fmt(statistics.median(oc) if oc else None, 2), fmt(statistics.median(se) if se else None, 3), s.get("false_positives", 0),
            s["status"] + (" (%s)" % s["error"][:80].replace("|", "/").replace("\n", " ") if s.get("error") and
                           not s["status"].startswith("piazzi") else "")))
    return "\n".join(L) + "\n"


def git_rev():
    try:
        return subprocess.run(["git", "-C", ROOT, "rev-parse", "--short", "HEAD"], capture_output=True,
                              text=True).stdout.strip()
    except OSError:
        return "?"


def cmd_run(a):
    doc_f = json.load(open(a.fields))
    seqs = doc_f["sequences"]
    if a.only:
        ids = set(a.only.split(","))
        seqs = [s for s in seqs if s["id"] in ids]
    if a.max:
        seqs = seqs[:a.max]
    piazzi = os.path.abspath(a.piazzi)
    for kv in a.set or []:
        k, _, v = kv.partition("=")
        EXTRA_SETTINGS[k] = v
    work = os.path.join(CACHE, "out", a.label)
    results = []
    for i, s in enumerate(seqs):
        log("[%d/%d] %s (%s V=%.1f)" % (i + 1, len(seqs), s["id"], s["target"]["name"], s["target"]["V"]))
        t0 = time.time()
        try:
            r = run_sequence(s, piazzi, work)
        except Exception as e:  # noqa: BLE001 — log and continue
            r = {"id": s["id"], "target": s["target"], "detected": False, "status": "data_error",
                 "error": str(e)[:400], "tracks": [], "false_positives": 0, "seconds": round(time.time() - t0, 2)}
        log("   -> %s det=%s refs=%s FP=%d" % (r["status"], r["detected"], r.get("nrefs"), r["false_positives"]))
        results.append(r)
    doc = {"label": a.label, "piazzi": rel_path(piazzi), "git": git_rev(), "date": time.strftime("%Y-%m-%d %H:%M"),
           "fields": rel_path(a.fields),
           "settings": dict(doc_f.get("settings") or ZTF_SETTINGS, **EXTRA_SETTINGS), "results": results}
    doc["metrics"] = metrics(doc)
    os.makedirs(RESULTS, exist_ok=True)
    jp = os.path.join(RESULTS, a.label + ".json")
    json.dump(doc, open(jp, "w"), indent=1, ensure_ascii=False)
    md = summary_md(doc)
    open(os.path.join(RESULTS, a.label + ".md"), "w").write(md)
    print(md)
    log("saved %s" % jp)


def cmd_rescore(a):
    """Reapplies association to all SkyBoT objects in the field to an old result (reads the ADESReport.psv files from
    cache/out/<label>), recomputes the O−C against Horizons (now with the ephemeris uncertainty) and recomputes the metrics;
    writes updated results/<label>.json/.md."""
    jp = os.path.join(RESULTS, a.label + ".json")
    doc = json.load(open(jp))
    seqs = {s["id"]: s for s in json.load(open(os.path.join(ROOT, doc["fields"])))["sequences"]}
    for r in doc["results"]:
        seq = seqs.get(r["id"])
        path = os.path.join(CACHE, "out", a.label, r["id"], "ADESReport.psv")
        if not seq or not os.path.exists(path) or r["status"] == "data_error":
            continue
        tracks = list(ades_tracks(path).values())
        if len(tracks) != len(r["tracks"]):
            log("%s: %d tracks in ADES, %d in result — skipped" % (r["id"], len(tracks), len(r["tracks"])))
            continue
        fobj = field_objects(seq)
        r["false_positives_target_only"] = sum(1 for t in r["tracks"] if not t["assoc"])
        r["false_positives"] = 0
        known = {o["name"]: o for o in seq["known"]}
        for t, rr in zip(r["tracks"], tracks):
            if t["assoc"]:
                # recomputed O−C (Horizons ephemeris position and uncertainty, cached)
                if t["assoc"] in known:
                    try:
                        t["meas"] = score_track(rr, known[t["assoc"]], seq.get("mpc", "I41"))
                        t.pop("horizons_error", None)
                    except RuntimeError as e:
                        t["horizons_error"] = str(e)[:300]
                continue
            fm = match_field(rr, fobj, seq["jd_mid"])
            if fm:
                t["assoc_field"], t["assoc_field_sep"] = fm
            else:
                r["false_positives"] += 1
    doc["metrics"] = metrics(doc)
    json.dump(doc, open(jp, "w"), indent=1, ensure_ascii=False)
    md = summary_md(doc)
    open(os.path.join(RESULTS, a.label + ".md"), "w").write(md)
    print(md)


def cmd_selftest(a):
    """Tests for association to field objects (offline)."""
    jref = 2452910.78
    o = {"name": "X", "ra": 10.0, "dec": 0.0, "dra_h": 36.0, "ddec_h": 0.0, "err": 2.0}
    far = {"name": "Y", "ra": 10.1, "dec": 0.0, "dra_h": 0.0, "ddec_h": 0.0, "err": 50.0}

    def row(jd, dra_arcsec):
        ra, de = predicted(o, jd, jref)
        yy, mo, dd, hh = jd_to_iso(jd)
        sec = hh * 3600
        t = "%04d-%02d-%02dT%02d:%02d:%06.3f" % (yy, mo, dd, int(sec // 3600), int(sec % 3600 // 60), sec % 60)
        return {"obsTime": t, "ra": ra + dra_arcsec / 3600, "dec": de}
    rr = [row(jref + k / 1440, 3.0) for k in range(4)]
    ok = True
    for what, cond in [("radius: floor MATCH_ARCSEC", field_match_radius(o) == MATCH_ARCSEC),
                       ("radius: 3σ of ephemeris, cap FIELD_MATCH_MAX", field_match_radius(far) == FIELD_MATCH_MAX),
                       ("track 3″ from prediction is associated", (match_field(rr, [far, o], jref) or ("",))[0] == "X"),
                       ("track 30″ away is not associated", match_field([row(jref + k / 1440, 30.0) for k in range(4)], [o], jref) is None),
                       ("ADES with ms", abs(iso_to_jd("2003-09-28T02:57:46.667Z") - 2452910.623456) < 1e-5)]:
        print("%s %s" % ("ok  " if cond else "FAIL", what))
        ok = ok and cond
    sys.exit(0 if ok else 1)


def cmd_summary(a):
    print(summary_md(json.load(open(a.file))))


def cmd_compare(a):
    A, B = json.load(open(a.a)), json.load(open(a.b))
    ma, mb = A["metrics"], B["metrics"]
    keys = ["n_run", "n_plate_fail", "n_detected", "detection_rate", "n_meas_target", "oc_median", "oc_rms", "oc_p95",
            "oc_rms_ra", "oc_rms_dec", "sigma_median", "ratio_median", "sig_eph_median", "ratio_c_median", "n_meas_clean",
            "oc_median_clean", "oc_rms_clean", "oc_p95_clean", "oc_mean_ra_clean", "oc_mean_dec_clean", "ratio_median_clean",
            "ratio_c_median_clean", "z_ra_median_clean", "z_dec_median_clean", "ratio_within_median", "z_ra_within_median",
            "z_dec_within_median", "fp_per_seq", "fp_per_seq_target_only",
            "piazzi_seconds_total"]
    print("| metric | %s | %s | Δ |\n|---|---|---|---|" % (A["label"], B["label"]))
    for k in keys:
        x, y = ma.get(k), mb.get(k)
        d = (y - x) if isinstance(x, (int, float)) and isinstance(y, (int, float)) else None
        print("| %s | %s | %s | %s |" % (k, fmt(x), fmt(y), fmt(d)))
    ra = {s["id"]: s for s in A["results"]}
    print("\nChanges per sequence:")
    for s in B["results"]:
        o = ra.get(s["id"])
        if o and (o["detected"] != s["detected"] or o["status"] != s["status"]):
            print("  %s: det %s -> %s, status %s -> %s" % (s["id"], o["detected"], s["detected"], o["status"],
                                                          s["status"]))



# ------------------------------------------------------------------ synthetic object injection (piazzi --inject)
# Scenarios on the real development images (ZTF dev, SDSS, NEAs dev; NEVER fields_test.json): the truth is exact
# and piazzi prints one line per object (inj|...) and per recovered measurement (injm|...: O−C against truth and ADES σ).
INJECT_SCENARIOS = {
    "slow": ["--inject-rate", "0.15", "0.8", "--inject-mag", "18", "21.5"],
    "trail": ["--inject-rate", "1", "30", "--inject-mag", "16", "19.5"],
    "faint": ["--inject-rate", "0.2", "1.5", "--inject-snr", "2", "4"],
    "cross": ["--inject-rate", "1", "15", "--inject-mag", "16", "19", "--inject-cross"],
    "color": ["--inject-rate", "0.3", "10", "--inject-mag", "17", "19.5", "--inject-color", "0.5"],
    "edge": ["--inject-rate", "1", "10", "--inject-mag", "16", "19", "--inject-edge"],
}
RATE_BINS = [(0, 1), (1, 3), (3, 10), (10, 30.01)]


def parse_inject(text):
    objs, meas, other = {}, [], None
    for ln in text.splitlines():
        f = ln.split("|")
        if f[0] == "inj" and len(f) >= 10:
            objs[int(f[1])] = {"mag": float(f[2]), "rate": float(f[3]), "pa": float(f[4]), "snr": float(f[5]),
                               "len_px": float(f[6]), "n_in": int(f[7]), "rec": f[8] == "1", "nmeas": int(f[9]),
                               "remark": f[10] if len(f) > 10 else "", "meas": []}
        elif f[0] == "injm" and len(f) >= 11:
            meas.append((int(f[1]), {"frame": int(f[2]), "dra": float(f[3]), "dde": float(f[4]), "era": float(f[5]),
                                     "ede": float(f[6]), "cra": float(f[7]), "cde": float(f[8]), "len": float(f[9]),
                                     "snr": float(f[10]), "remark": f[11] if len(f) > 11 else ""}))
        elif f[0] == "injsum" and len(f) >= 4:
            other = int(f[3])
    for i, m in meas:
        if i in objs:
            objs[i]["meas"].append(m)
    return list(objs.values()), other


def inject_one(job):
    seq, piazzi, scen, n, seed = job
    t0 = time.time()
    try:
        sets, center, imgs = seq_inputs(seq)
        catp, _ = prepare_catalog(piazzi, sets, center, imgs)
        pr = subprocess.run([piazzi] + sets + ["--catalog", catp] + center + ["--inject", str(n), "--inject-seed", str(seed)]
                            + INJECT_SCENARIOS[scen] + imgs, capture_output=True, text=True, timeout=1800)
        objs, other = parse_inject(pr.stdout)
        st = "ok" if pr.returncode == 0 else "rc_%d" % pr.returncode
        if "infeasible" in pr.stdout:
            st = "infeasible"
    except Exception as e:  # noqa: BLE001
        objs, other, st = [], None, "error: %s" % str(e)[:200]
    return {"id": seq["id"], "scenario": scen, "seed": seed, "status": st, "objects": objs, "other_tracks": other,
            "seconds": round(time.time() - t0, 1)}


def inject_metrics(runs, scen):
    """Completeness by rate and by magnitude/SNR; O−C (″) and per-axis |O−C|/σ (total ADES σ and centroid-only σ)."""
    objs = [o for r in runs if r["scenario"] == scen for o in r["objects"]]
    full = max([o["n_in"] for o in objs] or [0])
    elig = [o for o in objs if o["n_in"] >= 3] if scen == "edge" else objs
    m = {"n_obj": len(elig), "n_rec": sum(o["rec"] for o in elig)}
    m["completeness"] = m["n_rec"] / m["n_obj"] if m["n_obj"] else None
    m["by_rate"] = {}
    for lo, hi in RATE_BINS:
        b = [o for o in elig if lo <= o["rate"] < hi]
        if b:
            m["by_rate"]["%g-%g" % (lo, round(hi))] = [sum(o["rec"] for o in b), len(b)]
    key, edges = ("snr", [2, 2.5, 3, 3.5, 4.01]) if scen == "faint" else ("mag", [15, 17, 18, 19, 19.5, 20, 20.5, 21, 22])
    m["by_" + key] = {}
    for lo, hi in zip(edges, edges[1:]):
        b = [o for o in elig if lo <= o[key] < hi]
        if b:
            m["by_" + key]["%g-%g" % (lo, round(hi, 1))] = [sum(o["rec"] for o in b), len(b)]
    ms = [x for o in elig for x in o["meas"]]
    oc = [math.hypot(x["dra"], x["dde"]) for x in ms]
    m["n_meas"] = len(ms)
    m["oc_median"] = statistics.median(oc) if oc else None
    m["oc_p95"] = pct(oc, 0.95) if oc else None
    zr = [abs(x["dra"]) / x["era"] for x in ms if x["era"] > 0]
    zd = [abs(x["dde"]) / x["ede"] for x in ms if x["ede"] > 0]
    m["z_ra_median"] = statistics.median(zr) if zr else None
    m["z_dec_median"] = statistics.median(zd) if zd else None
    # pure centroid: the truth was drawn via the plate of the 1st reduction, so the O−C has almost no plate error
    zc = [abs(x["dra"]) / x["cra"] for x in ms if x["cra"] > 0] + [abs(x["dde"]) / x["cde"] for x in ms if x["cde"] > 0]
    m["z_centroid_median"] = statistics.median(zc) if zc else None
    # per axis against σ_c: the injection O−C has no plate error (the truth is drawn via the plate), so total |O−C|/σ
    # falls below 0.674 when σ_sys ⊕ σ_pred dominates (bright trails) without σ being wrong on the sky
    for ax, k, s in (("ra", "dra", "cra"), ("dec", "dde", "cde")):
        z = [abs(x[k]) / x[s] for x in ms if x[s] > 0]
        m["z_c_%s_median" % ax] = statistics.median(z) if z else None
        z = [abs(x[k]) / x[s] for x in ms if x[s] > 0 and x["len"] > 0]
        m["z_c_trail_%s_median" % ax] = statistics.median(z) if z else None
    tr = [x for x in ms if x["len"] > 0]
    m["n_meas_trail"] = len(tr)
    zt = [abs(x["dra"]) / x["era"] for x in tr if x["era"] > 0] + [abs(x["dde"]) / x["ede"] for x in tr if x["ede"] > 0]
    m["z_trail_median"] = statistics.median(zt) if zt else None
    m["oc_trail_median"] = statistics.median([math.hypot(x["dra"], x["dde"]) for x in tr]) if tr else None
    oth = [r["other_tracks"] for r in runs if r["scenario"] == scen and r["other_tracks"] is not None]
    m["other_tracks"] = sum(oth)
    m["n_runs"] = len(oth)
    m["full_frames"] = full
    return m


def inject_md(doc):
    L = ["# Injection — %s" % doc["label"], "",
         "piazzi `%s` (%s) · %s · %s · seed(s) %s" % (doc["piazzi"], doc["git"], doc["date"], doc["fields"],
                                                         ",".join(map(str, doc["seeds"]))), "",
         "| scenario | objects | recovered | completeness | by rate (\"/min) | by mag/SNR | |O−C| med/p95 (\") | "
         "|O−C|/σ RA/Dec (0.674) | /σ_c (RA/Dec) | trails: n, |O−C|/σ | other tracks |", "|" + "---|" * 11]
    for sc, m in doc["metrics"].items():
        br = " ".join("%s:%d/%d" % (k, *v) for k, v in m["by_rate"].items())
        bm = " ".join("%s:%d/%d" % (k, *v) for k, v in (m.get("by_mag") or m.get("by_snr") or {}).items())
        L.append("| %s | %d | %d | %s | %s | %s | %s/%s | %s/%s | %s | %d, %s | %d |" % (
            sc, m["n_obj"], m["n_rec"], fmt(m["completeness"], 3), br, bm, fmt(m["oc_median"]), fmt(m["oc_p95"]),
            fmt(m["z_ra_median"]), fmt(m["z_dec_median"]), fmt(m["z_centroid_median"]) + " (%s/%s)" % (fmt(m.get("z_c_ra_median")), fmt(m.get("z_c_dec_median"))), m["n_meas_trail"],
            fmt(m["z_trail_median"]), m["other_tracks"]))
    return "\n".join(L) + "\n"


def cmd_inject(a):
    fields = os.path.abspath(a.fields)
    if os.path.basename(fields) == "fields_test.json":
        sys.exit("fields_test.json is the blind set: injection forbidden")
    seqs = json.load(open(fields))["sequences"]
    if a.only:
        ids = set(a.only.split(","))
        seqs = [s for s in seqs if s["id"] in ids]
    if a.max:
        seqs = seqs[:a.max]
    for kv in a.set or []:
        k, _, v = kv.partition("=")
        EXTRA_SETTINGS[k] = v
    piazzi = os.path.abspath(a.piazzi)
    scens = list(INJECT_SCENARIOS) if a.scenario == "all" else a.scenario.split(",")
    seeds = [int(x) for x in str(a.seeds).split(",")]
    for s in seqs:  # cache and catalog serially (network), before the parallel stage
        try:
            sets, center, imgs = seq_inputs(s)
            prepare_catalog(piazzi, sets, center, imgs)
        except Exception as e:  # noqa: BLE001
            log("%s: %s" % (s["id"], str(e)[-120:]))
    jobs = [(s, piazzi, sc, a.n, sd) for sc in scens for sd in seeds for s in seqs]
    from concurrent.futures import ThreadPoolExecutor
    runs = []
    with ThreadPoolExecutor(max_workers=a.jobs) as ex:
        for r in ex.map(inject_one, jobs):
            log("%s %s s%d: %s, %d/%d recovered, %s other (%.0f s)" % (
                r["scenario"], r["id"], r["seed"], r["status"], sum(o["rec"] for o in r["objects"]), len(r["objects"]),
                r["other_tracks"], r["seconds"]))
            runs.append(r)
    doc = {"label": a.label, "piazzi": rel_path(piazzi), "git": git_rev(), "date": time.strftime("%Y-%m-%d %H:%M"),
           "fields": rel_path(fields), "n": a.n, "seeds": seeds, "settings": EXTRA_SETTINGS, "runs": runs,
           "metrics": {sc: inject_metrics(runs, sc) for sc in scens}}
    os.makedirs(RESULTS, exist_ok=True)
    jp = os.path.join(RESULTS, "inject-" + a.label + ".json") if not a.label.startswith("inject") else \
        os.path.join(RESULTS, a.label + ".json")
    json.dump(doc, open(jp, "w"), indent=1, ensure_ascii=False)
    md = inject_md(doc)
    open(jp[:-5] + ".md", "w").write(md)
    print(md)
    log("saved %s" % jp)

# ------------------------------------------------------------------ vs-mpc: the survey's own official astrometry at the MPC
MPC_GET_OBS = "https://data.minorplanetcenter.net/api/get-obs"
VSMPC_DT_MAX = 60.0  # s: pairing window Piazzi vs MPC observation from the same station


def mpc_obs(obj):
    """Optical observations (ADES_DF) of the object from the MPC public get-obs API, cached in cache/mpcobs/.
    The API is a GET with a JSON body {"desigs": [...], "output_format": ["ADES_DF"]}."""
    des = obj["num"] if obj.get("num") and obj["num"] not in ("-", "") else obj["name"]
    path = os.path.join(CACHE, "mpcobs", "%s.json" % re.sub(r"[^0-9A-Za-z_]+", "_", des))
    if os.path.exists(path) and os.path.getsize(path) > 0:
        return json.load(open(path))
    body = json.dumps({"desigs": [des], "output_format": ["ADES_DF"]}).encode()
    host = urllib.parse.urlparse(MPC_GET_OBS).netloc
    err = None
    for k in range(6):
        wait = _last_call.get(host, 0) + 2.0 - time.time()
        if wait > 0:
            time.sleep(wait)
        _last_call[host] = time.time()
        try:
            req = urllib.request.Request(MPC_GET_OBS, data=body, method="GET",
                                         headers={"User-Agent": UA, "Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=300) as r:
                d = json.loads(r.read().decode("utf-8", "replace"))
            rows = d[0]["ADES_DF"] if d and isinstance(d, list) else []
            os.makedirs(os.path.dirname(path), exist_ok=True)
            json.dump(rows, open(path + ".tmp", "w"))
            os.replace(path + ".tmp", path)
            return rows
        except (urllib.error.URLError, OSError, ValueError, KeyError, IndexError) as e:
            err = e
            code = getattr(e, "code", None)
            if code is not None and 400 <= code < 500 and code != 429:
                break
            time.sleep(min(60, 5 * 2 ** k))
    raise RuntimeError("MPC get-obs failed for %s: %s" % (des, err))


def cmd_vs_mpc(a):
    """Compares Piazzi's target measurements with the astrometry the station itself (I41 = ZTF, 645 = SDSS) submitted to
    the MPC for the same exposures. For each measurement, the station's MPC observation closest in time (|Δt| < 60 s);
    both are reduced against Horizons at their own times, so Δ = (P−H)(t_P) − (Z−H)(t_Z) already corrects
    for motion when the times differ."""
    jp = os.path.join(RESULTS, a.label + ".json")
    doc = json.load(open(jp))
    fields = a.fields or doc["fields"]
    seqs = {s["id"]: s for s in json.load(open(os.path.join(ROOT, fields) if not os.path.isabs(fields) else fields))["sequences"]}
    pairs, per_seq = [], []
    for r in doc["results"]:
        seq = seqs.get(r["id"])
        if not seq:
            continue
        tgt, stn = seq["target"], seq.get("mpc", "I41")
        meas = [m for t in r.get("tracks", []) if t["assoc"] == tgt["name"] for m in t.get("meas", [])]
        info = {"id": r["id"], "target": tgt["name"], "stn": stn, "n_piazzi": len(meas), "n_pairs": 0}
        per_seq.append(info)
        if not meas:
            info["note"] = "target not measured by Piazzi"
            continue
        try:
            obs = mpc_obs(tgt)
        except RuntimeError as e:
            info["note"] = str(e)[:200]
            continue
        st = [o for o in obs if o.get("stn") == stn and o.get("obstime") and o.get("ra") and o.get("dec")]
        jps = [iso_to_jd(m["obsTime"]) for m in meas]
        near = [o for o in st if min(abs(iso_to_jd(o["obstime"].rstrip("Z")) - j) for j in jps) * 86400 < 3600]
        info["n_mpc_stn"], info["n_mpc_night"] = len(st), len(near)
        cand = []
        for m, jpz in zip(meas, jps):
            best = None
            for o in st:
                dt = (iso_to_jd(o["obstime"].rstrip("Z")) - jpz) * 86400
                if abs(dt) < VSMPC_DT_MAX and (best is None or abs(dt) < abs(best[0])):
                    best = (dt, o)
            if best:
                cand.append((m, jpz, best[0], best[1]))
        if not cand:
            info["note"] = "no %s observation at the MPC within %d s (%d from the station the same night ±1 h)" % (
                stn, VSMPC_DT_MAX, len(near))
            continue
        eph = horizons(tgt, [jpz + dt / 86400 for _, jpz, dt, _ in cand], site=stn)
        # MPC timestamp diagnostic: the same MPC position against Horizons at Piazzi's time
        eph_p = horizons(tgt, [jpz for _, jpz, _, _ in cand], site=stn)
        for (m, jpz, dt, o), ep, epp in zip(cand, eph, eph_p):
            zra, zde = float(o["ra"]), float(o["dec"])
            zdra = (zra - ep[0]) * 3600 * math.cos(math.radians(zde))
            if abs(zdra) > 648000:
                zdra -= math.copysign(1296000 * math.cos(math.radians(zde)), zdra)
            zdde = (zde - ep[1]) * 3600
            p = {"id": r["id"], "target": tgt["name"], "obsTime": m["obsTime"], "mpc_obstime": o["obstime"],
                 "dt": round(dt, 3), "p_dra": m["dra"], "p_ddec": m["ddec"], "p_abs": m["abs"],
                 "z_dra": round(zdra, 3), "z_ddec": round(zdde, 3), "z_abs": round(math.hypot(zdra, zdde), 3),
                 "d_ra": round(m["dra"] - zdra, 3), "d_dec": round(m["ddec"] - zdde, 3),
                 "p_rmsRA": m.get("rmsRA"), "p_rmsDec": m.get("rmsDec"),
                 "z_rmsRA": fnum(o.get("rmsra")), "z_rmsDec": fnum(o.get("rmsdec")),
                 "z_precRA": fnum(o.get("precra")), "z_precDec": fnum(o.get("precdec")),
                 "z_astcat": o.get("astcat"), "z_ref": o.get("ref"), "z_subfmt": o.get("subfmt"),
                 "sig_eph_ra": m.get("sig_eph_ra"), "sig_eph_dec": m.get("sig_eph_dec")}
            p["d_abs"] = round(math.hypot(p["d_ra"], p["d_dec"]), 3)
            a_ra = (zra - epp[0]) * 3600 * math.cos(math.radians(zde))
            if abs(a_ra) > 648000:
                a_ra -= math.copysign(1296000 * math.cos(math.radians(zde)), a_ra)
            p["z_abs_at_tp"] = round(math.hypot(a_ra, (zde - epp[1]) * 3600), 3)
            # direct position comparison (same image, same Piazzi time: ignores the MPC timestamp)
            p["d_abs_pos"] = round(math.hypot(m["dra"] - a_ra, m["ddec"] - (zde - epp[1]) * 3600), 3)
            pairs.append(p)
        info["n_pairs"] = len(cand)
    rms = lambda v: math.sqrt(sum(x * x for x in v) / len(v)) if v else None
    mean = lambda v: statistics.mean(v) if v else None

    def stats(pp):
        if not pp:
            return {"n": 0}
        return {"n": len(pp), "n_seq": len({p["id"] for p in pp}),
                "d_median": pct([p["d_abs"] for p in pp], 0.5), "d_rms": rms([p["d_abs"] for p in pp]),
                "d_p95": pct([p["d_abs"] for p in pp], 0.95),
                "d_mean_ra": mean([p["d_ra"] for p in pp]), "d_mean_dec": mean([p["d_dec"] for p in pp]),
                "d_rms_ra": rms([p["d_ra"] for p in pp]), "d_rms_dec": rms([p["d_dec"] for p in pp]),
                "dt_median_abs": pct([abs(p["dt"]) for p in pp], 0.5),
                "p_oc_median": pct([p["p_abs"] for p in pp], 0.5), "z_oc_median": pct([p["z_abs"] for p in pp], 0.5),
                "p_oc_rms": rms([p["p_abs"] for p in pp]), "z_oc_rms": rms([p["z_abs"] for p in pp]),
                "p_mean_ra": mean([p["p_dra"] for p in pp]), "p_mean_dec": mean([p["p_ddec"] for p in pp]),
                "z_mean_ra": mean([p["z_dra"] for p in pp]), "z_mean_dec": mean([p["z_ddec"] for p in pp]),
                "d_pos_median": pct([p["d_abs_pos"] for p in pp], 0.5), "d_pos_rms": rms([p["d_abs_pos"] for p in pp]),
                "d_pos_p95": pct([p["d_abs_pos"] for p in pp], 0.95),
                "n_dt_gt1": sum(1 for p in pp if abs(p["dt"]) > 1),
                "dt_gt1_median": pct([p["dt"] for p in pp if abs(p["dt"]) > 1], 0.5),
                "z_oc_median_dt_gt1": pct([p["z_abs"] for p in pp if abs(p["dt"]) > 1], 0.5),
                "z_oc_median_dt_gt1_at_tp": pct([p["z_abs_at_tp"] for p in pp if abs(p["dt"]) > 1], 0.5),
                "piazzi_closer": sum(1 for p in pp if p["p_abs"] < p["z_abs"]),
                "z_rms_median": pct([math.hypot(p["z_rmsRA"], p["z_rmsDec"]) for p in pp
                                     if p["z_rmsRA"] and p["z_rmsDec"]], 0.5),
                "p_rms_median": pct([math.hypot(p["p_rmsRA"], p["p_rmsDec"]) for p in pp
                                     if p["p_rmsRA"] and p["p_rmsDec"]], 0.5),
                "n_z_rms": sum(1 for p in pp if p["z_rmsRA"]),
                "z_prec": sorted({"%s/%s" % (p["z_precRA"], p["z_precDec"]) for p in pp}),
                "z_astcat": sorted({str(p["z_astcat"]) for p in pp})}
    clean = [p for p in pairs if p["sig_eph_ra"] is not None and p["sig_eph_dec"] is not None
             and max(p["sig_eph_ra"], p["sig_eph_dec"]) < SIG_EPH_CLEAN]
    out = {"label": a.label, "fields": rel_path(fields) if os.path.isabs(fields) else fields, "date": time.strftime("%Y-%m-%d %H:%M"), "dt_max_s": VSMPC_DT_MAX,
           "metrics": stats(pairs), "metrics_clean": stats(clean), "per_seq": per_seq, "pairs": pairs}
    json.dump(out, open(os.path.join(RESULTS, a.label + "-vsmpc.json"), "w"), indent=1, ensure_ascii=False)
    md = vsmpc_md(out)
    open(os.path.join(RESULTS, a.label + "-vsmpc.md"), "w").write(md)
    print(md)


def vsmpc_md(out):
    def block(m, title):
        if not m.get("n"):
            return ["", "**%s**: no pairs." % title]
        return ["", "**%s** — %d pairs in %d sequences (median |Δt| %.2f s)" % (title, m["n"], m["n_seq"], m["dt_median_abs"]), "",
                "| metric | value |", "|---|---|",
                "| \\|Piazzi − MPC\\| median / rms / p95 (\") | %s / %s / %s |" % (fmt(m["d_median"]), fmt(m["d_rms"]), fmt(m["d_p95"])),
                "| \\|Piazzi − MPC\\| positions only (both at Piazzi's time) median / rms / p95 (\") | %s / %s / %s |" % (fmt(m["d_pos_median"]), fmt(m["d_pos_rms"]), fmt(m["d_pos_p95"])),
                "| mean bias Piazzi − MPC RA·cosδ / Dec (\") | %s / %s |" % (fmt(m["d_mean_ra"]), fmt(m["d_mean_dec"])),
                "| rms Piazzi − MPC RA·cosδ / Dec (\") | %s / %s |" % (fmt(m["d_rms_ra"]), fmt(m["d_rms_dec"])),
                "| \\|O−C\\| Horizons median Piazzi / MPC (\") | %s / %s |" % (fmt(m["p_oc_median"]), fmt(m["z_oc_median"])),
                "| \\|O−C\\| Horizons rms Piazzi / MPC (\") | %s / %s |" % (fmt(m["p_oc_rms"]), fmt(m["z_oc_rms"])),
                "| O−C mean RA·cosδ / Dec Piazzi (\") | %s / %s |" % (fmt(m["p_mean_ra"]), fmt(m["p_mean_dec"])),
                "| O−C mean RA·cosδ / Dec MPC (\") | %s / %s |" % (fmt(m["z_mean_ra"]), fmt(m["z_mean_dec"])),
                "| Piazzi closer to Horizons | %d of %d |" % (m["piazzi_closer"], m["n"]),
                "| pairs with \\|Δt\\| > 1 s: n, median Δt (s, MPC − Piazzi) | %d, %s |" % (m["n_dt_gt1"], fmt(m["dt_gt1_median"], 2)),
                "| &nbsp;&nbsp;\\|O−C\\| median MPC at MPC time / at Piazzi time (\") | %s / %s |" % (fmt(m["z_oc_median_dt_gt1"]), fmt(m["z_oc_median_dt_gt1_at_tp"])),
                "| median ADES σ Piazzi / MPC (\", quadrature) | %s / %s (%d with rms at the MPC) |" % (fmt(m["p_rms_median"]), fmt(m["z_rms_median"]), m["n_z_rms"]),
                "| MPC declared precision (precRA s / precDec \") | %s |" % ", ".join(m["z_prec"]),
                "| MPC catalog | %s |" % ", ".join(m["z_astcat"])]
    L = ["# Piazzi vs official astrometry at the MPC — %s" % out["label"], "",
         "Δ = (Piazzi − Horizons)(t_P) − (MPC − Horizons)(t_MPC); pair = observation from the same station with |Δt| < %g s." % out["dt_max_s"]]
    L += block(out["metrics"], "all pairs")
    L += block(out["metrics_clean"], "clean pairs (σ_eph < %g\" per axis)" % SIG_EPH_CLEAN)
    L += ["", "| sequence | target | Piazzi measurements | pairs | Δ med (\") | \\|O−C\\| med Piazzi / MPC (\") | closer | note |",
          "|---|---|---|---|---|---|---|---|"]
    for s in out["per_seq"]:
        pp = [p for p in out["pairs"] if p["id"] == s["id"]]
        if pp:
            pm, zm = statistics.median(p["p_abs"] for p in pp), statistics.median(p["z_abs"] for p in pp)
            L.append("| %s | %s | %d | %d | %.3f | %.3f / %.3f | %s | |" % (
                s["id"], s["target"], s["n_piazzi"], len(pp), statistics.median(p["d_abs"] for p in pp), pm, zm,
                "Piazzi" if pm < zm else "MPC"))
        else:
            L.append("| %s | %s | %d | 0 | — | — | — | %s |" % (s["id"], s["target"], s["n_piazzi"], s.get("note", "")))
    return "\n".join(L) + "\n"


# ------------------------------------------------------------------ check-ades: format validation (offline)
# Rules written from the ADES 2022 (PSV) specification and the MPC 80-column format; standard library only.
ADES_REQUIRED = ["mode", "stn", "obsTime", "ra", "dec", "astCat"]           # required in an optical observation
ADES_ID_COLS = ["permID", "provID", "trkSub"]                                # at least one
ADES_ASTCAT = {"Gaia3", "Gaia3E", "Gaia2", "Gaia1", "UCAC4", "UCAC5", "PPMXL", "URAT1", "ATLAS2", "PS1_DR1", "PS1_DR2",
               "2MASS", "Tycho2", "NOMAD", "USNOB1", "CMC15", "SDSS8", "SDSS9", "SDSS12"}
ADES_BANDS = {"Vj", "Rc", "Ic", "Bj", "Uj", "C", "V", "R", "B", "U", "I", "J", "H", "K", "W", "Y", "z", "y", "w", "o", "c",
              "G", "Gb", "Gr", "g", "r", "i", "u", "L", "Sg", "Sr", "Si", "Sz", "Su", "Pg", "Pr", "Pi", "Pz", "Pw", "Py",
              "Ao", "Ac", "Lu", "Lg", "Lr", "Li", "Lz", "Ly"}
ADES_MODES = {"CCD", "CMO", "VID", "PHO", "ENC", "PMT", "MIC", "MEA", "TDI", "UNK"}
MPC_BANDS = set(" BVRIJWUCLHKgrizyowuvGc")
RE_OBSTIME = re.compile(r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}(\.\d{1,6})?Z$")
RE_STN = re.compile(r"^([0-9A-Z]\d{2}|C\d{2}|\d{3})$")
RE_PROV = re.compile(r"^(1[89]|20)\d{2} [A-HJ-Y][A-HJ-Z](\d{1,3})?$|^[PT]-L$|^\d{4} P-L$|^\d{4} T-[123]$")
RE_PERM = re.compile(r"^[1-9]\d{0,6}$")


def _decimals(s):
    return len(s.split(".", 1)[1]) if "." in s else 0


def check_ades_text(text):
    """List of problems (str) in an ADESReport.psv."""
    bad = []
    lines = text.replace("\r\n", "\n").split("\n")
    lines = [l for l in lines if l.strip()]
    if not lines or lines[0].strip() != "# version=2022":
        bad.append("first line is not '# version=2022'")
    blocks, cur, hdr_i = {}, None, None
    for i, l in enumerate(lines):
        if l.startswith("# "):
            cur = l[2:].strip()
            blocks.setdefault(cur, {})
        elif l.startswith("! "):
            k, _, v = l[2:].partition(" ")
            if cur is None:
                bad.append("'!' line outside a block: %r" % l)
            else:
                blocks[cur].setdefault(k, []).append(v.strip())
        else:
            hdr_i = i
            break
    for blk, keys in [("observatory", ["mpcCode"]), ("submitter", ["name"]), ("measurers", ["name"]),
                      ("telescope", ["design", "aperture", "detector"])]:
        if blk not in blocks:
            bad.append("header without # %s" % blk)
            continue
        for k in keys:
            vals = blocks[blk].get(k) or []
            if not vals or not all(vals):
                bad.append("# %s without a filled-in ! %s" % (blk, k))
    tel = blocks.get("telescope", {})
    try:
        if float((tel.get("aperture") or ["0"])[0]) <= 0:
            bad.append("# telescope: aperture must be > 0 (m)")
    except ValueError:
        bad.append("# telescope: non-numeric aperture")
    code = (blocks.get("observatory", {}).get("mpcCode") or [""])[0]
    if code and not RE_STN.match(code):
        bad.append("invalid mpcCode: %r" % code)
    if hdr_i is None:
        bad.append("no column header line")
        return bad
    cols = [c.strip() for c in lines[hdr_i].split("|")]
    if len(set(cols)) != len(cols):
        bad.append("duplicate columns")
    for c in ADES_REQUIRED:
        if c not in cols:
            bad.append("missing required column: %s" % c)
    if not any(c in cols for c in ADES_ID_COLS):
        bad.append("no identification column (permID/provID/trkSub)")
    if "rmsRA" not in cols or "rmsDec" not in cols:
        bad.append("no rmsRA/rmsDec")
    if "mag" in cols and "band" not in cols:
        bad.append("mag without band")
    trk_ids, prev = {}, None
    for n, l in enumerate(lines[hdr_i + 1:], hdr_i + 2):
        f = [x.strip() for x in l.split("|")]
        if len(f) != len(cols):
            bad.append("line %d: %d fields, header has %d" % (n, len(f), len(cols)))
            continue
        r = dict(zip(cols, f))
        w = "linha %d: " % n
        for c in ADES_REQUIRED:
            if c in r and not r[c]:
                bad.append(w + "%s vazio" % c)
        if not any(r.get(c) for c in ADES_ID_COLS):
            bad.append(w + "no permID/provID/trkSub")
        trk = r.get("trkSub", "")
        if len(trk) > 8 or (trk and not re.match(r"^[A-Za-z0-9]+$", trk)):
            bad.append(w + "invalid trkSub %r (≤ 8 alphanumerics)" % trk)
        perm, prov = r.get("permID", ""), r.get("provID", "")
        if perm and not RE_PERM.match(perm):
            bad.append(w + "spurious permID %r" % perm)
        if prov and not RE_PROV.match(prov):
            bad.append(w + "spurious provID %r" % prov)
        key = (perm, prov)
        if trk:
            # trkSub unique per track: a new track starts when time does not advance
            new = prev is None or prev[0] != trk or r.get("obsTime", "") <= prev[1]
            if new and trk in trk_ids:
                bad.append(w + "trkSub %s repeated in another track" % trk)
            if not new and trk_ids.get(trk) != key:
                bad.append(w + "trkSub %s with inconsistent identification" % trk)
            trk_ids.setdefault(trk, key)
            prev = (trk, r.get("obsTime", ""))
        if r.get("mode") and r["mode"] not in ADES_MODES:
            bad.append(w + "invalid mode %r" % r["mode"])
        if r.get("mode") and r["mode"] != "CCD":
            bad.append(w + "mode %r (expected CCD)" % r["mode"])
        if r.get("stn") and (not RE_STN.match(r["stn"]) or (code and r["stn"] != code)):
            bad.append(w + "stn %r invalid/different from mpcCode" % r["stn"])
        if r.get("obsTime") and not RE_OBSTIME.match(r["obsTime"]):
            bad.append(w + "obsTime not ISO 8601 with Z: %r" % r["obsTime"])
        try:
            ra, de = float(r["ra"]), float(r["dec"])
            if not (0 <= ra < 360) or not (-90 <= de <= 90):
                bad.append(w + "ra/dec out of range")
            if not re.match(r"^[+-]\d", r["dec"]):
                bad.append(w + "dec without sign")
            sra, sde = float(r["rmsRA"]), float(r["rmsDec"])
            if not (sra > 0 and sde > 0):
                bad.append(w + "rmsRA/rmsDec must be > 0")
            else:
                stp_ra = 10 ** -_decimals(r["ra"]) * 3600 * math.cos(math.radians(de))
                stp_de = 10 ** -_decimals(r["dec"]) * 3600
                # rounding (step/√12) ≤ 1/3 of σ and non-spurious precision (step ≥ σ/1000)
                for ax, st, sg in (("ra", stp_ra, sra), ("dec", stp_de, sde)):
                    if st / math.sqrt(12) > sg / 3:
                        bad.append(w + "%s precision too coarse for rms %.3f″ (step %.4f″)" % (ax, sg, st))
                    if st < sg / 1000:
                        bad.append(w + "%s has spurious precision for rms %.3f″" % (ax, sg))
                for ax, s, sg in (("rmsRA", r["rmsRA"], sra), ("rmsDec", r["rmsDec"], sde)):
                    if 10 ** -_decimals(s) > 0.2 * sg:
                        bad.append(w + "%s over-rounded (%s)" % (ax, s))
        except (KeyError, ValueError):
            bad.append(w + "non-numeric ra/dec/rmsRA/rmsDec")
        ac = r.get("astCat", "")
        if ac and ac not in ADES_ASTCAT:
            bad.append(w + "astCat not allowed %r" % ac)
        if r.get("mag"):
            try:
                mg = float(r["mag"])
                if not (-5 < mg < 35):
                    bad.append(w + "mag out of range")
                if r.get("rmsMag") and float(r["rmsMag"]) <= 0:
                    bad.append(w + "rmsMag ≤ 0")
            except ValueError:
                bad.append(w + "non-numeric mag")
            if r.get("band") not in ADES_BANDS:
                bad.append(w + "invalid band %r" % r.get("band"))
            if not r.get("photCat"):
                bad.append(w + "mag without photCat")
            elif r["photCat"] not in ADES_ASTCAT:
                bad.append(w + "photCat not allowed %r" % r["photCat"])
        elif r.get("band") or r.get("rmsMag"):
            bad.append(w + "band/rmsMag without mag")
        for c, lo in (("exp", 0), ("seeing", 0), ("photAp", 0), ("nStars", 0)):
            if r.get(c):
                try:
                    if float(r[c]) <= lo:
                        bad.append(w + "%s ≤ %g" % (c, lo))
                except ValueError:
                    bad.append(w + "%s non-numeric" % c)
    return bad


def check_mpc_text(text):
    """List of problems in an MPCReport.txt (COD/OBS/... header + 80-column lines)."""
    bad = []
    code = None
    for n, l in enumerate(text.replace("\r\n", "\n").split("\n"), 1):
        if not l:
            continue
        w = "linha %d: " % n
        if re.match(r"^(COD|CON|OBS|MEA|TEL|NET|ACK|AC2|BND|COM|NUM)( |$)", l):
            if l.startswith("COD"):
                code = l[4:].strip()
            if len(l) > 80:
                bad.append(w + "header longer than 80 columns")
            if l[:3] in ("OBS", "MEA", "TEL", "NET") and not l[4:].strip():
                bad.append(w + "%s vazio" % l[:3])
            continue
        if len(l) != 80:
            bad.append(w + "%d columns (expected 80)" % len(l))
            continue
        num, prov = l[0:5], l[5:12]
        if not (num.strip() or prov.strip()):
            bad.append(w + "no designation (columns 1–12)")
        if num.strip() and not re.match(r"^([0-9A-Za-z]\d{4}|~[0-9A-Za-z]{4})$", num):
            bad.append(w + "invalid packed number %r" % num)
        if prov.strip() and not re.match(r"^[A-Za-z0-9]{1,7} *$", prov):
            bad.append(w + "invalid provisional designation %r" % prov)
        if l[13] not in " AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQqRrSsTtUuVvWwXxYyZz0123456789!\"#%&'()*+,-./:;<=>?@[]^_`{|}~$":
            bad.append(w + "invalid note 1")
        if l[14] not in "CBVPTMAXx":
            bad.append(w + "invalid note 2 (column 15) %r" % l[14])
        elif l[14] != "C":
            bad.append(w + "column 15 should be C (CCD)")
        if not re.match(r"^\d{4} \d{2} \d{2}\.\d{5,6} ?$", l[15:32]):
            bad.append(w + "date outside columns 16–32: %r" % l[15:32])
        if not re.match(r"^\d{2} \d{2} \d{2}\.\d{2,3} *$", l[32:44]):
            bad.append(w + "RA outside columns 33–44: %r" % l[32:44])
        if not re.match(r"^[+-]\d{2} \d{2} \d{2}\.\d{1,2} *$", l[44:56]):
            bad.append(w + "Dec outside columns 45–56: %r" % l[44:56])
        if l[56:65].strip():
            bad.append(w + "columns 57–65 should be blank")
        for ax, fld in (("AR", l[32:44]), ("Dec", l[44:56])):
            pp = fld.split()
            if len(pp) == 3 and re.match(r"^[+-]?\d+$", pp[1]) and re.match(r"^\d+(\.\d+)?$", pp[2]) and \
                    (int(pp[1]) >= 60 or float(pp[2]) >= 60 or (ax == "AR" and int(pp[0]) >= 24)):
                bad.append(w + "%s has a sexagesimal field out of range: %r" % (ax, fld))
        mag, band = l[65:70], l[70]
        if mag.strip():
            if not re.match(r"^\d{1,2}\.\d{1,2} *$", mag.lstrip()) or not mag[0].isdigit():
                bad.append(w + "magnitude outside columns 66–70: %r" % mag)
            if band not in MPC_BANDS or band == " ":
                bad.append(w + "invalid band (column 71) %r" % band)
        elif band.strip():
            bad.append(w + "band without magnitude")
        stn = l[77:80]
        if not RE_STN.match(stn) or (code and stn != code):
            bad.append(w + "observatory code (78–80) invalid/different from COD: %r" % stn)
    if code is None:
        bad.append("no COD line")
    return bad


def _latest_label(fields):
    rel = rel_path(fields)
    best = None
    for fn in os.listdir(RESULTS):
        if not fn.endswith(".json") or fn.startswith("inject"):
            continue
        try:
            d = json.load(open(os.path.join(RESULTS, fn)))
        except ValueError:
            continue
        if d.get("fields") == rel and os.path.isdir(os.path.join(CACHE, "out", d.get("label", ""))):
            t = os.path.getmtime(os.path.join(RESULTS, fn))
            if best is None or t > best[0]:
                best = (t, d["label"])
    return best and best[1]


def cmd_check_ades(a):
    files = []
    if a.ades or a.mpc:
        files = [(p, check_ades_text) for p in a.ades or []] + [(p, check_mpc_text) for p in a.mpc or []]
    else:
        label = a.label or _latest_label(a.fields)
        if not label:
            sys.exit("no result with --fields %s in results/ (use --label)" % a.fields)
        base = os.path.join(CACHE, "out", label)
        ids = [s["id"] for s in json.load(open(a.fields))["sequences"]]
        for i in ids:
            for fn, fx in (("ADESReport.psv", check_ades_text), ("MPCReport.txt", check_mpc_text)):
                p = os.path.join(base, i, fn)
                if os.path.exists(p):
                    files.append((p, fx))
        log("check-ades: %s (%d files)" % (label, len(files)))
    nbad = 0
    for p, fx in files:
        probs = fx(open(p, encoding="utf-8", errors="replace").read())
        nbad += bool(probs)
        for x in probs[:a.max_per_file]:
            print("%s: %s" % (rel_path(p), x))
        if len(probs) > a.max_per_file:
            print("%s: ... %d more" % (rel_path(p), len(probs) - a.max_per_file))
    print("check-ades: %d/%d files with problems" % (nbad, len(files)))
    sys.exit(1 if nbad else 0)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    p = sp.add_parser("select")
    p.add_argument("--n", type=int, default=30)
    p.add_argument("--per-bin", type=int, default=7, help="maximum sequences per integer V bin")
    p.add_argument("--probes", type=int, default=400)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--out", default=FIELDS)
    p.add_argument("--append", action="store_true", help="append to the existing --out instead of overwriting")
    p.add_argument("--exclude", action="append", help="fields file whose sequences (field+night) to avoid")
    p.add_argument("--vmin", type=float, default=15.0)
    p.add_argument("--vmax", type=float, default=20.0)
    p.add_argument("--span-min-h", type=float, default=0.5)
    p = sp.add_parser("inject", help="synthetic object injection into the development images (completeness, O−C)")
    p.add_argument("--fields", default=FIELDS)
    p.add_argument("--scenario", default="all", help="|".join(INJECT_SCENARIOS) + ", comma-separated list or all")
    p.add_argument("--n", type=int, default=6, help="objects per sequence")
    p.add_argument("--seeds", default="1")
    p.add_argument("--label", default="inject")
    p.add_argument("--piazzi", default=os.path.join(ROOT, "cpp", "build", "piazzi"))
    p.add_argument("--only")
    p.add_argument("--max", type=int)
    p.add_argument("--jobs", type=int, default=6)
    p.add_argument("--set", action="append")
    p = sp.add_parser("select-targets", help="sequences from known targets (NEAs, bright objects) via Horizons")
    p.add_argument("--objects", required=True, help="comma-separated designations, or 'cad' (JPL close approaches)")
    p.add_argument("--out", required=True)
    p.add_argument("--append", action="store_true")
    p.add_argument("--exclude", action="append")
    p.add_argument("--n", type=int, default=10)
    p.add_argument("--per-object", type=int, default=1)
    p.add_argument("--vmin", type=float, default=15.0)
    p.add_argument("--vmax", type=float, default=19.5)
    p.add_argument("--rate-min", type=float, default=0.0, help='"/min')
    p.add_argument("--rate-max", type=float, default=30.0)
    p.add_argument("--span-min-h", type=float, default=0.33)
    p.add_argument("--windows", type=int, default=3, help="date windows tried per object")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--label", default="targets")
    p.add_argument("--exclude-objects", help="designations/numbers to avoid, comma-separated")
    p.add_argument("--date-min", default="2019-03-01", help="start of the close-approach search (--objects cad)")
    p.add_argument("--date-max", default="2023-10-01")
    p = sp.add_parser("select-sdss", help="independent TEST set (SDSS DR17, MPC 645) via SkyServer + SkyBoT")
    p.add_argument("--out", default=FIELDS_SDSS)
    p.add_argument("--n", type=int, default=35)
    p.add_argument("--per-bin", type=int, default=6, help="maximum sequences per integer V bin")
    p.add_argument("--per-run", type=int, default=2, help="maximum sequences per SDSS run")
    p.add_argument("--probes", type=int, default=600)
    p.add_argument("--seed", type=int, default=7)
    p.add_argument("--vmin", type=float, default=15.0)
    p.add_argument("--vmax", type=float, default=21.0)
    p.add_argument("--beta-max", type=float, default=12.0, help="maximum |ecliptic latitude| of the field (degrees)")
    p.add_argument("--filters", default=SDSS_FILTERS, help="subset of riuzg (temporal order is enforced)")
    p.add_argument("--margin", type=float, default=60, help="margin (px) from target to frame edge in every filter")
    p.add_argument("--modulus", type=int, default=40, help="Field table sample: field %% modulus = residue")
    p.add_argument("--residue", type=int, default=7)
    p = sp.add_parser("run")
    p.add_argument("--label", required=True)
    p.add_argument("--piazzi", default=os.path.join(ROOT, "cpp", "build", "piazzi"))
    p.add_argument("--fields", default=FIELDS)
    p.add_argument("--only")
    p.add_argument("--max", type=int)
    p.add_argument("--set", action="append", help="piazzi setting K=V on top of the set's settings (experiments)")
    p = sp.add_parser("fetch")
    p.add_argument("--fields", default=FIELDS)
    p.add_argument("--passes", type=int, default=20)
    p.add_argument("--sleep", type=float, default=60, help="pause (s) between passes")
    p.add_argument("--piazzi", help="if given, also downloads the Gaia catalog for each sequence (only the step that prints the URL)")
    p = sp.add_parser("compare")
    p.add_argument("a")
    p.add_argument("b")
    p = sp.add_parser("rescore", help="reapplies association to all field SkyBoT objects in results/<label>.json")
    p.add_argument("label")
    p = sp.add_parser("vs-mpc", help="compares with the station's official astrometry (I41/645) at the MPC (get-obs API)")
    p.add_argument("--label", required=True)
    p.add_argument("--fields", help="fields file (default: the one recorded in results/<label>.json)")
    p = sp.add_parser("check-ades", help="validates the format of generated ADESReport.psv/MPCReport.txt (offline)")
    p.add_argument("--fields", default=FIELDS)
    p.add_argument("--label", help="cache/out/<label> (default: the most recent result with --fields)")
    p.add_argument("--ades", action="append", help="standalone ADES file (no --fields/--label needed)")
    p.add_argument("--mpc", action="append", help="standalone MPC 80-column file")
    p.add_argument("--max-per-file", type=int, default=5)
    sp.add_parser("selftest", help="association tests (offline)")
    p = sp.add_parser("summary")
    p.add_argument("file")
    a = ap.parse_args()
    {"select": cmd_select, "select-sdss": cmd_select_sdss, "select-targets": cmd_select_targets, "run": cmd_run, "fetch": cmd_fetch, "compare": cmd_compare, "rescore": cmd_rescore, "selftest": cmd_selftest, "summary": cmd_summary, "inject": cmd_inject, "vs-mpc": cmd_vs_mpc, "check-ades": cmd_check_ades}[a.cmd](a)


if __name__ == "__main__":
    main()
