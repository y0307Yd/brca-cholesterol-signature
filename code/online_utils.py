# -*- coding: utf-8 -*-
"""Shared helpers for the online finishing pipeline.

Run with the analysis python:
    C:\\Users\\Y\\.codex\\py311\\python.exe work\\finish_*.py
"""
import gzip
import json
import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

UA = {
    "User-Agent": ("Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
                   "research-bioinformatics/1.0 (+https://www.ncbi.nlm.nih.gov/)")
}

OUT = Path(__file__).resolve().parent.parent / "outputs" / "chol_metab_signature"
FIG = OUT / "figures"
WORK = Path(__file__).resolve().parent
DATA = Path(__file__).resolve().parent.parent / "data"


def check_net(urls=None, timeout=10):
    if urls is None:
        urls = [
            "https://ftp.ncbi.nlm.nih.gov/",
            "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/",
            "https://www.cbioportal.org/",
            "https://example.com/",
        ]
    res = []
    for u in urls:
        ok, status = False, None
        r = subprocess.run(
            ["curl.exe", "--ssl-no-revoke", "-sS", "-o", "NUL", "-w", "%{http_code}",
             "--max-time", str(timeout), u],
            capture_output=True, text=True)
        if r.returncode == 0 and r.stdout.strip() not in ("", "000"):
            ok, status = True, r.stdout.strip()
        else:
            status = f"curl rc={r.returncode} {r.stderr.strip()[:80]}"
        res.append({"url": u, "ok": ok, "status": status})
    return res


def download(url, dest, max_try=3, timeout=14400, chunk=1 << 20):
    """Resumable download via curl (--ssl-no-revoke). Returns dest path."""
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    last = None
    for attempt in range(1, max_try + 1):
        cmd = ["curl.exe", "--ssl-no-revoke", "-f", "-L", "--retry", "5",
               "--retry-delay", "3", "--retry-all-errors", "-C", "-",
               "--max-time", str(timeout), "-o", str(tmp), url]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode == 0:
            tmp.replace(dest)
            return str(dest)
        last = f"curl rc={r.returncode}: {r.stderr.strip()[-300:]}"
        print(f"[download] attempt {attempt}/{max_try}: {url}\n    {last}")
        if attempt == max_try:
            break
        time.sleep(3 * attempt)
    raise RuntimeError(f"download failed: {url} ({last})")


def read_text_maybe_gzip(path):
    p = Path(path)
    with open(p, "rb") as f:
        magic = f.read(2)
    if magic == b"\x1f\x8b":
        with gzip.open(p, "rt", encoding="utf-8", errors="replace") as f:
            return f.read()
    return p.read_text(encoding="utf-8", errors="replace")


def int_transform(x):
    """Rank-based inverse-normal transform (per-gene, cross-platform comparable)."""
    x = np.asarray(x, dtype=float)
    n = len(x)
    ranks = stats.rankdata(x, method="average")
    return stats.norm.ppf((ranks - 0.5) / n)


def coxph(time, event, Xcov, labels, max_iter=120):
    """Newton-Raphson Cox PH (Breslow) -> dict of HR/CI/p per covariate."""
    time = np.asarray(time, dtype=float)
    event = np.asarray(event, dtype=float)
    Xcov = np.asarray(Xcov, dtype=float)
    n, k = Xcov.shape
    idx = np.argsort(time, kind="mergesort")
    t, d, X = time[idx], event[idx], Xcov[idx]
    beta = np.zeros(k)
    for it in range(max_iter):
        et = np.exp(X @ beta)
        rev = np.cumsum(et[::-1])[::-1]
        grad = np.zeros(k)
        H = np.zeros((k, k))
        for u in np.unique(t):
            fi = int(np.searchsorted(t, u, side="left"))
            ev = (t == u) & (d == 1)
            ne = int(ev.sum())
            if ne == 0:
                continue
            xev = X[ev]
            sx = xev.sum(axis=0)
            risk = rev[fi]
            if risk <= 0:
                continue
            xr = X[fi:]
            w = et[fi:] / risk
            xw = (xr * w[:, None]).sum(axis=0)
            grad += sx - ne * xw
            xxw = np.einsum("ij,ik,i->jk", xr, xr, w)
            H += ne * (xxw - np.outer(xw, xw))
        try:
            step = np.linalg.solve(H + 1e-9 * np.eye(k), grad)
        except np.linalg.LinAlgError:
            step = np.linalg.lstsq(H + 1e-9 * np.eye(k), grad, rcond=None)[0]
        beta = beta + step
        if np.max(np.abs(step)) < 1e-8:
            break
    try:
        I = np.linalg.inv(H + 1e-9 * np.eye(k))
        se = np.sqrt(np.maximum(np.diag(I), 1e-12))
    except np.linalg.LinAlgError:
        se = np.full(k, np.nan)
    out = {}
    for j, lab in enumerate(labels):
        hr = float(np.exp(beta[j]))
        lo = float(np.exp(beta[j] - 1.96 * se[j]))
        hi = float(np.exp(beta[j] + 1.96 * se[j]))
        z = beta[j] / se[j] if se[j] > 0 else np.nan
        p = float(2 * stats.norm.sf(abs(z))) if not np.isnan(z) else np.nan
        out[lab] = {"HR": hr, "CI95": [lo, hi], "beta": float(beta[j]),
                    "se": float(se[j]), "p": p}
    return out


def logrank(time, event, group):
    t = np.asarray(time, dtype=float)
    d = np.asarray(event, dtype=float)
    g = np.asarray(group, dtype=float)
    o = np.argsort(t, kind="mergesort")
    t, d, g = t[o], d[o], g[o]
    O1 = E = V = 0.0
    for u in np.unique(t):
        atrisk = t >= u
        n = int(atrisk.sum())
        n1 = int(g[atrisk].sum())
        ev = d[t == u]
        ev1 = ev[g[t == u] == 1].sum()
        e = ev.sum()
        if n > 1 and e > 0 and 0 < n1 < n:
            O1 += ev1
            E += e * n1 / n
            V += (n1 * (n - n1) * e * (n - e)) / (n * n * (n - 1))
    z = (O1 - E) / np.sqrt(V) if V > 0 else np.nan
    chi = z * z if not np.isnan(z) else np.nan
    p = float(stats.chi2.sf(chi, 1)) if not np.isnan(chi) else np.nan
    return {"chi2": float(chi) if not np.isnan(chi) else np.nan, "p": p}


def auc_ci(y, p, n_boot=2000, seed=20260809):
    from sklearn.metrics import roc_auc_score
    y = np.asarray(y, dtype=int)
    p = np.asarray(p, dtype=float)
    if len(np.unique(y)) < 2:
        return np.nan, (np.nan, np.nan)
    auc = float(roc_auc_score(y, p))
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(n_boot):
        ii = rng.integers(0, len(y), len(y))
        if len(np.unique(y[ii])) < 2:
            continue
        boots.append(roc_auc_score(y[ii], p[ii]))
    if boots:
        lo, hi = np.percentile(boots, [2.5, 97.5])
    else:
        lo = hi = np.nan
    return auc, (float(lo), float(hi))


def save_patch(patch, key):
    """Merge patch dict into work/v13_patch.json."""
    path = WORK / "v13_patch.json"
    data = json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}
    data[key] = patch
    path.write_text(json.dumps(data, ensure_ascii=False, indent=1),
                    encoding="utf-8")
    print("[patch]", key, "->", path)


def load_patch():
    path = WORK / "v13_patch.json"
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def fmt_p(p, decimals=3):
    if p is None or (isinstance(p, float) and np.isnan(p)):
        return "NA"
    p = float(p)
    if p < 0.001:
        return f"{p:.2e}"
    return f"{p:.{decimals}f}"
