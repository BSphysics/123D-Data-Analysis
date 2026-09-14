# -*- coding: utf-8 -*-
"""
pshg_gof.py -- per-pixel goodness-of-fit and trust metrics for
polarisation-resolved SHG (single input-rotation model).

v2: robust to LOW-SIGNAL / baseline-subtracted data. The noise model is now
affine,  var = gain*max(signal,0) + read_noise2,  with BOTH terms estimated
from the data by a residual photon-transfer fit. This replaces the old pure
shot-noise self-calibration, which detonated on signed/near-zero data (pixels
with fit<=0 hit the inverse-variance floor and blew the gain up to ~1e11).
On bright data read_noise2 -> ~0 and the behaviour matches v1.
"""
import numpy as np
from scipy import stats


def pshg_gof(pshg, angles_rad, gain=None, read_noise2=None,
             sat_level=32600, signal_percentile=60):
    pshg = np.asarray(pshg, float)
    N, H, W = pshg.shape
    t = np.asarray(angles_rad, float)
    Y = pshg.reshape(N, -1)

    X = np.column_stack([np.ones(N), np.cos(2*t), np.sin(2*t),
                                     np.cos(4*t), np.sin(4*t)])
    beta, *_ = np.linalg.lstsq(X, Y, rcond=None)
    fit = X @ beta
    resid = Y - fit
    RSS = (resid**2).sum(0)

    a0 = beta[0]
    a2, b2, a4, b4 = beta[1], beta[2], beta[3], beta[4]
    R2 = np.hypot(a2, b2)
    R4 = np.hypot(a4, b4)

    saturated = (pshg >= sat_level).any(0).ravel()
    sig_mask = (a0 > np.percentile(a0, signal_percentile)) & (~saturated)

    dof = N - 5
    # ---- affine noise model: var = gain*max(signal,0) + read_noise2 -----
    if gain is None or read_noise2 is None:
        gain, read_noise2 = _estimate_affine_noise(a0, RSS, dof)

    var = np.clip(gain * np.clip(fit, 0, None) + read_noise2, 1e-9, None)
    chi2 = (resid**2 / var).sum(0)

    # gentle recentre so p-values are calibrated (should be ~1 if fit is good)
    calibration = stats.chi2.median(dof) / np.median(chi2[sig_mask])
    chi2 = chi2 * calibration
    chi2_red = chi2 / dof
    pvalue = stats.chi2.sf(chi2, dof)

    # per-component amplitude noise, INCLUDING the read/offset floor -------
    g_eff = gain / calibration
    c_eff = read_noise2 / calibration
    var_amp = g_eff * np.clip(a0, 0, None) + c_eff          # per-frame var at DC
    s_amp = np.sqrt(2.0 * np.clip(var_amp, 1e-12, None) / N)
    sigma_phi2 = np.rad2deg(0.5 * s_amp / np.clip(R2, 1e-9, None))

    Phi2 = np.rad2deg(0.5 * np.arctan2(b2, a2)) % 180.0
    Phi4 = np.rad2deg(0.25 * np.arctan2(b4, a4)) % 45.0
    dphi = (Phi4 - Phi2 + 22.5) % 45.0 - 22.5
    phase_mismatch = np.abs(dphi)

    trust = sig_mask & (pvalue > 1e-3) & (sigma_phi2 < 10.0)

    r = lambda x: x.reshape(H, W)
    return dict(
        a0=r(a0), R2=r(R2), R4=r(R4), I2=r(R2/np.clip(a0, 1e-9, None)),
        Phi2=r(Phi2), Phi4=r(Phi4),
        chi2=r(chi2), chi2_red=r(chi2_red), pvalue=r(pvalue),
        sigma_phi2=r(sigma_phi2), phase_mismatch=r(phase_mismatch),
        saturated=r(saturated), trust=r(trust), sig_mask=r(sig_mask),
        gain_used=float(gain), read_noise2_used=float(read_noise2),
        calibration_factor=float(calibration), dof=dof,
    )


def _estimate_affine_noise(a0, RSS, dof, min_per_bin=200, nbins=30):
    """Fit var = g*max(signal,0) + c from the residual photon-transfer curve.
    Uses the per-bin MEDIAN of RSS/dof (robust to model-violating pixels) and
    de-biases the median-of-chi2 factor. Falls back to a pure-shot estimate if
    there are too few populated bins."""
    rssd = RSS / dof
    factor = stats.chi2.median(dof) / dof        # median(rssd) = factor * var
    lo, hi = np.percentile(a0, 5), np.percentile(a0, 99)
    edges = np.linspace(lo, hi, nbins)
    ctr, vb = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        m = (a0 >= a) & (a0 < b)
        if m.sum() > min_per_bin:
            ctr.append(0.5 * (a + b))
            vb.append(np.median(rssd[m]) / factor)
    ctr, vb = np.asarray(ctr), np.asarray(vb)
    if ctr.size < 3:                              # fallback: pure shot noise
        g = np.median(rssd) / factor / max(np.median(np.clip(a0, 1, None)), 1.0)
        return float(max(g, 1e-6)), 0.0
    A = np.column_stack([np.clip(ctr, 0, None), np.ones_like(ctr)])
    (g, c), *_ = np.linalg.lstsq(A, vb, rcond=None)
    return float(max(g, 1e-6)), float(max(c, 0.0))
