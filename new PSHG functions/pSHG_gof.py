# -*- coding: utf-8 -*-
"""
pshg_gof.py -- per-pixel goodness-of-fit and trust metrics for
polarisation-resolved SHG (single input-rotation model).

Rebuild of the module validated 2026-09-10 on 12-frame RTT data.
Handles arbitrary N and angle set; dof = N - 5. Self-calibrates the
shot-noise gain so no detector characterisation is needed.
"""
import numpy as np
from scipy import stats


def pshg_gof(pshg, angles_rad, gain=None, read_noise2=0.0,
             sat_level=32600, signal_percentile=60, weighted_fit=False):
    pshg = np.asarray(pshg, float)
    N, H, W = pshg.shape
    t = np.asarray(angles_rad, float)
    Y = pshg.reshape(N, -1)

    X = np.column_stack([np.ones(N), np.cos(2*t), np.sin(2*t),
                                     np.cos(4*t), np.sin(4*t)])
    beta, *_ = np.linalg.lstsq(X, Y, rcond=None)
    fit = X @ beta

    if weighted_fit:
        var = np.clip(gain_var(fit, 1.0, read_noise2), 1e-9, None)
        beta, fit = _wls(X, Y, 1.0/var)

    resid = Y - fit
    RSS = (resid**2).sum(0)

    a0 = beta[0]
    a2, b2, a4, b4 = beta[1], beta[2], beta[3], beta[4]
    R2 = np.hypot(a2, b2)
    R4 = np.hypot(a4, b4)

    saturated = (pshg >= sat_level).any(0).ravel()
    sig_mask = (a0 > np.percentile(a0, signal_percentile)) & (~saturated)

    dof = N - 5
    if gain is None:
        w0 = 1.0 / np.clip(gain_var(fit, 1.0, read_noise2), 1e-9, None)
        chi2_raw = (w0 * resid**2).sum(0)
        gain = np.median(chi2_raw[sig_mask]) / stats.chi2.median(dof)

    var = np.clip(gain * np.clip(fit, 0, None) + read_noise2, 1e-9, None)
    chi2 = (resid**2 / var).sum(0)

    calibration = stats.chi2.median(dof) / np.median(chi2[sig_mask])
    chi2 = chi2 * calibration
    chi2_red = chi2 / dof
    pvalue = stats.chi2.sf(chi2, dof)

    g_eff = gain / calibration
    s_amp = np.sqrt(2.0 * g_eff * np.clip(a0, 1.0, None) / N)
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
        gain_used=float(gain), calibration_factor=float(calibration), dof=dof,
    )


def gain_var(fit, g, read2):
    return g * np.clip(fit, 0, None) + read2

def _wls(X, Y, W):
    N, P = X.shape
    Npix = Y.shape[1]
    beta = np.empty((P, Npix)); CH = 200000
    for s in range(0, Npix, CH):
        e = min(s+CH, Npix); w = W[:, s:e]; y = Y[:, s:e]
        XtWX = np.moveaxis(np.einsum('ki,kp,kj->ijp', X, w, X), 2, 0)
        XtWy = np.einsum('ki,kp,kp->ip', X, w, y).T
        beta[:, s:e] = np.linalg.solve(XtWX, XtWy[..., None])[..., 0].T
    return beta, X @ beta
