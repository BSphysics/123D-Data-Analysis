# -*- coding: utf-8 -*-
"""
snr_mask.py -- robust, tunable, smoothed trust mask for pSHG I2/phi2.

v2: reads the affine noise model (gain + read_noise2) from pshg_gof, so the
per-component noise includes the read/baseline floor. This is what lets it work
on low-signal (read-noise-limited) data instead of masking the whole frame.

Reliability metric = propagated orientation error sigma(phi2), which also
bounds the I2 error (sigma_I2/I2 = 2*sigma_phi2[rad]). It is in physical units
(degrees) so the same threshold transfers between acquisitions.

Tune first:  sigma_phi2_max (strictness), pool_sigma (pooling/smoothness)
Set&forget:  min_feature (smallest region/hole), close_radius

sigma_phi2_max is an error-bar threshold, in real degrees. At each pixel the fit gives you a fibre orientation Φ₂, and — from that pixel's modulation depth and its shot noise
— a formal uncertainty on that orientation. sigma_phi2_max says "throw away any pixel whose orientation I can't pin down to better than this many degrees." 
So sigma_phi2_max = 3 means "keep only pixels where the fibre angle is known to ±3°." It's not an intensity cut dressed up — it directly asks how trustworthy the answer is. 
And because the same signal-to-noise that blurs the angle also blurs I₂ (they're locked together, σI₂/I₂ ≈ 2·σΦ₂ in radians), that one number bounds both parameters at once. 
A bright pixel with weak modulation still fails, correctly, because a strong DC signal with a shallow polarisation response gives you a poorly-determined angle.

pool_sigma is a soft binning radius, in pixels. Rather than judging each pixel in isolation, it lets a small neighbourhood vote together before the reliability is assessed
— essentially a gentle version of hardware binning done in software. Averaging over more pixels beats down the noise (roughly as the square root of the area pooled),
so a bigger pool_sigma raises the effective SNR and lets marginal-but-real regions survive, while also smoothing the mask boundary. 
The cost is spatial resolution: you're trading the ability to resolve fine, few-pixel features for a cleaner, more contiguous, more robust mask
— which is exactly the trade you said you wanted. So pool_sigma sets how much neighbours are allowed to vouch for each other,
and sigma_phi2_max sets the reliability bar they have to clear together.

"""
import numpy as np
from scipy.ndimage import gaussian_filter, binary_closing, generate_binary_structure
try:
    from skimage.morphology import remove_small_objects, remove_small_holes, disk
    _HAVE_SK = True
except Exception:
    _HAVE_SK = False


def snr_trust_mask(res, sigma_phi2_max=6.0, pool_sigma=1.5,
                   min_feature=150, close_radius=2,
                   exclude_saturated=True, p_min=0.0):
    """
    res           : dict from pshg_gof (needs a0, R2, Phi2, gain_used,
                    read_noise2_used, calibration_factor, dof; optionally
                    saturated, pvalue)
    sigma_phi2_max: keep pixels whose pooled orientation error <= this (deg).
                    PRIMARY knob. Low-signal note: at 10-count signal the
                    per-pixel error is huge, so lean on pool_sigma; ~6-10 deg
                    is a sensible strictness once pooled.
    pool_sigma    : Gaussian pooling radius (px). Soft-bins the evidence before
                    thresholding: bigger = smoother boundary AND higher effective
                    SNR (~4*pi*s^2 pixels pooled). On weak samples this is the
                    knob that recovers structure -- push it to 3-5.
    min_feature   : remove kept regions and fill holes smaller than this (px^2).
    close_radius  : binary closing radius (px) for a smooth boundary. 0 = off.
    exclude_saturated : hard-drop any pixel that clipped in any frame.
    p_min         : optional goodness-of-fit gate (0 = off).

    Returns dict: mask, sigma_phi2_eff (deg, pooled+debiased continuous field),
    snr (pooled modulation SNR), kept_fraction.
    """
    N = res['dof'] + 5
    cal = res['calibration_factor']
    g_eff = res['gain_used'] / cal
    c_eff = res.get('read_noise2_used', 0.0) / cal            # read/baseline floor
    a0 = res['a0']
    # per-component noise variance INCLUDING the floor c_eff (key fix)
    s_amp2 = 2.0 * (g_eff * np.clip(a0, 0, None) + c_eff) / N

    R2 = res['R2']
    ang2 = np.deg2rad(2.0 * res['Phi2'])
    a2, b2 = R2 * np.cos(ang2), R2 * np.sin(ang2)

    # ---- soft spatial pooling = weighted local binning -----------------
    if pool_sigma and pool_sigma > 0:
        A = gaussian_filter(a2, pool_sigma)
        B = gaussian_filter(b2, pool_sigma)
        ksum2 = 1.0 / (4.0 * np.pi * pool_sigma**2)
        var_pool = ksum2 * gaussian_filter(s_amp2, pool_sigma / np.sqrt(2.0))
    else:
        A, B, var_pool = a2, b2, s_amp2

    sig = np.sqrt(np.clip(var_pool, 1e-12, None))
    R2p = np.hypot(A, B)
    R2p_true = np.sqrt(np.clip(R2p**2 - 2.0 * var_pool, 0.0, None))   # debias
    snr = R2p_true / sig
    with np.errstate(divide='ignore'):
        sigma_phi2_eff = np.rad2deg(0.5 * sig / R2p_true)

    keep = np.isfinite(sigma_phi2_eff) & (sigma_phi2_eff <= sigma_phi2_max)
    if exclude_saturated and 'saturated' in res:
        keep &= ~res['saturated'].astype(bool)
    if p_min > 0 and 'pvalue' in res:
        keep &= res['pvalue'] >= p_min

    keep = _clean(keep, min_feature, close_radius)
    return dict(mask=keep, sigma_phi2_eff=sigma_phi2_eff, snr=snr,
                kept_fraction=float(keep.mean()))


def _clean(mask, min_feature, close_radius):
    import warnings
    if min_feature and min_feature > 0:
        if _HAVE_SK:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', FutureWarning)
                mask = remove_small_objects(mask, int(min_feature))
                mask = remove_small_holes(mask, int(min_feature))
        else:
            mask = _remove_small_scipy(mask, int(min_feature))
            mask = ~_remove_small_scipy(~mask, int(min_feature))
    if close_radius and close_radius > 0:
        st = disk(int(close_radius)) if _HAVE_SK else generate_binary_structure(2, 1)
        mask = binary_closing(mask, structure=st)
    return mask


def _remove_small_scipy(mask, min_size):
    from scipy.ndimage import label
    lab, n = label(mask)
    if n == 0:
        return mask
    sizes = np.bincount(lab.ravel())
    keep_lab = np.where(sizes >= min_size)[0]
    keep_lab = keep_lab[keep_lab != 0]
    return np.isin(lab, keep_lab)
