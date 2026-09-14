# -*- coding: utf-8 -*-
"""
snr_mask.py -- robust, tunable, smoothed trust mask for pSHG I2/phi2.

Operates on the dict returned by pshg_gof(). The reliability metric is the
propagated orientation error sigma(phi2), which also bounds the I2 error
(sigma_I2/I2 = 2*sigma_phi2[rad]). Because it is in physical units (degrees)
it is invariant to signal level -> the same threshold transfers between
acquisitions, unlike a raw-intensity cut.

Pipeline:  per-pixel sigma(phi2) with debiased R2   (noise-rectification removed)
        -> soft spatial pooling (weighted binning, tunable radius)
        -> threshold at sigma_phi2_max
        -> drop small islands / fill small holes / optional edge smoothing

Tune once:   sigma_phi2_max (strictness), pool_sigma (pooling/smoothness)
Set&forget:  min_feature (smallest region/hole you care about), close_radius

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
                    calibration_factor, dof; optionally saturated, pvalue)
    sigma_phi2_max: keep pixels whose pooled orientation error <= this (deg).
                    PRIMARY knob. ~3 strict, ~6 moderate, ~10 permissive.
                    Equivalent I2 bound: sigma_I2/I2 <= 2*sigma_phi2_max[rad].
    pool_sigma    : Gaussian pooling radius (px) for soft binning of the
                    evidence before thresholding. 0 = per-pixel. Larger =
                    smoother boundary AND higher effective SNR (bins ~4*pi*s^2
                    pixels), so you can raise the strictness when you raise this.
    min_feature   : remove kept regions and fill holes smaller than this (px^2).
    close_radius  : binary closing radius (px) for a smooth boundary. 0 = off.
    exclude_saturated : hard-drop any pixel that clipped in any frame.
    p_min         : optional goodness-of-fit gate (0 = off). e.g. 1e-3 also
                    removes model-violating pixels (motion, mixed fibres).

    Returns dict: mask (bool), sigma_phi2_eff (deg, pooled/debiased continuous
    field for your own re-thresholding), snr (pooled modulation SNR),
    kept_fraction.
    """
    N = res['dof'] + 5
    g_eff = res['gain_used'] / res['calibration_factor']
    a0 = np.clip(res['a0'], 1.0, None)
    s_amp2 = 2.0 * g_eff * a0 / N                      # per-component noise var
    R2 = res['R2']
    ang2 = np.deg2rad(2.0 * res['Phi2'])
    a2, b2 = R2 * np.cos(ang2), R2 * np.sin(ang2)      # recover Fourier coeffs

    # ---- soft spatial pooling = weighted local binning -----------------
    if pool_sigma and pool_sigma > 0:
        A = gaussian_filter(a2, pool_sigma)
        B = gaussian_filter(b2, pool_sigma)
        # Var of a normalised-Gaussian-weighted mean: sum(g_i^2) * <var>,
        # and normalised g^2 is a Gaussian of sigma pool_sigma/sqrt(2).
        ksum2 = 1.0 / (4.0 * np.pi * pool_sigma**2)
        var_pool = ksum2 * gaussian_filter(s_amp2, pool_sigma / np.sqrt(2.0))
    else:
        A, B, var_pool = a2, b2, s_amp2

    sig = np.sqrt(np.clip(var_pool, 1e-12, None))
    R2p = np.hypot(A, B)
    R2p_true = np.sqrt(np.clip(R2p**2 - 2.0 * var_pool, 0.0, None))  # debias
    snr = R2p_true / sig
    with np.errstate(divide='ignore'):
        sigma_phi2_eff = np.rad2deg(0.5 * sig / R2p_true)   # -> inf on noise-only

    # ---- reliability threshold (bounds phi2 AND I2 error) --------------
    keep = np.isfinite(sigma_phi2_eff) & (sigma_phi2_eff <= sigma_phi2_max)
    if exclude_saturated and 'saturated' in res:
        keep &= ~res['saturated'].astype(bool)
    if p_min > 0 and 'pvalue' in res:
        keep &= res['pvalue'] >= p_min

    # ---- morphological cleanup: chunky, contiguous regions -------------
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
        if _HAVE_SK:
            st = disk(int(close_radius))
        else:
            st = generate_binary_structure(2, 1)
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
