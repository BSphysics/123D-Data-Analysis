# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 15:42:52 2026

@author: bs426
"""

import os
scriptDir = os.getcwd()
import sys
sys.path.append(os.path.join(scriptDir,"new pSHG functions" ))
import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt


def density_plot(I4a, I4s, I2=None, draw_arc=False, ratio=1.4,
                 rho_for_arc=(30, 60), out_path="pshg_density.png"):
    """
    I4a, I4s, (I2): 2-D arrays straight from your .npy files.
    Background pixels are assumed stored as exactly 0 and are masked out.
 
    Returns
    -------
    fig      : the matplotlib Figure
    centroid : (cx, cy) tuple = density-weighted centroid in (I4a, I4s) space
    """
    # --- masking: keep only real signal -------------------------------
    if I2 is not None:
        valid = (I2 > 0.02) & (np.abs(I4s) > 1e-6)
    else:
        valid = (np.abs(I4a) > 1e-6) | (np.abs(I4s) > 1e-6)
 
    a = I4a[valid].ravel()
    s = I4s[valid].ravel()
    print(f"{valid.sum()} valid pixels of {I4a.size} "
          f"({100 * valid.mean():.1f}%)")
    print(f"  I4a median = {np.median(a):+.3f}   "
          f"frac>0 = {np.mean(a > 0):.2f}")
    print(f"  I4s median = {np.median(s):+.3f}")
 
    fig, ax = plt.subplots(figsize=(6.5, 6))
 
    # --- the density itself -------------------------------------------
    hb = ax.hexbin(a, s, gridsize=60, cmap="inferno",
                   bins="log", mincnt=1)
    fig.colorbar(hb, ax=ax, label="log$_{10}$ pixel count")
 
    # --- density-weighted centroid ------------------------------------
    # The density-weighted centroid in (I4a, I4s) space, weighting every
    # hexagon by how many pixels fall in it, is mathematically identical to
    # the plain mean of the underlying pixel values -- so compute it directly
    # from the raw pixels. (Doing it via hexbin.get_array() is fragile: with
    # bins='log' different matplotlib versions return either log10(count) or
    # raw counts, which silently corrupts the weighting.)
    #
    # NOTE: this is the MEAN, which is what "density-weighted centroid" means.
    # It is pulled toward the dense core but still feels the sparse wings. If
    # you want a value robust to outlier wings, swap to np.median below.
    cx = float(np.median(a))
    cy = float(np.median(s))
    centroid = (cx, cy)
    print(f"  density-weighted centroid (mean): I4a={cx:+.4f}  I4s={cy:+.4f} \n")
 
    # plot it as a green disc
    ax.plot(cx, cy, "o", color="lime", mec="k", mew=1.2, ms=14,
            zorder=5, label="density centroid ")
    ax.annotate(f"({cx:+.3f}, {cy:+.3f})", (cx, cy),
                textcoords="offset points", xytext=(10, 8),
                fontsize=9, color="k",
                bbox=dict(boxstyle="round,pad=0.2", fc="lime", ec="k", alpha=0.7))
 
    # --- model reference arcs (optional) ------------------------------
    # if draw_arc:
    #     for rho, col in zip(rho_for_arc, ("cyan", "lime")):
    #         xa, xs = birefringence_arc(ratio, rho)
    #         ax.plot(xa, xs, "-", color=col, lw=2.2,
    #                 label=fr"biref. arc $\rho={rho}°$")
    #     ax.plot(0, c6_response(ratio, 60, 0)[2], "o",
    #             color="white", mec="k", ms=9, label="0$\\lambda$ baseline")
 
    ax.axvline(0, color="0.6", lw=0.6)
    ax.set_xlabel(r"$I_{4a}$")
    ax.set_ylabel(r"$I_{4s}$")
    ax.set_xlim(-0.4, 0.4)
    ax.set_ylim(-0.45, 0.30)
    ax.set_title("$I_{4a}$ vs $I_{4s}$")
    ax.legend(fontsize=8, loc="lower left")
 
    fig.tight_layout()
    fig.savefig(out_path, dpi=350)
    # print(f"saved -> {out_path}")
    return fig, centroid
