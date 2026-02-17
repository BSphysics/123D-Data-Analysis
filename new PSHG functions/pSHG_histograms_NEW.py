# -*- coding: utf-8 -*-
"""
Created on Fri Feb  6 14:04:46 2026

@author: bs426
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import xlsxwriter


# ---------------------------
# Gaussian model
# ---------------------------

def gaussian(x, a, mu, sigma, offset=0):
    return a * np.exp(-(x - mu)**2 / (2 * sigma**2)) + offset


# ---------------------------
# Robust Gaussian fit helper
# ---------------------------

def fit_gaussian(bin_centers, counts):
    
    mask = counts > 0
    x = bin_centers[mask]
    y = counts[mask]

    if len(x) < 5:
        return None   # not enough data to fit

    # initial guesses from moments
    mu0 = np.average(x, weights=y)
    sigma0 = np.sqrt(np.average((x - mu0)**2, weights=y))
    a0 = np.max(y)
    offset0 = np.min(y)

    p0 = [a0, mu0, sigma0, offset0]

    try:
        popt, _ = curve_fit(gaussian, x, y, p0=p0, maxfev=5000)
        return popt
    except RuntimeError:
        return None


# ---------------------------
# Main histogram function
# ---------------------------

def pSHGhistogramsNEW(data, mask, name, data_path):

    # Apply mask safely
    data = data.astype(float)
    data[mask == 0] = np.nan
    data = data[~np.isnan(data)]

    #print(f"\nNumber of pixels = {len(data)}")

    # Parameter-specific bins
    bin_settings = {
        'Phi2': np.arange(0, 180, 2),
        'I2'  : np.arange(0, 1.5, 0.02),
        'I4a' : np.arange(-0.5, 0.5, 0.01),
        'I4s' : np.arange(-0.5, 0.5, 0.01)
    }

    if name not in bin_settings:
        print("Name must be Phi2, I2, I4a or I4s")
        return None, None, data

    bins = bin_settings[name]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12,4))

    counts, bins, _ = ax1.hist(
        data.ravel(),
        bins=bins,
        edgecolor='black',
        alpha=0.5
    )

    bin_centers = bins[:-1] + np.diff(bins)/2
    ax1.set_title(f"{name} histogram")

    # Plot lighter version for fit overlay
    width = (bins[1]-bins[0])*0.9
    ax2.bar(bin_centers, counts, width=width, alpha=0.2, edgecolor='black')

    # ---- Gaussian fit ----
    popt = fit_gaussian(bin_centers, counts)

    mu = np.mean(data)
    sigma = np.std(data)
    fwhm_moment = 2.355 * sigma

    if popt is not None:

        xfit = np.linspace(bins[0], bins[-1], 1000)
        yfit = gaussian(xfit, *popt)
        ax2.plot(xfit, yfit, '--r')

        fwhm_fit = 2.355 * abs(popt[2])
        centre_fit = popt[1]

        ax2.set_title(
            f"FWHM = {fwhm_fit:.3g} , Centre = {centre_fit:.3g}"
        )

        print(f"\n\n{name} Gaussian centre = {centre_fit:.4g}")
        print(f"{name} Gaussian FWHM   = {fwhm_fit:.4g}")

    else:
        ax2.set_title("Gaussian fit failed (using moments)")
        print(f"{name} Gaussian fit failed")

    print(f"{name} mean = {mu:.4g} ± {sigma:.4g}")

    # ---------------- Save figure ----------------

    plt.savefig(f"{data_path}\\{name} histogram.png",
                bbox_inches='tight', pad_inches=0)

    # ---------------- Save Excel ----------------

    wb = xlsxwriter.Workbook(f"{data_path}\\{name} histogramData.xlsx")
    ws = wb.add_worksheet(name)

    ws.write_row(1,4, ['Histogram bins','Histogram counts', '', 'Pixel values'])

    ws.write_column(2,4, bins[:-1])
    ws.write_column(2,5, counts)
    ws.write_column(2,7, data)

    wb.close()

    return bin_centers, counts, data
