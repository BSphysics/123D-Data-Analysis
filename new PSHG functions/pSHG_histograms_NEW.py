# Created on Wed Apr 22 14:04:46 2026

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator
import xlsxwriter


def orientation_stats(angles_deg):
    """
    Compute mean and spread for orientation data (0-180 degrees).
    Handles the wrap-around at 0/180 correctly.
    """
    angles_rad = np.deg2rad(2 * angles_deg)
    C = np.mean(np.cos(angles_rad))
    S = np.mean(np.sin(angles_rad))
    mean_angle = np.rad2deg(np.arctan2(S, C)) / 2
    if mean_angle < 0:
        mean_angle += 180
    R = np.sqrt(C**2 + S**2)
    circular_std = np.rad2deg(np.sqrt(-2 * np.log(R))) / 2
    return mean_angle, circular_std, R


def shift_angles_for_fitting(angles_deg):
    """
    Shift orientation angles so the distribution is centred
    away from the 0/180 wraparound boundary.
    """
    mean_angle, _, _ = orientation_stats(angles_deg)
    shift = 90 - mean_angle
    shifted = (angles_deg + shift) % 180
    return shifted, shift

from scipy.optimize import curve_fit
import xlsxwriter

# ---------------------------
# Main histogram function
# ---------------------------

def pSHGhistogramsNEW(data, mask, name, data_path):

    # Parameter-specific bins
    bin_settings = {
        'Phi2': 2,
        'I2'  : np.arange(0, 1.5, 0.02),
        'I4a' : np.arange(-0.5, 0.5, 0.01),
        'I4s' : np.arange(-0.5, 0.5, 0.01)
    }

    if name not in bin_settings:
        print("Name must be Phi2, I2, I4a or I4s")
        return None, None, data


    # Apply mask and clean data
    data = data.astype(float)
    data[mask == 0] = np.nan
    data[data == 0] = np.nan
    data = data[~np.isnan(data)]
    data = data[np.isfinite(data)]

    # fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))

    if name == 'Phi2':
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))
        bin_width = bin_settings['Phi2']
        mean_angle, circ_std, R = orientation_stats(data)
        data_shifted, shift = shift_angles_for_fitting(data)
    
        bins = np.arange(0, 180 + bin_width, bin_width)
    
        # ax1 - raw data, no shifting, straight 0-180 axis (sanity check)
        counts_raw, bin_edges_raw, _ = ax1.hist(data.ravel(), bins=bins, edgecolor='black', alpha=0.5)
        ax1.set_xlim(0, 180)
        ax1.set_title(f'Phi2 raw histogram')
    
        # ax2 - shifted data, axis relabelled back to true angles
        counts, bin_edges, _ = ax2.hist(data_shifted.ravel(), bins=bins, edgecolor='black', alpha=0.5)
        bin_centers = bin_edges[:-1] + np.diff(bin_edges) / 2
    
        tick_locs = np.arange(0, 181, 20)
        ax2.xaxis.set_major_locator(FixedLocator(tick_locs))
        tick_spacing = 20
        shift_rounded = round(shift / tick_spacing) * tick_spacing
        
        ax2.set_xticklabels([f'{int((t - shift_rounded) % 180):.0f}' for t in tick_locs])
        ax2.set_xlim(0, 180)
        ax2.set_title(f'Phi2 peak centred histogram')
    
        # Use shifted counts for Excel output
        counts = counts
        bin_edges = bin_edges
        fwhm = 2.355 * circ_std

        print(f"\nPhi2 mean = {mean_angle:.4g} degrees")
        print(f"Phi2 std  = {circ_std:.4g} degrees")
        print(f"Phi2 R (alignment) = {R:.4f}")

    else:
        fig, ax1 = plt.subplots(1, 1, figsize=(12, 4))
        bins = bin_settings[name]
        counts, bin_edges, _ = ax1.hist(data.ravel(), bins=bins, edgecolor='black', alpha=0.8, color='#FFD700')
        bin_centers = bin_edges[:-1] + np.diff(bin_edges) / 2

        # width = (bins[1] - bins[0]) * 0.9
        # ax2.bar(bin_centers, counts, width=width, alpha=0.5, edgecolor='black')

        mu = np.mean(data)
        sigma = np.std(data)
        ax1.set_title(f'{name} histogram  (mean = {mu:.3g} , std = {sigma:.3g})')
        # ax2.set_title(f'std = {sigma:.3g}')

        print(f"\n{name} mean  = {mu:.4g}")
        print(f"{name} std   = {sigma:.4g}")

    # ---------------- Save figure ----------------
    plt.savefig(f"{data_path}\\{name} histogram.png", bbox_inches='tight', pad_inches=0)

    # ---------------- Save Excel ----------------
    wb = xlsxwriter.Workbook(f"{data_path}\\{name} histogramData.xlsx")
    ws = wb.add_worksheet(name)
    ws.write_row(1, 4, ['Histogram bins', 'Histogram counts', '', 'Pixel values'])
    ws.write_column(2, 4, bin_edges[:-1].tolist())
    ws.write_column(2, 5, counts.tolist())
    ws.write_column(2, 7, data.tolist())
    wb.close()

    return bin_centers, counts, data
