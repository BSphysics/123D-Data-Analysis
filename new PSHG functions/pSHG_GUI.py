# -*- coding: utf-8 -*-
"""
Created on Wed May 25 15:42:25 2022
@author: Ben

Updated:
  * Two-way "masking method" toggle button:  Intensity threshold  <->  GoF/SNR trust mask.
    The parameter fields for each method show/hide with the selection.
  * Added inputs for the trust-mask parameters (sigma_phi2_max, pool_sigma,
    min_feature, close_radius).
  * The GUI now remembers every setting between runs (JSON file next to this
    script) EXCEPT the data-folder box, so batch analyses need no re-clicking.

Return (positional; original 9 items, then 5 new trailing items):
  data_path, input1, input2, input3, input4, input5, input6, threshold,
  polygonal_roi_image, mask_mode, sigma_phi2_max, pool_sigma, min_feature, close_radius
"""

import os
import json
import PySimpleGUI as sg

# --------------------------------------------------------------------------
# Persistent settings: remember everything except the folder box between runs
# --------------------------------------------------------------------------
_SETTINGS_FILE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "pSHG_gui_settings.json")

_PERSIST_KEYS = ["input1", "input2", "input3", "input4", "input5", "input6",
                 "roi_A", "roi_B", "roi_C", "thresh",
                 "sigma_phi2_max", "pool_sigma", "min_feature", "close_radius"]


def _load_settings():
    try:
        with open(_SETTINGS_FILE, "r") as f:
            return json.load(f)
    except Exception:
        return {}


def _save_settings(values, mask_mode):
    out = {k: values.get(k) for k in _PERSIST_KEYS if values and k in values}
    out["mask_mode"] = mask_mode
    try:
        with open(_SETTINGS_FILE, "w") as f:
            json.dump(out, f, indent=2)
    except Exception as e:
        print("\n Could not save GUI settings:", e)


def pSHGGUI(initialDir, defaultDir, defaultThreshold, mode):

    S = _load_settings()                         # last-used values (if any)

    def g(key, fallback):                        # last value, else fallback
        v = S.get(key, fallback)
        return fallback if v is None else v

    if mode.lower() == 'transmission':
        sg.theme("DarkGreen5")
    else:
        sg.theme("DarkBlue15")

    mask_mode = g("mask_mode", "intensity")      # "intensity" or "trust"

    def mode_button_text(m):
        a = ("\u25B6 " if m == "intensity" else "   ") + "Intensity threshold"
        b = ("\u25B6 " if m == "trust" else "   ") + "GoF / SNR trust mask"
        return a + "      |      " + b + "      (click to switch)"

    # ROI radio defaults (exactly one True)
    roiA, roiB, roiC = bool(g("roi_A", True)), bool(g("roi_B", False)), bool(g("roi_C", False))
    if not (roiA or roiB or roiC):
        roiA = True

    file_list_column = [
        [sg.Text("Select data folder"),
         sg.In(size=(25, 5), enable_events=True, key="folder"),
         sg.FolderBrowse(initial_folder=initialDir)],
    ]

    # --- intensity-threshold parameter group (shown in intensity mode) ---
    intensity_group = sg.Column([
        [sg.Checkbox('Use slider intensity threshold?',
                     default=bool(g("input6", True)), key="input6")],
        [sg.Text('Intensity threshold', size=(18, 1)),
         sg.InputText(g("thresh", ""), size=(10, 1), enable_events=True, key="thresh")],
    ], key="-INT_GROUP-", visible=(mask_mode == "intensity"), pad=(0, 0))

    # --- trust-mask parameter group (shown in trust mode) ---
    trust_group = sg.Column([
        [sg.Text('sigma_phi2_max (deg)', size=(18, 1)),
         sg.InputText(g("sigma_phi2_max", "3.0"), size=(8, 1), key="sigma_phi2_max")],
        [sg.Text('pool_sigma (px)', size=(18, 1)),
         sg.InputText(g("pool_sigma", "1.5"), size=(8, 1), key="pool_sigma")],
        [sg.Text('min_feature (px^2)', size=(18, 1)),
         sg.InputText(g("min_feature", "150"), size=(8, 1), key="min_feature")],
        [sg.Text('close_radius (px)', size=(18, 1)),
         sg.InputText(g("close_radius", "2"), size=(8, 1), key="close_radius")],
    ], key="-TRUST_GROUP-", visible=(mask_mode == "trust"), pad=(0, 0))

    config_column = [
        [sg.Text("Select options for the pSHG analysis")],
        [sg.Checkbox('Whole image histograms?', default=bool(g("input1", False)), key="input1")],
        [sg.Checkbox('Show single pixel fits?', default=bool(g("input2", False)), key="input2")],
        [sg.Checkbox('Polar histogram plot?', default=bool(g("input3", False)), key="input3")],
        [sg.Checkbox('Show Arrow plot?', default=bool(g("input4", True)), key="input4")],
        [sg.Checkbox('Use polygonal ROI?', default=bool(g("input5", False)), key="input5")],
        [sg.Radio('Phi 2', "roi_mode", default=roiA, key="roi_A"),
         sg.Radio('I2', "roi_mode", default=roiB, key="roi_B"),
         sg.Radio('Total SHG', "roi_mode", default=roiC, key="roi_C")],
        [sg.HorizontalSeparator()],
        [sg.Text("Masking method :")],
        [sg.Button(mode_button_text(mask_mode), key="-MODE-")],
        [sg.pin(intensity_group), sg.pin(trust_group)],
        [sg.HorizontalSeparator()],
        [sg.Button("Run")],
    ]

    layout = [[sg.Column(file_list_column),
               sg.VSeperator(),
               sg.Column(config_column)]]

    window = sg.Window("pSHG setup", layout)

    # Event loop
    while True:
        event, values = window.read()
        if event == sg.WIN_CLOSED or event == "Run":
            break
        if event == "-MODE-":                    # flip the two-way switch
            mask_mode = "trust" if mask_mode == "intensity" else "intensity"
            window["-MODE-"].update(mode_button_text(mask_mode))
            window["-INT_GROUP-"].update(visible=(mask_mode == "intensity"))
            window["-TRUST_GROUP-"].update(visible=(mask_mode == "trust"))

    if values is None:                           # window closed with the X
        values = {k: S.get(k) for k in _PERSIST_KEYS}

    _save_settings(values, mask_mode)            # remember for next time
    window.close()

    # ---- resolve directory ----
    data_path = values.get("folder")
    if not data_path:
        data_path = defaultDir
        print('\n Using default directory (see below)\n' + str(data_path))
    else:
        print('\n Using user defined directory (see below)\n' + str(data_path))

    # ---- resolve intensity threshold (unchanged logic) ----
    thresh = values.get("thresh")
    if not thresh and values.get("input6") == False:
        threshold = defaultThreshold
        print('\n Using default threshold i.e. ' + str(threshold))
    elif not thresh and values.get("input6") == True:
        threshold = 0.0
        print('\n Using slider threshold')
    else:
        try:
            threshold = float(thresh)
            print('\n Using user defined threshold i.e. ' + str(threshold))
        except (TypeError, ValueError):
            threshold = defaultThreshold
            print('\n Threshold unreadable, using default i.e. ' + str(threshold))

    # ---- ROI image choice ----
    if values.get("roi_A"):
        polygonal_roi_image = "Phi2"
    elif values.get("roi_B"):
        polygonal_roi_image = "I2"
    elif values.get("roi_C"):
        polygonal_roi_image = "All SHG"
    else:
        polygonal_roi_image = None

    # ---- trust-mask parameters (safe fallbacks) ----
    def _num(key, cast, default):
        try:
            return cast(float(values.get(key)))
        except (TypeError, ValueError):
            return default
    sigma_phi2_max = _num("sigma_phi2_max", float, 3.0)
    pool_sigma = _num("pool_sigma", float, 1.5)
    min_feature = _num("min_feature", int, 150)
    close_radius = _num("close_radius", int, 2)

    print('\n Masking method: ' + ('GoF / SNR trust mask' if mask_mode == "trust"
                                    else 'Intensity threshold'))
    if mask_mode == "trust":
        print('   sigma_phi2_max={}, pool_sigma={}, min_feature={}, close_radius={}'.format(
            sigma_phi2_max, pool_sigma, min_feature, close_radius))

    if values.get("input1"):
        print('\n Plotting whole image histograms')
    if values.get("input2"):
        print('\n Use the pSHG fit viewer')
    if values.get("input3"):
        print('\n Plotting the polar histogram')

    return (data_path, values.get("input1"), values.get("input2"),
            values.get("input3"), values.get("input4"), values.get("input5"),
            values.get("input6"), threshold, polygonal_roi_image,
            mask_mode, sigma_phi2_max, pool_sigma, min_feature, close_radius)
