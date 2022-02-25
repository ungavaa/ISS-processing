# ******************************************************************************
#                      Make discrete inventory from ISS images
#
# Auteur : Julien-Pierre Houle
# Date :   June 2021
# ******************************************************************************

import os

import numpy as np
import pandas as pd
import pytools as pt
import yaml

with open("params") as f:
    p = yaml.safe_load(f)

# PARAMETERS
focal_dist = p["focal_dist"]  # focal distance (nm)
obs_angle = p["obs_angle"]  # observer_angle
vband_width = p["vband_width"]  # nm
inventory_name = p["inventory_name"]  # filename of the discrete inventory

# PATH to the arrays load from extract_from_tif.py
PATH_ARRAYS = p["PATH_ARRAY"]
PATH_ZONES_INVENTORY = p["PATH_ZONES_INVENTORY"]
# PATH to a directory containing the different class spectrum
PATH_SPCT = p["PATH_SPCT"]
# Path to the LIGHTS directory produce by Illumina when executing illum inputs.
PATH_LIGHTS = p["PATH_LIGHTS"]

# Note: The LIGTHS directory contain the spct and lop files. You should add the
# desired lop file in the LIGTHS directory if the lop value wanted is missing.

# Loading np arrays files
files = {
    f: np.load(f"{PATH_ARRAYS}/np_{f}.npy")
    for f in ["domain", "coord", "intensity", "tech"]
}


print("Calculating pow..")
abc = pd.read_csv(["spectrum"])
spectrum = abc["class_spct"]
Uplight = abc["lop"]
tech_equiv = {
    0: [np.nan, 0],  # spct | lop
    1: ["ClassI", "5"],
    2: ["ClassII", "5"],
    3: ["ClassIII", "15"],
    4: ["ClassIV", "1"],
    5: ["ClassV", "1"],
}  # Modification 2 pour 1 (car pas de 2.lop dans Lights)

# Arrays valeur pour chaque pixel
arr_spcts = np.array([tech_equiv[v][0] for v in files["tech"]])
arr_lops = np.array([tech_equiv[v][1] for v in files["tech"]])
keys = list(zip(arr_lops, arr_spcts))


wav, norm_spectrum = np.loadtxt(
    f"{PATH_LIGHTS}/Lights/photopic.dat", skiprows=1
).T
sptc = np.loadtxt(f"{PATH_LIGHTS}/Lights/scotopic.dat", skiprows=1)
wl, refl = np.loadtxt(f"{PATH_LIGHTS}/Lights/asphalt.aster").T
wl *= 1000.0
refl /= 100.0
asphalt = np.interp(wav, wl, refl)
angles = np.arange(181, dtype=float)

# Load pre-normalize spcts & lops
list_spcts = os.listdir(PATH_SPCT)
spcts = {
    spct.split("_")[0]: pt.load_spct(wav, norm_spectrum, f"{PATH_SPCT}/{spct}")
    for spct in list_spcts
}
lops = {
    lop: pt.load_lop(angles, f"{PATH_LIGHTS}/Lights/{lop}_pcUPLIGHT.lop")
    for lop in ["1", "5", "15"]
}

S = ((408000 * 8.4e-6) / focal_dist) ** 2  # Taille du pixel au sol (unit/??)
a = np.deg2rad(angles)
mids = np.concatenate([[a[0]], np.mean([a[1:], a[:-1]], 0), [a[-1]]])
sinx = 2 * np.pi * (np.cos(mids[:-1]) - np.cos(mids[1:]))


# Calculating integral for each lop-spct combinaisons
# Atmospheric correction pre-calculated by Alejandro
# Intensity: (nw/sr/cm^2/angstrom) We suppose const value on photopic bandwidth
# Formula: I = DNB * S / integral( R(lambd) * T(lambd) *
# (1/pi* p(lambd)* F(lambd) + G(lambd))) dlambd)
integral = dict()
for lop in lops:
    F = np.sum(lops[lop][angles > 90] * sinx[angles > 90])
    for spct in spcts:
        integral[lop, spct] = S / (
            (
                sum(
                    spcts[spct]
                    * (1 / np.pi * asphalt * F + lops[lop][31])
                    * norm_spectrum
                )
            )
            * (wl[1] - wl[0])
        )


# Attributing POW to pixels
POW = np.zeros(len(keys))
for idx in range(len(keys)):
    POW[idx] = (
        integral[keys[idx]] * files["intensity"][idx] * 1e-5 * vband_width * 10
    )  # factor 10 to convert nm to angstrom


print("Creating discrete inventory..")
df = pd.DataFrame(files["coord"], columns=["#lat", "lon"])
df["pow"] = POW.tolist()
df["hobs"] = 0
df["dobs"] = 0
df["fobs"] = 0
df["hlamp"] = 0
df["spct"] = arr_spcts.tolist()
df["lop"] = arr_lops.tolist()
df.to_csv(inventory_name, index=False, sep=" ")
print("Done")
