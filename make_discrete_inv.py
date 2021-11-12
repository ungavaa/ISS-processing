# ******************************************************************************
#                      Make discrete inventory from ISS images
#
# Auteur : Julien-Pierre Houle
# Date :   June 2021
# ******************************************************************************

import os
import sys
import numpy as np
import pandas as pd
import progressbar
sys.path.append("/home/git/illumina") # path to Illumina model
import pyproj
import pytools as pt


# PARAMETERS
focal_dist = 0.4   # focal distance (nm)
obs_angle = 31     # observer_angle
vband_width = 100  # nm
inventory_name = "discrete_inventory.txt"   # filename of the discrete inventory

PATH_ARRAYS = "sources"   # PATH to the arrays load from extract_from_tif.py
PATH_ZONES_INVENTORY = "zones_inventory.csv"
PATH_SPCT = "Datas/spct"   # PATH to a directory containing the different class spectrum
PATH_LIGHTS = "."          # Path to the LIGHTS directory produce by Illumina when executing illum inputs.

# Note: The LIGTHS directory contain the spct and lop files. You should 
# add the desired lop file in the LIGTHS directory if the lop value wanted is missing. 

# Loading np arrays files
files = { f: np.load(f'{PATH_ARRAYS}/np_{f}.npy') for f in
                        ['domain', 'coord', 'intensity', 'tech'] }

# Read values from zones inventory
df_zones = pd.read_csv(PATH_ZONES_INVENTORY)
df_zones = df_zones.sort_values('radius', ascending=False).reset_index(drop=True)
zones_p = { p: df_zones[f'{p}'].to_numpy() for p in
                    ['lat', 'lon', 'radius', 'do', 'ho', 'fo', 'hl'] }


# Convert lat, lon zones to meters (2949)
proj = pyproj.Transformer.from_crs( 4326, 2949, always_xy=False )
i, j = proj.transform( zones_p['lat'], zones_p['lon'] )
zones_center = np.array((i, j)).T
zones = np.zeros((len(files['domain'])))
hobs  = np.zeros((len(files['domain'])))
dobs  = np.zeros((len(files['domain'])))
fobs  = np.zeros((len(files['domain'])))
hlamp = np.zeros((len(files['domain'])))

print("Create circulars zones masks..")
for idx, center in enumerate(zones_center):
    dist_from_center = np.sqrt((files['domain'][:,0] - center[0])**2 + (files['domain'][:,1] - center[1])**2)
    mask = dist_from_center <= zones_p['radius'][idx]

    zones[mask] = idx
    hobs[zones==idx] = zones_p['ho'][idx]
    dobs[zones==idx] = zones_p['do'][idx]
    fobs[zones==idx] = zones_p['fo'][idx]
    hlamp[zones==idx] = zones_p['hl'][idx]


print('Calculating pow..')
tech_equiv = {0: [np.nan,      0],   # spct | lop
              1: ['ClassI',   '5'],
              2: ['ClassII',  '5'],
              3: ['ClassIII', '15'],
              4: ['ClassIV',  '1'],
              5: ['ClassV',   '1']}  # Modification 2 pour 1 (car pas de 2.lop dans Lights)

# Arrays valeur pour chaque pixel
arr_spcts= np.array([tech_equiv[v][0] for v in files['tech']])
arr_lops = np.array([tech_equiv[v][1] for v in files['tech']])
keys = list(zip(arr_lops, arr_spcts))


wav, norm_spectrum = np.loadtxt(f"{PATH_LIGHTS}/Lights/photopic.dat", skiprows=1).T
sptc = np.loadtxt(f"{PATH_LIGHTS}/Lights/scotopic.dat", skiprows=1)
wl, refl = np.loadtxt(f'{PATH_LIGHTS}/Lights/asphalt.aster').T
wl *= 1000.
refl /= 100.
asphalt = np.interp(wav, wl, refl)
angles = np.arange(181, dtype=float)

# Load pre-normalize spcts & lops
list_spcts = os.listdir(PATH_SPCT)
spcts = {spct.split("_")[0] : pt.load_spct(wav, norm_spectrum, f'{PATH_SPCT}/{spct}') for spct in list_spcts}
lops = {lop: pt.load_lop(angles, f'{PATH_LIGHTS}/Lights/{lop}_pcUPLIGHT.lop') for lop in ['1','5','15']}

S = ((408000 * 8.4e-6) / focal_dist)**2  # Taille du pixel au sol (unit/??)
a = np.deg2rad(angles)
mids = np.concatenate([[a[0]],np.mean([a[1:],a[:-1]],0),[a[-1]]])
sinx = 2*np.pi*(np.cos(mids[:-1])-np.cos(mids[1:]))


# Calculating integral for each lop-spct combinaisons
# Atmospheric correction pre-calculated by Alejandro
# Intensity: (nw/sr/cm^2/angstrom) - We suppose constant value on photopic band width
# Formula: I = DNB * S / integral( R(lambd) * T(lambd) * (1/pi* p(lambd)* F(lambd) + G(lambd))) dlambd)
integral = dict()
for lop in lops:
    F = np.sum(lops[lop][angles>90] * sinx[angles>90])
    for spct in spcts:
        integral[lop,spct] = S / ((sum(spcts[spct] *( 1/np.pi\
                                * asphalt * F + lops[lop][31]) * norm_spectrum)) * (wl[1]-wl[0]))


# Attributing POW to pixels
POW = np.zeros(len(keys))
for idx in range(len(keys)):
    POW[idx] = integral[keys[idx]] * files['intensity'][idx] * 1e-5 * vband_width * 10 # factor 10 to convert nm to angstrom


print("Creating discrete inventory..")
df = pd.DataFrame(files['coord'], columns = ['#lat', 'lon'])
df['pow'] = POW.tolist()
df['hobs'] = hobs.tolist()
df['dobs'] = dobs.tolist()
df['fobs'] = fobs.tolist()
df['hlamp'] = hlamp.tolist()
df['spct'] = arr_spcts.tolist()
df['lop'] = arr_lops.tolist()
df.to_csv(inventory_name, index=False, sep=' ')
print('Done')
