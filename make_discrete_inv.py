# ******************************************************************************
#                      Make discrete inventory from ISS images
#
# Auteur : Julien-Pierre Houle
# Date :   Avril 2021
# ******************************************************************************

import sys
import pandas as pd
import progressbar
import pytools
import pyproj
sys.path.append("/home/jhoule42/git/illumina")
import pytools as pt


# Parameters
focal_dist = 0.4   # focal distance (nm)
obs_angle = 31
vband_width = 100  # nm


# Loading np arrays files
files = { f: np.load(f'Intrusif/np_arrays/np_{f}.npy') for f in
                    ['domain', 'coord', 'intensity', 'tech'] }

# Extract values from zones
df_zones = pd.read_csv('Intrusif/zone_obst.csv')
df_zones = df_zones.sort_values('radius', ascending=False).reset_index()
zones_p = { p: df_zones[f'{p}'].to_numpy() for p in
                    ['lat', 'lon', 'radius', 'hobs', 'dobs'] }

# Convert to meters (2949)
proj = pyproj.Transformer.from_crs( 4326, 2949, always_xy=False )
i, j = proj.transform( zones_p['lat'], zones_p['lon'] )
zones_center = np.array((i, j)).T
zones = np.zeros((len(files['domain'])))
hobs =  np.zeros((len(files['domain'])))
dobs =  np.zeros((len(files['domain'])))


print("Create circulars zones masks..")
for idx, center in enumerate(zones_center):
    dist_from_center = np.sqrt((files['domain'][:,0] - center[0])**2 + (files['domain'][:,1] - center[1])**2)
    mask = dist_from_center <= zones_p['radius'][idx]*1000

    zones[mask] = idx
    hobs[zones==idx] = zones_p['hobs'][idx]
    dobs[zones==idx] = zones_p['dobs'][idx]


print('Calculating pow..')
tech_equiv = {0: [np.nan,    0,  0], # spct | hlamp | lop
              1: ['Class1',  8, '5'],
              2: ['Class2',  8, '5'],
              3: ['Class3',  2, '15'],
              4: ['Class4',  8, '2'],
              5: ['Class5',  8, '2']}

# Arrays valeur pour chaque pixel
arr_spcts= np.array([tech_equiv[v][0] for v in files['tech']])
hlamp    = np.array([tech_equiv[v][1] for v in files['tech']])
arr_lops = np.array([tech_equiv[v][2] for v in files['tech']])
keys = list(zip(arr_lops, arr_spcts))


wav, norm_spectrum = np.loadtxt("Validation_Direct/Lights/photopic.dat", skiprows=1).T
sptc = np.loadtxt('Intrusif/Lights/scotopic.dat', skiprows=1)
wl, refl = np.loadtxt('Intrusif/Lights/asphalt.aster').T
wl *= 1000.
refl /= 100.
asphalt = np.interp(wav, wl, refl)
angles = np.arange(181, dtype=float)

# Load pre-normalize spcts & lops
list_spcts = os.listdir('Datas/spct')
spcts = {spct[:6]: pt.load_spct(wav, norm_spectrum, f'Datas/spct/{spct}') for spct in list_spcts}
lops = {lop: pt.load_lop(angles, f'Intrusif/Lights/{lop}_pcUPLIGHT.lop') for lop in ['2','5','15']}

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
bar = progressbar.ProgressBar(maxval=len(keys)).start()
for idx in range(len(keys)):
    POW[idx] = integral[keys[idx]] * files['intensity'][idx] * 1e-5 * vband_width * 10 # factor 10 to convert nm to angstrom
    bar.update(idx)



print("\nCreating discrete inventory..")
df = pd.DataFrame(files['coord'], columns = ['lat', 'lon'])
df['pow'] = POW.tolist()
df['hobs'] = hobs.tolist()
df['dobs'] = dobs.tolist()
df['fobs'] = 0.7
df['hlamp'] = hlamp.tolist()
df['spct'] = arr_spcts.tolist()
df['lop'] = arr_lops.tolist()
df.to_csv('Intrusif/discrete_inventory.txt', index=False)
