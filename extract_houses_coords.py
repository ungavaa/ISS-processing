import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
import yaml
import osmnx as ox
from pyproj import Proj, transform

with open("params") as f:
	p = yaml.safe_load(f)

seuil = p['seuil']
nom_fichier = p['filename']

hdul = fits.open(nom_fichier)
im = hdul[0].data

plt.imshow(im,vmin=0,vmax=1.5)
plt.colorbar()
plt.savefig("image.png")


Y,X = np.where(im>seuil)

with open("xy.dat",'w') as f:
	for x,y in zip(X,Y):
		f.write(f"{x},{y}\n")

CRPIX1 = hdul[0].header['CRPIX1']
CDELT1 = hdul[0].header['CDELT1']
CRVAL1 = hdul[0].header['CRVAL1']
CRPIX2 = hdul[0].header['CRPIX2']
CDELT2 = hdul[0].header['CDELT2']
CRVAL2 = hdul[0].header['CRVAL2']


lons = ( X - CRPIX1 ) * CDELT1 + CRVAL1
lats = ( Y - CRPIX2 ) * CDELT2 + CRVAL2


with open("coords.dat",'w') as f:
	for lat,lon in zip(lats,lons):
		f.write(f"{lat},{lon}\n")

def street_orientation(lats, lons, srs):
    print("    Loading graph")
    Graph = ox.graph_from_bbox(
        north=max(lats) + 0.01,
        south=min(lats) - 0.01,
        east=max(lons) + 0.01,
        west=min(lons) - 0.01,
        network_type="drive",
        simplify=False,
        retain_all=True,
        truncate_by_edge=True,
        clean_periphery=True,
    )

    Graph = ox.utils_graph.get_undirected(Graph)
    Graph = ox.bearing.add_edge_bearings(Graph, precision=0)
    Graph = ox.projection.project_graph(Graph, to_crs=srs)
    nodes, edges = ox.graph_to_gdfs(Graph)
    df_routes = edges.filter(["name", "bearing", "geometry"], axis=1)

    inProj = Proj("epsg:4326")
    outProj = Proj(srs)
    X, Y = transform(inProj, outProj, lons, lats, always_xy=True)

    print("    Get nearest edges")
    edges_ID = ox.distance.nearest_edges(Graph, X, Y)
    nearest_edges = df_routes.loc[map(tuple, edges_ID)]

    print("    Compute bearings")
    coords = np.array([x.coords.xy for x in nearest_edges["geometry"]])
    lon_c, lat_c = transform(
        outProj, inProj, coords[:, 0, 0], coords[:, 1, 0], always_xy=True
    )

    distance_deg = ox.distance.euclidean_dist_vec(lat_c, lon_c, lats, lons)
    distance_m = 6371000 * np.deg2rad(distance_deg)
    
    lats_proche = lats[distance_m < p['distance']]
    lons_proche = lons[distance_m < p['distance']]
    
    with open("coords_filtered.dat",'w') as f:
	for lat,lon in zip(lats_proche,lons_proche):
		f.write(f"{lat},{lon}\n")
    
    

