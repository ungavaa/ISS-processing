#!/usr/bin/env python3

import numpy as np
import osmnx as ox
import pandas as pd
import yaml
from osgeo import gdal
from pyproj import Proj, transform

with open("params") as f:
    p = yaml.safe_load(f)

threshold = p["threshold"]
filename = "Image_Vrad.tiff"


def open_tiff(filename, dtype=np.float32):
    # Load file, and access the band and get a NumPy array
    src = gdal.Open(filename, gdal.GA_Update)
    band = src.GetRasterBand(1)
    ar = band.ReadAsArray()
    return src, ar


src, im = open_tiff(filename)

df = pd.DataFrame()
df["Y"], df["X"] = np.where(im > threshold * im.max())
df["val"] = im[df["Y"], df["X"]]

CRPIX1 = hdul[0].header["CRPIX1"]
CDELT1 = hdul[0].header["CDELT1"]
CRVAL1 = hdul[0].header["CRVAL1"]
CRPIX2 = hdul[0].header["CRPIX2"]
CDELT2 = hdul[0].header["CDELT2"]
CRVAL2 = hdul[0].header["CRVAL2"]

df["lons"] = (df["X"] - CRPIX1) * CDELT1 + CRVAL1
df["lats"] = (df["Y"] - CRPIX2) * CDELT2 + CRVAL2

if p["distance"]:
    srs = (
        "epsg:32"
        + ("6" if np.mean(df["lats"]) >= 0 else "7")
        + "%02d" % (np.mean(df["lons"]) / 6 + 31)
    )  # WGS84/UTM

    print("    Loading graph")
    Graph = ox.graph_from_bbox(
        north=max(df["lats"]) + 0.01,
        south=min(df["lats"]) - 0.01,
        east=max(df["lons"]) + 0.01,
        west=min(df["lons"]) - 0.01,
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
    X, Y = transform(inProj, outProj, df["lons"], df["lats"], always_xy=True)

    print("    Get nearest edges")
    edges_ID = ox.distance.nearest_edges(Graph, X, Y)
    nearest_edges = df_routes.loc[map(tuple, edges_ID)]

    print("    Compute distance")
    coords = np.array([x.coords.xy for x in nearest_edges["geometry"]])
    lon_c, lat_c = transform(
        outProj, inProj, coords[:, 0, 0], coords[:, 1, 0], always_xy=True
    )

    distance_deg = ox.distance.euclidean_dist_vec(
        lat_c, lon_c, df["lats"], df["lons"]
    )
    distance_m = 6.371e6 * np.deg2rad(distance_deg)

    df = df[distance_m < p["distance"]]

df.to_csv(
    "xyz.csv", columns=["lons", "lats", "val"], header=False, index=False
)
