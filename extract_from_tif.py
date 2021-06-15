# ******************************************************************************
#                         Convert georef file to np.arrays
#
# Auteur : Julien-Pierre Houle
# Date :   June 2021
# ******************************************************************************

import gdal
from geopy import distance
import pyproj, osr
import numpy as np
import pandas as pd
from osgeo import gdal
import os

# Create XYZ file from raster (use in terminal)
# gdal_translate -of XYZ PATH_RASTER.tiff PATH_XYZ.csv

PATH_RASTER = "Datas/Iss_img/Georef" # PATH to directory containing raster files
PATH_XYZ    = "Datas/Iss_img/XYZ"    # PATH to directory containing XYZ.csv files
PATH_ARRAYS = "sources/np_arrays"    # PATH to directory where np arrays will be write

params = dict()
files = dict()

print("Extracting values from geotiff..")
for idx, img in enumerate(os.listdir(PATH_RASTER)):

    raster = gdal.Open(f'{PATH_RASTER}/{img}')
    x, y, val = np.loadtxt(f'{PATH_XYZ}/{img[:-4]}.csv', delimiter=' ').T
    val[val == 0.] = np.nan  # convert 0 to nan

    # Extracting epsg from tiff
    proj = osr.SpatialReference(wkt=raster.GetProjection())
    init_proj = proj.GetAttrValue('AUTHORITY',1)
    proj = pyproj.Transformer.from_crs(int(init_proj), 4326, always_xy=False)
    lat, lon = proj.transform(x, y)

    params.update({ img.split('_')[0]: (np.around(lat, 5), np.around(lon, 5), val)})


# Create mask to remove nan
mask  = ~np.isnan(params['tech'][2]) &  ~np.isnan(params['intensity'][2])
lat   = params['tech'][0][mask]
lon   = params['tech'][1][mask]
files.update({ 'tech': params['tech'][2][mask] })
files.update({ 'intensity': params['intensity'][2][mask] })

print("Create domain image..")
proj_dom = pyproj.Transformer.from_crs(4326, 2949, always_xy=False)
i, j = proj_dom.transform(lat, lon)  # convert to meters
files.update({ 'domain': np.array((i, j)).T})
files.update({ 'coord': np.array((lat, lon)).T})


# Save np arrays
for arr in files:
    np.save(f'{PATH_ARRAYS}/np_{arr}', files[arr])
