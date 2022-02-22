#!/usr/bin/env python3

import sys

import astropy.io.fits as pyfits
import numpy
import numpy as np
from osgeo import gdal

filename = sys.argv[1]

hdulist = pyfits.open(filename)
other = np.asarray(hdulist[0].data)
earth = other.copy()
earth = np.float32(earth)
nrows, ncols = earth.shape[0], earth.shape[1]
dst_ds = gdal.GetDriverByName("GTiff").Create(
    filename[:-5] + ".tiff", ncols, nrows, 1, gdal.GDT_Float32
)
dst_ds.SetProjection("epsg:4326")
dst_ds.SetGeoTransform(
    (
        hdulist[0].header["CRVAL1"],
        hdulist[0].header["CRPIX1"],
        0,
        hdulist[0].header["CRVAL2"],
        0,
        hdulist[0].header["CRPIX2"],
    )
)
dst_ds.GetRasterBand(1).WriteArray(earth)
