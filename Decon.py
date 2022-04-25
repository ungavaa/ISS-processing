#!/usr/bin/env python3

import numpy as np
from osgeo import gdal
from skimage import restoration

print("L'opération peut prendre jusqu'à 30 minutes...")


def open_tiff(filename, dtype=np.float32):
    # Load file, and access the band and get a NumPy array
    src = gdal.Open(filename, gdal.GA_Update)
    band = src.GetRasterBand(1)
    ar = band.ReadAsArray()
    return src, ar


src, psf = open_tiff("psf.tif")
src, ville = open_tiff("Vrad.tiff")


psf1 = psf / np.sum(psf)

ville[np.isnan(ville)] = 0
psf1[np.isnan(psf1)] = 0
ville += ville.max() * 1e-5 * np.random.random(ville.shape)
ville += ville.max() * 1e-5 * np.random.random(ville.shape)

deconvolved = restoration.richardson_lucy(ville, psf1, 1000, clip=False)


def save_geotiff(filename, data):
    nband = 1
    nrow, ncol = data.shape
    driver = gdal.GetDriverByName("GTiff")
    dst_dataset = driver.Create(
        filename + ".tiff", ncol, nrow, nband, gdal.GDT_Float32
    )
    dst_dataset.SetGeoTransform(src.GetGeoTransform())
    dst_dataset.SetProjection(src.GetProjection())
    dst_dataset.GetRasterBand(1).WriteArray(data.astype(float))
    dst_dataset = None


b = np.sort(np.nan_to_num(deconvolved).flat)
seuil = b[np.where(np.cumsum(b) / np.sum(b) > 0.01)[0][0]]
deconvolved[deconvolved < seuil] = np.nan

# Sauvegarder l'intensité

save_geotiff("Image_deconv", deconvolved)
