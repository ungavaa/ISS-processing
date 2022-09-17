#!/usr/bin/env python3

import os
from glob import glob

import astropy.convolution as apconv
import astropy.io.fits as pyfits
import numpy as np
import yaml
from osgeo import gdal, osr


def fits2tiff(filename):
    hdulist = pyfits.open(filename)
    other = np.asarray(hdulist[0].data)
    earth = other.copy()
    earth = np.float32(earth)
    nrows, ncols = earth.shape[0], earth.shape[1]
    filename = filename[:-5] + ".tiff"
    dst_ds = gdal.GetDriverByName("GTiff").Create(
        filename, ncols, nrows, 1, gdal.GDT_Float32
    )
    sr = osr.SpatialReference()
    sr.ImportFromEPSG(4326)
    dst_ds.SetProjection(sr.ExportToWkt())
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
    return filename


def open_tiff(filename, dtype=np.float32):
    if filename.endswith("fits"):
        filename = fits2tiff(filename)

    # Load file, and access the band and get a NumPy array
    src = gdal.Open(filename, gdal.GA_Update)
    band = src.GetRasterBand(1)
    ar = band.ReadAsArray()
    return src, ar


def convolve(im, window):
    return apconv.convolve(im, apconv.Box2DKernel(width=window))


def convolve_freq(im, window):
    vals = np.unique(im)
    freqs = [convolve(im == i, window) for i in vals]
    conv = vals[np.argmax(freqs, axis=0)]
    conv[~np.any(freqs, axis=0)] = np.nan
    return conv


def fill_voids(im, window, n_keep, threshold=3):
    frac_valid = convolve(im > 0, window)
    mask = (frac_valid > n_keep / window**2) & (np.nan_to_num(im) == 0)
    im_grow = convolve(im, window)
    mask &= im_grow > np.nanmean(im) * threshold
    return np.where(mask, im_grow, im), mask


def remove_isolated(im, window, n_keep):
    frac_valid = convolve(im > 0, window)
    mask = (frac_valid < n_keep / window**2) & (np.nan_to_num(im) != 0)
    return np.where(mask, np.nan, im), mask


with open("iss_params.in") as f:
    p = yaml.safe_load(f)

exp_folds = glob("iss???e*/")
if not exp_folds:
    print("ERROR: Data folder not found. Please make sure it was not renamed.")
    exit()
if len(exp_folds) > 1:
    print(f"ERROR: Too many data folders found. {exp_folds}")
    exit()

exp_fold = exp_folds[0]
basename, ts = exp_fold.split("_")
print(f"Data `{basename}` ({ts[:4]}-{ts[4:6]}-{ts[6:8]}) found.")

# Vrad = u.load_fits(f"{exp_fold}/popular/Corr_{basename}_ImpactVlG_GR.fits")[0]
# tech = u.load_fits(f"{exp_fold}/popular/Corr_{basename}CompositeW.fits")[0]

# open intensity and technology (spectral class) images
src, Vrad = open_tiff(
    f"{exp_fold}/popular/Corr_{basename}_ImpactVlG_GR.tiff"  # fits
)

src, tech = open_tiff(f"{exp_fold}/popular/Corr_{basename}CompositeW.fits")

Vrad -= np.nanmax(Vrad) * p["threshold"]
Vrad[Vrad < 0] = np.nan

n_changed = 1
changed = np.zeros_like(Vrad, dtype=np.bool8)
while n_changed:
    Vrad, mask = fill_voids(Vrad, window=3, n_keep=3, threshold=3)
    n_changed = np.sum(mask)
    changed |= mask
    print("Void pixel filled:", n_changed)
tech[changed] = convolve_freq(tech, 3)[changed]

n_changed = 1
changed = np.zeros_like(Vrad, dtype=np.bool8)
while n_changed:
    Vrad, mask = remove_isolated(Vrad, window=5, n_keep=3)
    n_changed = np.sum(mask)
    changed |= mask
    print("Isolated pixel removed:", n_changed)
tech[changed] = np.nan

if not os.path.isdir(p["wd"]):
    os.makedirs(p["wd"])

# rst = rio.open(f"{exp_fold}/quality/{basename}RGBdisto2_rect.tiff")


def save_geotiff(filename, data):
    nband = 1
    nrow, ncol = data.shape
    driver = gdal.GetDriverByName("GTiff")
    dst_dataset = driver.Create(
        filename + ".tiff", ncol, nrow, nband, gdal.GDT_Float32
    )
    # sets same geotransform as input
    dst_dataset.SetGeoTransform(src.GetGeoTransform())
    # sets same projection as input
    dst_dataset.SetProjection(src.GetProjection())
    dst_dataset.GetRasterBand(1).WriteArray(data.astype(float))
    dst_dataset = None


save_geotiff(f"{p['wd']}/Vrad", Vrad)
save_geotiff(f"{p['wd']}/tech", tech)
# save_geotiff(f"{p['wd']}/domain", domain)
