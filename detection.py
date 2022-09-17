#!/usr/bin/env python3

import os

import astropy.io.fits as pyfits
import numpy as np
import yaml
from astropy.convolution import Box2DKernel, convolve
from osgeo import gdal, osr
from scipy import stats
from scipy.ndimage import gaussian_laplace, median_filter, minimum_filter


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


with open("iss_params.in") as f:
    p = yaml.safe_load(f)

src, image_intensity = open_tiff("./processing/Vrad.tiff")


def threshold_detection(image_intensity, window):
    im = np.zeros(image_intensity.shape)
    image_nonan = np.nan_to_num(image_intensity, nan=0)
    local_emission = median_filter(
        image_nonan, window
    )  # The value of the pixels is changed by the median of the window
    # local_emission_sum = (window**2)*convolve(local_emission,Box2DKernel(width=window)) # as the convolve Box2Dkernel adds a factor (1/window^2)
    # factor=4/window
    factor = 1
    threshold = factor * local_emission
    im[image_intensity > threshold] = 1

    return im


def rings(image, num_rings):
    im = np.zeros(image.shape)
    imaslope = np.zeros(image.shape)

    for i in range(0, image.shape[0]):
        for j in range(0, image.shape[1]):
            if image[i, j] > 0:
                # We work with 1 ring (3x3)
                yf = np.arange(num_rings) * 1.0  # array([0., 1.])
                xf = np.arange(num_rings) + 1.0  # array([1., 2.])

                yf0 = np.nansum(image[i - 1 : i + 2, j - 1 : j + 2])
                yf[0] = (
                    yf0 - image[i, j]
                ) / 8.0  # mean of the first ring (3x3 minus center)

                if yf[0] < image[i, j]:
                    if num_rings == 2:
                        yf1 = np.nansum(image[i - 2 : i + 3, j - 2 : j + 3])
                        yf[1] = (
                            yf1 - yf0
                        ) / 16.0  # mean of the second ring (5x5 minus (first ring + center))
                        if yf[1] < image[i, j] and (yf[1] < yf[0]):
                            im[i, j] = 1
                            output_linreg = stats.linregress(xf, yf)
                            imaslope[i, j] = output_linreg[
                                0
                            ]  # slope of the linear regression
                            # imaslope_r[i,j] = output_linreg[2] #coeficients of the linear regression
                    else:
                        im[i, j] = 1
                        # output_linreg = stats.linregress(xf,yf)
                        # imaslope[i,j] = output_linreg[0] #slope of the linear regression
                        # imaslope_r[i,j] = output_linreg[2] #coeficients of the linear regression

    return im, imaslope


def gaussian(image_intensity, sigma, window):
    def local_min(laplacian_map):
        min_im = minimum_filter(laplacian_map, size=(window, window))
        min_x, min_y = np.where((laplacian_map - min_im) == 0)
        return min_x, min_y

    image_nonan = np.nan_to_num(image_intensity, nan=0)
    laplacian = gaussian_laplace(
        image_nonan, sigma
    )  # second derivative (change of the slope)
    min_x, min_y = local_min(laplacian)
    maxmap = np.zeros(image_intensity.shape)
    maxmap[min_x, min_y] = 1

    return maxmap


def rad_conservation(original_image_intensity, image_binary_lamps):
    # im=np.ones(original_image_intensity.shape)
    width = 5
    image_nonan = np.nan_to_num(original_image_intensity, nan=0)
    local_emission = (convolve(image_nonan, Box2DKernel(width))) * (width**2)

    image_final = local_emission * image_binary_lamps

    return image_final


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


if not os.path.isdir(p["wd"]):
    os.makedirs(p["wd"])

im_binary_threshold = threshold_detection(image_intensity, 3)
# save_geotiff(f"{p['wd']}/threshold_lamps", im_binary_threshold)
im_binary_rings, image_slope = rings(image_intensity, 1)
# save_geotiff(f"{p['wd']}/rings_lamps", im_binary_rings)
im_binary_gaus = gaussian(image_intensity, 1, 3)
# save_geotiff(f"{p['wd']}/gauss_lamps", im_binary_gaus)

im_binary_lamps = im_binary_threshold * im_binary_rings * im_binary_gaus
save_geotiff(f"{p['wd']}/binary_lamps", im_binary_lamps)

# im_conservation=rad_conservation(image_intensity,im_binary_lamps)
