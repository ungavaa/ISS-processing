#!/usr/bin/python3
#
# ====================
import getopt
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from datetime import date, timedelta
from astropy.convolution import Box2DKernel, Gaussian2DKernel, convolve
from astropy.stats import SigmaClip, sigma_clipped_stats
from astropy.visualization import SqrtStretch, simple_norm
from astropy.visualization.mpl_normalize import ImageNormalize
from numpy import inf
from photutils.aperture import CircularAperture
from photutils.background import Background2D, MedianBackground
from photutils.detection import DAOStarFinder
from scipy import ndimage, stats
from scipy.ndimage import gaussian_filter
from scipy.spatial.distance import cdist
from scipy.optimize import curve_fit
from osgeo import gdal, osr
import statistics as st

# reading command line parameters
def input(argv):
    Ifile = "undefined"
    try:
        opts, args = getopt.getopt(
            argv,
            "h:i:d:b:",
            [
                "help=",
                "ifile=",
            ],
        )
    except getopt.GetoptError:
        print("ProcessNighMon.py -i <Ifile>")
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print("ProcessNighMon.py -i <Ifile>")
            sys.exit()
        elif opt in ("-i", "--ifile"):
            Ifile = arg
    print("Image file is :", Ifile)
    return Ifile


def open_tiff(filename, dtype=np.float32):
    if filename.endswith("fits"):
        filename = fits2tiff(filename)

    # Load file, and access the band and get a NumPy array
    src = gdal.Open(filename, gdal.GA_Update)
    band = src.GetRasterBand(1)
    ar = band.ReadAsArray()
    return src, ar


def save_geotiff(filename, data):
    nband = 1
    nrow, ncol = data.shape
    driver = gdal.GetDriverByName("GTiff")
    dst_dataset = driver.Create(filename + ".tiff", ncol, nrow, nband, gdal.GDT_Float32)
    # sets same geotransform as input
    dst_dataset.SetGeoTransform(src.GetGeoTransform())
    # sets same projection as input
    dst_dataset.SetProjection(src.GetProjection())
    dst_dataset.GetRasterBand(1).WriteArray(data.astype(float))
    dst_dataset = None


# ================================================
# MAIN
# default Parameters
FWHM = 13
IntegHalfSize = 1
minflux = 0.01  # W/sr/m^2 Photopic
# load command line parameters
Ifile = input(sys.argv[1:])
outname = Ifile + ".csv"
# create the data file if it do not exists
if os.path.exists(outname) == False:
    o = open(outname, "w")
    first_line = "Xpos,Ypos,Flux \n"
    second_line = "(pixel),(pixel),(W/sr/m^2) \n"
    o.write(first_line)
    o.write(second_line)
    o.close()
src, imag = open_tiff(f"{Ifile}")
ny = np.shape(imag)[0]
nx = np.shape(imag)[1]
imag[np.isnan(imag)] = 0
# Search for sources
mean, median, std = sigma_clipped_stats(imag, sigma=3.0)
daofind = DAOStarFinder(fwhm=1.5, threshold=5 * std)
sources = daofind(imag)

# print(sources)
positions = np.transpose((sources["xcentroid"], sources["ycentroid"]))
# positions = (np.rint(positions)).astype(int)
print("Detected points : ", np.shape(positions)[0])
Flux = np.zeros(np.shape(positions)[0])
Back = np.zeros(np.shape(positions)[0])
for nd in range(np.shape(positions)[0]):
    xsa = round(positions[nd, 0])
    ysa = round(positions[nd, 1])
    Flux[nd] = np.sum(
        imag[
            ysa - IntegHalfSize : ysa + 2 * IntegHalfSize,
            xsa - IntegHalfSize : xsa + 2 * IntegHalfSize,
        ]
    )
sources = sources[Flux > minflux]
positions = positions[Flux > minflux]
norm = ImageNormalize(stretch=SqrtStretch())
plt.figure()
apertures = CircularAperture(positions, r=IntegHalfSize + 0.5)
plt.imshow(imag, cmap="magma", origin="lower", norm=norm, interpolation="nearest")
plt.colorbar()
apertures.plot(color="white", lw=1.0, alpha=0.3)
Flux = Flux[Flux > minflux]
# TODO : transform x, y into lat,lon if necessary and save it in output csv
imagout = np.zeros([ny, nx])
imagout[:] = np.nan
# write the sources to the output csv
print("Number of light points : ", np.shape(Flux)[0])
for no in range(np.shape(Flux)[0]):
    o = open(outname, "a")
    outputline = (
        str(positions[no, 0])
        + ","
        + str(positions[no, 1])
        + ","
        + str("{:6.3f}".format(Flux[no]))
        + "\n"
    )
    # print(no, outputline)
    o.write(outputline)
    o.close()
    imagout[round(positions[no, 1]), round(positions[no, 0])] = Flux[no]
# plt.show()
save_geotiff("Output", imagout)
