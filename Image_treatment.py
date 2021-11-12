import numpy as np
import matplotlib.pyplot as plt
import sys
from astropy.io import fits
from astropy.utils.data import download_file
from scipy import stats
from astropy.convolution import convolve, Gaussian2DKernel, Box2DKernel
from osgeo import gdal

basename = "./ISSparis12_image/Corr_iss032e017501"

def open_tiff(filename,dtype=np.float32):
	# Load file, and access the band and get a NumPy array
	src = gdal.Open(filename, gdal.GA_Update)
	band = src.GetRasterBand(1)
	ar = band.ReadAsArray()
	return src, ar

src, image_intensity = open_tiff(basename + '_ImpactVlG_GR.tiff')
src, image_tech = open_tiff(basename + 'Composite.tiff')

src, image_intensity = open_tiff(basename + '_ImpactVlG_GR.tiff')
src, image_tech = open_tiff(basename + 'Composite.tiff')
## 1. Start of treatment : elimination of noise
## 1a. We first eliminate negative data, as these are a mistake resulting of the pre-treatment
image_intensity[image_intensity<0] = np.nan

## 1b. Statistical data
def compute_stats(image):
	mean = np.nanmean(image)
	median = np.nanmedian(image)
	standard_deviation = np.nanstd(image)
	mode = stats.mode(image[~np.isnan(image)],axis=None)[0][0]
	print(f"{mode=} {mean=} {median=} {standard_deviation=}")
	return mean, median, standard_deviation, mode

mean, median, standard_deviation, mode = compute_stats(image_intensity)

## 1c. We create a temporary map with the pixels of less value since the noise should be
##     smaller or around the same value as the less valuable pixels.
##     We then extract statistical data of these small values
small_values = image_intensity[image_intensity<standard_deviation]
mean2, median2, standard_deviation2, mode2 = compute_stats(small_values)

## 1d. Using the statistical data on hand, we estimate the value of the noise.
noise = mode + standard_deviation2
## this value can change. Modifie accordingly in order to achieve the best results (see step 1.f)
image_intensity -= noise

## 1e. Eliminating negative pixels created by our treatment of the noise in the image
image_intensity[image_intensity<0] = np.nan

## 1f. Compare your image to the ones in ReadMe for reference as well as Google maps to know which zones should emit light or not
##     We want to eliminate a certain quantity of small value pixels in aeras that shouldn't emit light
##     At the same time, we want to limit the values eliminated in aeras where light should be emitted
plt.figure()
plt.imshow(image_intensity, cmap="rainbow")
plt.colorbar()
plt.title('Intensity without noise')
plt.show()


## 2. Elimination of remaining values in dark aeras
##    If your image is already free of valued pixels in dark aeras following step 1,
##    comment this section and skip to step 3


## 2a. We define our convolution fonction
def image_to_binary(image):
	im = image.copy()
	im[im>=0] = 1
	im[np.isnan(image)] = 0
	return im

def convolution_nb_nan(image, window, keep_value):
	im = image.copy()
	threshold = keep_value / window**2
	nb_nan = convolve(image_to_binary(image), Box2DKernel(width=window))
	im[window**2 * nb_nan < keep_value] = np.nan
	return np.nan_to_num(im)

## 2b. We create a binary image where value pixels are equal to 1 and NaN pixels are equal to 0
image_intensity = convolution_nb_nan( image_intensity, window=3, keep_value=4 )

## 2c. Second convolution if necessary. If not, skip to step 3
image_intensity = convolution_nb_nan( image_intensity, window=5, keep_value=12 )

## 2d. Compare your image to the examples in the ReadMe
plt.figure()
plt.imshow(image_intensity, cmap="rainbow")
plt.colorbar()
plt.title('Intensity with clean dark aeras')
plt.show()


## 3. Concordance...

## 3a. We create a new binary map
image_binary = image_to_binary(image_intensity)

## 3b. By multiplying our technology map with the binary of intensity,
##     we eliminate data where we have none in the intensity image
image_tech = image_tech * image_binary
image_tech[image_tech==0] = np.nan


## 4. Creation of MSI map

MSI_array = [ 0.62, 0.43, 0.35, 0.08, 0.118, 0.017 ]
image_MSI = np.zeros_like(image_tech)
for i, msi in enumerate(MSI_array, 1):
	image_MSI[image_tech == i] = msi


## 5. Bluring

## 5a. Defining standard deviation
##     All measurements are in meters
pix_size = 8e-6 # pixel size
focal_distance = 400e-3 # focal distance
ISS_altitude = 400e3 # altutide of ISS
distance_dev = 25 # influence standard deviation of a single lamp
deviation = (distance_dev * focal_distance)/(pix_size*ISS_altitude) # influence standard deviation in pixel

## 5b. Defining the bluring fonction
def Convolution_without_zero(image,stdev=deviation):
	im = image.copy()
	mask = im==0
	im[mask] = np.nan
	blurred = convolve(im, Gaussian2DKernel(x_stddev=stdev))
	blurred[mask] = np.nan
	return blurred

## 5c. Bluring intensity image
Blur_intensity = Convolution_without_zero(image_intensity)

## 5d. Bluring MSI image
Blur_MSI = Convolution_without_zero(image_MSI)

## 6. Creation of impact MSI image
Impact_MSI = Blur_MSI*Blur_intensity

# Uncomment to show final images
plt.figure()
plt.imshow(Blur_MSI, cmap='rainbow')
plt.colorbar()
plt.title('Blured MSI')

plt.figure()
plt.imshow(Impact_MSI, cmap="rainbow")
plt.colorbar()
plt.title('Impact MSI')

plt.figure()
plt.imshow(Blur_intensity, cmap="rainbow")
plt.colorbar()
plt.title("Blured Intensity")
plt.show()

plt.figure()
plt.imshow(image_tech, cmap="rainbow")
plt.colorbar()
plt.title("Technology")
plt.show()

## 7. Saving results
def save_geotiff( filename, data ):
	nband = 1
	nrow, ncol = data.shape
	driver = gdal.GetDriverByName('GTiff')
	dst_dataset = driver.Create(filename+".tiff", ncol, nrow, nband, gdal.GDT_Float32)
	dst_dataset.SetGeoTransform(src.GetGeoTransform())  ##sets same geotransform as input
	dst_dataset.SetProjection(src.GetProjection())  ##sets same projection as input
	dst_dataset.GetRasterBand(1).WriteArray(data.astype(float))
	dst_dataset = None


## 7a. Saving intensity

save_geotiff('Image_Vrad',Blur_intensity)
np.save('Image_Vrad',Blur_intensity)

## 7b. Saving MSI
save_geotiff('Image_MSI',Blur_MSI)
np.save('Image_MSI',Blur_MSI)

## 7c. Saving Impact MSI
save_geotiff('Impact_MSI',Impact_MSI)
np.save('Impact_MSI',Impact_MSI)

## 7d. Saving Technology
save_geotiff('tech_image',image_tech)
np.save('tech_image',image_tech)
