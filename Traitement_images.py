import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.utils.data import download_file
from scipy import stats

image_file = fits.open('Corr_iss035e017088_ImpactVlG_GR.fits')
# image_file.info()



image_intensity =image_file[0].data
print(type(image_intensity))
print(image_intensity.shape)
image_file.close()
image_file2 = fits.open('Corr_iss035e017088Composite.fits')
image_tech = image_file2[0].data

ncol=np.size(image_tech,1)
nrow=np.size(image_tech,0)

# Élimination du bruit dans l'image
mode= stats.mode (image_intensity[~np.isnan(image_intensity)],axis=None)[0][0]
image_intensity-=mode
image_intensity [image_intensity<0] = np.nan

# Élimination des NaN isolés sur la terre
from astropy.convolution import Gaussian2DKernel, interpolate_replace_nans
kernel = Gaussian2DKernel(x_stddev=1)
image_intensity_interp = interpolate_replace_nans(image_intensity, kernel)


# Copie de l'image initiale pour manipulations ultérieures
image_intensity_interp_num = np.nan_to_num(image_intensity_interp)

# Élimination des pixels isolés dans le fleuve
	# Création de l'image binaire
image_intensity_NaN = np.isnan(image_intensity)
image_intensity[image_intensity_NaN] = 0
image_intensity [image_intensity>0] = 1
Image_binaire = image_intensity

	# Filtre Gaussien sur image binaire (avec seuil de bascule)
from scipy.ndimage import gaussian_filter
mask = gaussian_filter(image_intensity, sigma=2)
image_intensity_binaire_finale = (mask>0.83).astype(float)

# Image finale d'intensité: image binaire modifiée multipliée par image interpolée pour éliminer les NaN isolés 
image_intensity_treated = image_intensity_binaire_finale * image_intensity_interp_num

# Concordance entre les deux images (0 aux mêmes endroits)
image_tech2 = image_intensity_binaire_finale * image_tech

# Copie de l'image des technologie initiale pour manipulations ultérieures
image_tech_final = image_tech2.copy() 
plt.figure()
plt.imshow(image_tech_final, cmap='rainbow')
plt.colorbar()
plt.title('Technologies, final')

# Création image des technologies binaire
image_tech2[image_tech2>0] = 1

# Concordance entre les deux images 
image_intensity_finale = image_tech2 * image_intensity_treated
plt.figure()
plt.imshow(image_intensity_finale, cmap='rainbow')
plt.colorbar()
plt.title('Intensité lumineuse, final')

# Fonction de floutage
from astropy.convolution import convolve
def Convolution_sans_zero(image,stdev=2):
	image=image.copy()
	mask = image==0
	image[mask]=np.nan
	kernel = Gaussian2DKernel(x_stddev=stdev)
	Blur = convolve(image, kernel)
	Blur[mask]=np.nan
	return Blur

# Floutage de l'image de l'intensité
Blur_intensity=Convolution_sans_zero(image_intensity_finale,stdev=2)
plt.figure()
plt.imshow(Blur_intensity, cmap= 'rainbow')
plt.colorbar()
plt.title('Intensité flouée')

# Variables des valeurs de MSI associées aux catégories de technologies
MSI_1=0.539
MSI_2=0.446
MSI_3=0.274
MSI_4=0.043
MSI_5=0.118
MSI_6=0.017

# Création carte MSI par remplacement des valeurs de 1 à 6 par les valeurs de MSI associées
image_MSI = np.zeros_like(image_tech_final)
image_MSI[image_tech_final == 1] = MSI_1
image_MSI[image_tech_final == 2] = MSI_2
image_MSI[image_tech_final == 3] = MSI_3
image_MSI[image_tech_final == 4] = MSI_4
image_MSI[image_tech_final == 5] = MSI_5
image_MSI[image_tech_final == 6] = MSI_6
plt.figure()
plt.imshow(image_MSI, cmap='rainbow')
plt.colorbar()
plt.title('Carte MSI')

# Floutage de la carte des MSI
# stdev dépend de la résolution de pixels
Blur_MSI = Convolution_sans_zero(image_MSI,stdev=2)
plt.figure()
plt.imshow(Blur_MSI, cmap = 'rainbow')
plt.colorbar()
plt.title('MSI flouée')

# Création de la carte de l'impact final (MSI et intensité combinés)
Impact_MSI=Blur_MSI*Blur_intensity
plt.figure()
plt.imshow(Impact_MSI, cmap="rainbow")
plt.colorbar()
plt.title('Impact MSI')
# plt.show()




from osgeo import gdal
 
# final_data is a 2-D Numpy array of the same dimensions as src_data
final_data = Impact_MSI
 
# get parameters
# geotransform = src_dataset.GetGeoTransform()
# spatialreference = src_dataset.GetProjection()

nband = 1
print(ncol,nrow)
print(type(ncol))
# create dataset for output
fmt = 'GTiff'
driver = gdal.GetDriverByName(fmt)
dst_dataset = driver.Create('C:/Users/nikki/OneDrive/Bureau/tota.tiff', ncol, nrow, nband, gdal.GDT_Float32)
# dst_dataset.SetGeoTransform(geotransform)
# dst_dataset.SetProjection(spatialreference)
dst_dataset.GetRasterBand(1).WriteArray(final_data)
dst_dataset = None

# np.save('Impact_MSI',Impact_MSI)
# np.save('Blur_MSI',Blur_MSI)
# np.save('Blur_Vrad',Blur_intensity)