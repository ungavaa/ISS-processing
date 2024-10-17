import numpy as np
import pandas as pd
import rasterio as rio
import skimage

rst = rio.open("Vrad_georef_modif.tif")
img = np.nan_to_num(rst.read()[0])
#img += 10
psf = np.ones((5, 5)) / 25
print('1')
#deconv = skimage.restoration.richardson_lucy(
#    img, psf, 30, clip=False, filter_epsilon=0.0001
#)
peaks = skimage.feature.peak_local_max(img)
print('2')
mask = np.zeros_like(img, dtype=int)
mask[tuple(peaks.T)] = np.arange(len(peaks))+1
print('3')
ws = skimage.segmentation.watershed(-img, mask)
print('4')
#img -= 10
df = pd.DataFrame(
    skimage.measure.regionprops_table(
        ws, img, properties=["centroid_weighted", "image_intensity"]
    )
)
print('5')
df["image_intensity-sum"] = df["image_intensity"].map(np.sum)
print('6')
df["longitude"], df["latitude"] = rio.transform.xy(
    rst.transform, df["centroid_weighted-0"], df["centroid_weighted-1"]
)
print('7')


df.to_csv(
    "Vrad_georef_modif.csv",
    index=False,
    columns=[
        "longitude",
        "latitude",
        "image_intensity-sum",
    ],
    header=["lon", "lat", "flux"],
)
print('8')