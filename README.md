# ISS-processing
Image processing based on SAVESTARS SL preprocessed ISS images
This work was made by Martin Aubé, Julien-Pierre Houle, Nikki Veilleux, Justine Desmarais, Emie Bordeleau and Alexandre Simoneau.

Outputs of the software are the Melatonin Suppression Index (MSI), the Visual radiance and the MSI impact 
defined as the product of the MSI with the visual radiance. The software can also produce inputs for the Illumina model.

### Obtain a discrete lamps inventory from ISS images
In order to use this sowftware to create discrete lights inventory, you need to first install the Illumina model. This can be done by following the installations step in the Wiki of Martin Aubé:https://lx02.cegepsherbrooke.qc.ca/~aubema/index.php/Prof/IlluminaGuide2021.
Make sure to execute: ```illum inputs``` to create the necessary files. Make sure to modify the parameters at the beginning of each execution script.

1. Convert the georeferencing ISS images (raster) to XYZ files  with the following command:
``` 
gdal_translate -of XYZ PATH_RASTER.tiff PATH_XYZ.csv
```
2. Execute extract-output-data.py to save XYZ files to Numpy arrays contening data on the domain, intensity and technology of the lamps.
3. Execute make_discrete_inv.py to create a discrete inventory of the ligths sources that can be used a an input to the Illumina model.
