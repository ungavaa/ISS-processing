#!/usr/bin/env python3

import geopandas as gpd
import illum.MultiScaleData as MSD
import rasterio
import rasterio.features
import yaml

gpd.io.file.fiona.drvsupport.supported_drivers["KML"] = "r"

with open("inputs_params.in") as f:
    p = yaml.safe_load(f)

polys = gpd.read_file("Polygones.kml")
params = polys.pop("Description").str.split("<br>", expand=True)
for i in params:
    split = params[i].str.split(":", expand=True)
    polys[split[0][0]] = split[1]
polys = polys.sort_values("ID")
polys.index = sorted(polys.index + 1)

zones = MSD.from_domain("domain.ini")
polys = polys.to_crs(zones._attrs["srs"])
for i, layer in enumerate(zones):
    min_x = zones._attrs["layers"][i]["xmin"]
    max_y = zones._attrs["layers"][i]["ymax"]
    p_size = zones._attrs["layers"][i]["pixel_size"]

    shapes = zip(polys.geometry, polys.index)
    transform = rasterio.transform.from_origin(min_x, max_y, p_size, p_size)
    rasterio.features.rasterize(
        shapes=shapes, fill=0, out=layer, transform=transform
    )
zones.save("zones")

equiv = [
    ["hlamp", "altlp"],
    ["hobs", "obsth"],
    ["dobs", "obstd"],
    ["opacity", "obstf"],
]
for key, name in equiv:
    ds = zones.copy()
    for i, layer in enumerate(ds):
        for z in polys.index:
            ds[i][layer == z] = float(polys[key][z])
    ds.save(f"Inputs/{p['exp_name']}_{name}")
