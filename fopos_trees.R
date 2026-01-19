library(lidR)
library(sf)

las      <- readLAS("data/USGS_LPC_NJ_MERCERCO_2009_000164.laz")
poly     <- st_read("data/riparian_east.shp")

# crs and bounding boxes are incompatible...
# check crs
st_crs(las)
st_crs(poly)
# check bounding boxes
st_bbox(las)
st_bbox(poly)

# usgs laz metadata file refers to NAVD88, which equals EPSG:8760, a compound coordinate reference system (CRS) for NJ
# set crs of las to EPSG:8760
las <- st_set_crs(las, 8760) 
# set crs of the polygon layer to match
poly <- st_transform(poly, crs = st_crs(las))

# re-check crs
st_crs(las)
st_crs(poly)
# re-check bounding boxes
st_bbox(las)
st_bbox(poly)

# now the las can be clipped to the area
las_clipped <- clip_roi(las, poly)

plot(las_clipped)
