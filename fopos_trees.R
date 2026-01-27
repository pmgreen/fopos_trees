# analyze tree heights and density at Mountain Lakes Preserve in Princeton, NJ
# main reference: https://lab.jonesctr.org/ -- straight up copied and adjusted
# Ubuntu people (lidRviewer dependencies):
#   sudo apt-get install libsdl2-dev freeglut3-dev
#   install.packages("lidRviewer", repos = c("https://r-lidar.r-universe.dev", "https://cloud.r-project.org"))   

library(dplyr)
library(lidR)
library(lidRviewer)
library(future)
library(ggplot2)
library(RCSF)
library(sf)
library(terra)

# select and download LiDAR from USGS LiDAR explorer
las <- readLAS("data/USGS_LPC_NJ_MERCERCO_2009_000164.laz")
colnames(las@data) # just examine the field names

# read in a polygon representing the area of interest (e.g. shapefile)
riparian_east <- st_read("data/riparian_east.shp")

# Classify all isolated returns
# ivf is one of the two LiDR noise segmentation algorithms (sor is the other)
#las = classify_noise(las, ivf(res = 3, n = 10))
# The classify_noise function sets Classification = 18 (LASNOISE = 18); filter it out
#las = filter_poi(las, Classification != LASNOISE)

plot(las) # plot it
view(las) # or open it lidRviewer

# Clip out the riparian east area, but first ...
# crs and bounding boxes of LiDAR and Riparian East polygon are incompatible
# check crs
st_crs(las)
st_crs(riparian_east)
# check the bounding boxes
st_bbox(las)
st_bbox(riparian_east)

# the USGS laz metadata file refers to 'NAVD88' which equals EPSG:8760, a compound crs for NJ so ...
# set crs of las to EPSG:8760
las <- st_set_crs(las, 8760) 
# now set crs of the polygon layer to match
riparian_east <- st_transform(riparian_east, crs = st_crs(las))

# re-check crs
st_crs(las)
st_crs(riparian_east)
# re-check bounding boxes to verify that they are compatible
st_bbox(las)
st_bbox(riparian_east)

# now the las can be clipped to the area of interest
riparian_east <- clip_roi(las, riparian_east)

plot(riparian_east) # plot it
view(riparian_east) # or open in lidRviewer

# now we can analyze the riparian east project area
# derive a canopy height model (chm)
riparian_east = classify_ground(riparian_east, algorithm = csf())
dtm = rasterize_terrain(riparian_east, res=1, algorithm=tin())
dsm = rasterize_canopy(riparian_east, algorithm = pitfree())
chm = dsm - dtm

par(mfrow = c(1,3))
plot(dsm)
plot(dtm)
plot(chm)

# estimate tree locations
treetops = locate_trees(chm, lmf(ws = 6, hmin=10))
# estimate areas without trees (< 10 feet)
nontrees = locate_trees(chm, lmf(ws = 6, hmin=0))
nontrees <- nontrees[nontrees$Z <= 10, ]

# just see how many of each category were found
nrow(treetops)
nrow(nontrees)
# ensure dataframes are identical
colnames(treetops)
colnames(nontrees)
# take a look at the heights
treetops$Z
nontrees$Z

# estimate crowns
plot(treetops$geometry, add = FALSE, pch = 16, cex = 0.2)
crown_delineation_algorithm = dalponte2016(chm*1, treetops, th_tree = 2)
crown_raster = crown_delineation_algorithm()

# convert crown raster to polygons
crowns = as.polygons(crown_raster)

# plot and make sure there are plenty of colors to tell the crowns apart
par(mfrow = c(1,2))
plot(crowns, col=pastel.colors(8000))
plot(chm); plot(crowns, border=grey(0.5), add = TRUE)

# convert treetops to vector object and write to disk
# writeVector(vect(treetops), 'data/treetops.shp')
# writeVector(crowns, 'data/crowns.shp')

# Find the number of crowns
n_trees = nrow(crowns)
n_trees

# Calculate the area of the scene. [Product of the resoultion (1x1)] * [product of the dimensions (1000 x 1000)]  
scene_area_m2 = prod(res(chm)) * prod(dim(chm)) 
scene_area_ha = scene_area_ha = scene_area_m2 / 10000

# Calculate tree density
print(n_trees/scene_area_ha) # There are ~68 trees per ha in the scene

# get area and radius of delineated crowns
crowns_sf = st_as_sf(crowns) # convert to sf object
crowns_sf$crown_area = st_area(crowns_sf)
crowns_sf$crown_radius = sqrt(crowns_sf$crown_area/pi)

# Make histograms of tree heights and crown area
par(mfrow = c(1,2), mar=c(4,4,2,2)) # make a 1x3 panel plot
## hist(crowns_sf$crown_radius, main='Histogram of crown radius', xlab='Crown radius (m)')
hist(treetops$Z, main = 'Histogram of tree height', xlab='Tree height (f)')
ggplot(treetops, aes(x = Z)) + geom_histogram(color = "#333", fill = "darkgreen", bins = 30, alpha = 0.5) + labs(x = 'heights') +
  theme_minimal()

# Make a pie chart of forested vs. non-forested area
slices <- c(nrow(treetops), nrow(nontrees))
labels <- c("forested","non-forested")
pie(slices, labels = labels, main = "forested (>10ft) vs. non-forested (<10ft)",col = c("#669900", "#F2F2F2"),border = "#999",cex = 1)

# plot crowns to heights
plot(crowns_sf$crown_radius ~ treetops$Z, xlab='Tree height (f)',
     ylab='Crown radius (m)', pch=16, col='#00000033')
