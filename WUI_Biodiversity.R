rm(list=ls())
# Load relevant packages

require(terra)
require(raster)
require(sf)
require(tidyverse)
require(data.table)
require(fossil)
require(zetadiv)

# Load users working directory
source('rvar/var.R')

# Set random number string
set.seed(1)

# Set working directory
# RVar_wd should be stored in rvar/var.R
setwd(RVar_wd)

#Read in GBIF Data
gbif_input <- fread(input='.fullDatasets/0064246-251120083545085.csv',sep='\t')

#Filter data to phylum Arthropoda (Insects)
gbif_input <- gbif_input[gbif_input$class=='Insecta',]

#From the filtered gbif data, keep each entry that has a specified species name.
gbif_input <- gbif_input[gbif_input$species!='',]

#Get spatial co-ords to from the gbif records, to be able to filter WUI data
gbif_spatial <- st_as_sf(gbif_input, coords=c('decimalLongitude','decimalLatitude'),crs=4326)

#Study Area
if(file.exists('data/Study_counties.shp')){
  Study_counties <- st_read('data/Study_counties.shp')
} else {
  #Create list of Study counties
  Study_counties <- c('Los Angeles','Orange','San Bernardino','San Diego','Santa Barbara','Ventura')
  
  #CA County boundaries: https://data.ca.gov/dataset/ca-geographic-boundaries
  Study_counties <- read_sf('.fullDatasets/CA_Counties.shp') %>% dplyr::filter(NAME %in% Study_counties)
  
  #Only need the name as an identifier
  Study_counties <- Study_counties %>% select(NAME)
  
  
  #Transform to 4326 co-ord system
  Study_counties <- st_transform(Study_counties,st_crs(gbif_spatial))
  
  st_write(Study_counties,'data/Study_counties.shp',append = FALSE)
}
  
#Plotting gbif data with CA Counties
#Plot GBIF Data & CA Counties
ggplot() +
  geom_sf(data=Study_counties, fill = 'lightgrey',color='black') +
  geom_sf(data=gbif_spatial, size=1,alpha=0.7, color='orange')
  theme_minimal()
  
#Get number of records per county
Species_by_county <- st_join(gbif_spatial,Study_counties %>% select('NAME'), left=FALSE) %>% pull(NAME) %>% table
# LA: 85, O: 35, SBer: 5387, SD: 1933, SBar:5886, V:42
# TODO: Does the geographic spread skew our analysis?
# Look at median WUI for the all the samples in each county?

#WHP
#WHP Data: https://www.fs.usda.gov/rds/archive/catalog/RDS-2015-0047-4
if(file.exists('data/WHP_SA.tif')){
  WHP_Data <- rast('data/WHP_SA.tif')
} else {
  WHP_Data <- rast('.fullDatasets/whp2023_GeoTIF/whp2023_cnt_conus.tif')
  #Convert to 4326 (Data is 5070)
  WHP_Data <- terra::project(WHP_Data,'EPSG:4326')
  
  #Crop data to Study Area
  WHP_Data <- WHP_Data %>% crop(Study_counties) %>% mask(Study_counties)
  
  #WHP uses Max 32-bit INT to denote NA. Checking layer in QGIS, actual max value is ~144,000
  #SOURCE: https://imagery.geoplatform.gov/iipp/rest/services/Fire_Aviation/USFS_EDW_RMRS_WildfireHazardPotentialContinuous/ImageServer
  WHP_Data[WHP_Data > 145000] <- NA
  
  #Write raster for later loading/github
  raster::writeRaster(WHP_Data,'data/WHP_SA.tif',overwrite=TRUE)

  # LOOK INTO HOW WHP IS IMPACTED BY ORNAMENTAL PLANTS
}

#FHSZ
#FHSZ Data: https://www.lab.data.ca.gov/dataset/fire-hazard-severity-zones-in-sra-effective-april-1-2024-with-lra-recommended-2007-2011
if(file.exists('data/FHSZ_Data.shp')){
  FHSZ_Data <- st_read('data/FHSZ_Data.shp')
} else {
  FHSZ_Data <- st_read('.fullDatasets/FHSZ_SRA_LRA_Combined.shp')
  #Convert to 4326 (Data is 3310)
  FHSZ_Data <- st_transform(FHSZ_Data,st_crs(gbif_spatial))
  
  FHSZ_Data <- st_make_valid(FHSZ_Data)
  
  FHSZ_Data <- st_intersection(FHSZ_Data,Study_counties)
  
  FHSZ_Data <- FHSZ_Data %>% select(FHSZ)
  
  st_write(FHSZ_Data,'data/FHSZ_Data.shp',append=FALSE)
}


#WUI
if(file.exists('data/WUI_SA.tif')){
  WUI_Data <- rast('data/WUI_SA.tif')
} else {
  #Download Wildland Urban Interface raster of North America from: https://geoserver.silvis.forest.wisc.edu/geodata/globalwui/NA.zip
  #Gather a list of all directories in the WUI folder that begin with X: (Coordinate System )
  WUI_tiles <- paste(list.dirs('.fullDatasets/NA')[grepl('.fullDatasets/NA/X',list.dirs('.fullDatasets/NA'))],'/WUI.tif',sep='')
  
  #Get CRS of the rasters
  WUI_tile <- rast(WUI_tiles[1])
  WUI_crs <- crs(WUI_tile)
  
  #Project Study Counties data into WUI:Faster
  Study_counties_WUI_projection <- st_transform(Study_counties,WUI_crs)
  
  #Build an index of all the tiles, including a polygon of the bounds of each
  WUI_tiles_index <- map_dfr(WUI_tiles,function(file){
    
    #Get the header only of a raster file
    WUI_tile_header <- rast(file)
    #Get the extent of that header
    WUI_tile_extent <- ext(WUI_tile_header)
    
    st_as_sf(as.polygons(WUI_tile_extent,crs=WUI_crs)) %>% mutate(filename=file)
  })
  
  #Find the tiles that overlap with Study counties
  WUI_tiles <- WUI_tiles_index %>% st_filter(Study_counties_WUI_projection) %>% pull(filename) %>% unique()
  
  
  #Combine into single WUI
  WUI_Data <- vrt(WUI_tiles)
  
  #Crop & Mask to counties
  WUI_Data <- crop(WUI_Data,Study_counties_WUI_projection)
  WUI_Data <- mask(WUI_Data,Study_counties_WUI_projection)
  
  #Reproject to 4326 ~1h 30mins to complete
  WUI_Data <- project(WUI_Data,'EPSG:4326')
  
  #Write out the raster
  raster::writeRaster(WUI_Data,'data/WUI_SA.tif',overwrite=TRUE)
  
  #Drop unnecessary vars
  rm(WUI_tile,WUI_crs,Study_counties_WUI_projection,WUI_tiles_index)
}

#Check for Correlations between FHSZ & WUI, WHP & WUI


#Combining FIRE Data: Do a FAMD to create an equivalent of PCA -> Single 'Fire Risk' statistic -> RF Model? -> Cluster

#0  non-WUI / background / no classification
#1	Forest/Shrubland/Wetland-dominated Intermix WUI — areas where buildings and wildland vegetation are intermixed and the wildland component is forest/shrub/wetland. 
#2	Forest/Shrubland/Wetland-dominated Interface WUI — buildings are near large patches of forest/shrub/wetland vegetation (but not mixed within). 
#3	Grassland-dominated Intermix WUI — buildings intermingle with grasslands. 
#4	Grassland-dominated Interface WUI — buildings near large grassland areas. 
#5	Non-WUI: Forest/Shrubland/Wetland-dominated — wildland vegetation area not classified as WUI. 
#6	Non-WUI: Grassland-dominated — grassland not in WUI. 
#7	Non-WUI: Urban — built/urban land not in a WUI context. 
#8	Non-WUI: Other — other land types (e.g., water, bare land) outside WUI definitions.