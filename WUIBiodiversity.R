rm(list=ls())

#load relevant packages
require(data.table)
require(sf)
require(raster)
require(tidyr)
require(lubridate)
require(zetadiv)
require(terra)
require(fossil)
require(dplyr)

#Set random number string
set.seed(1)

# Load users working directory
source('rvar/var.R')

# Set working directory
# RVar_wd should be stored in rvar/var.R
setwd(RVar_wd)

#Set sf settings to reduce risk of errors
sf_use_s2(FALSE)

#Read in GBIF data
gbif_input <- fread(input="0064246-251120083545085.zip",sep="\t")

#Retain species-level occurrences.
gbif_input <- gbif_input[gbif_input$species!="",]

#Filter occurrences data
gbif_input <- gbif_input[gbif_input$class=="Insecta",]

#Make a spatial points object
gbif_spatial <- st_as_sf(gbif_input,coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

#Wildfire Hazard Potential data from https://www.fs.usda.gov/rds/archive/catalog/RDS-2015-0047-4
WHP_input <- rast("whp2023_cnt_conus.tif")
WHP_input <- terra::project(WHP_input, "EPSG:4326")

#Extract WHP values
WHP_points <- terra::extract(WHP_input,gbif_spatial)

#Add in WHP values
gbif_spatial$WHP <- WHP_points$Band_1


#Get Fire Hazard Severity Zones map data from https://www.lab.data.ca.gov/dataset/fire-hazard-severity-zones-in-sra-effective-april-1-2024-with-lra-recommended-2007-2011?utm_source=chatgpt.com
FHSZ_input <- st_read("FHSZ_SRA_LRA_Combined.shp")
#Reproject
FHSZ_input <- st_transform(FHSZ_input, 4326)

#Add in FHSZ values to the sample locations.
FHSZ_spatial <- st_intersection(gbif_spatial,FHSZ_input[,"FHSZ"])
tmp <- st_coordinates(FHSZ_spatial)
colnames(tmp) <- c("decimalLongitude","decimalLatitude")
FHSZ_spatial <- st_drop_geometry(FHSZ_spatial)
FHSZ_spatial <- cbind(FHSZ_spatial,tmp)
#Remove duplicates
FHSZ_spatial <- data.table(FHSZ_spatial[!duplicated(FHSZ_spatial),])

#Set unique sample IDs
FHSZ_spatial[, sampleid := .GRP, by = .(eventDate,decimalLongitude,decimalLatitude)]
#Remove sites without FHSZ values.
FHSZ_spatial <- FHSZ_spatial[!is.na(FHSZ_spatial$FHSZ),]
#Remove undersampled sites
n_threshold <- 20
FHSZ_filtered <- FHSZ_spatial %>%
  group_by(sampleid) %>%
  filter(n() >= n_threshold) %>%
  ungroup()
#Set FHSZ values as factors
FHSZ_filtered$FHSZ <- as.factor(FHSZ_filtered$FHSZ)

#Set unique sample IDs
FHSZ_spatial[, sampleid := .GRP, by = .(eventDate,decimalLongitude,decimalLatitude)]
#Remove sites without FHSZ values.
FHSZ_spatial <- FHSZ_spatial[!is.na(FHSZ_spatial$FHSZ),]
#Remove undersampled sites
n_threshold <- 20
FHSZ_filtered <- FHSZ_spatial %>%
  group_by(sampleid) %>%
  filter(n() >= n_threshold) %>% ungroup()
#Set FHSZ values as factors
FHSZ_filtered$FHSZ <- as.factor(FHSZ_filtered$FHSZ)

#Create a presence/absence data table.
taxa_count <- FHSZ_filtered[,c("sampleid","taxonKey")]
taxa_count <- taxa_count[!duplicated(taxa_count),]
taxa_count <- taxa_count %>%
  dplyr::group_by(sampleid, taxonKey) %>%
  dplyr::mutate(taxonCount = n()) %>%
  ungroup()
taxa_count <- taxa_count[!duplicated(taxa_count),]
gbif_pa <- taxa_count %>%
  pivot_wider(
    names_from = taxonKey,
    values_from = taxonCount,
    values_fill = 0   # fill missing combinations with 0 (or NA if you prefer)
  )
gbif_pa <- data.frame(gbif_pa)
rownames(gbif_pa) <- gbif_pa$sampleid
gbif_pa$sampleid <- NULL


#Set environmental data set
env_sampled <- FHSZ_filtered[,c("sampleid","FHSZ","WHP")]
env_sampled$WHP <- as.factor(env_sampled$WHP)
env_sampled <- env_sampled[!duplicated(env_sampled),]
env_sampled <- data.frame(env_sampled[,c("FHSZ","WHP")])
#Set location data set
location_sampled <- FHSZ_filtered[,c("sampleid","decimalLongitude","decimalLatitude")]
location_sampled <- location_sampled[!duplicated(location_sampled),]
location_sampled$sampleid <- NULL
#Calculate how much zeta diversity varies with WUI and distance for common (High zeta order) and rare species (Low zeta order).
zeta_low <- 2
zeta_high <- 10
#Zeta.varpart returns a data frame with one column containing the variation explained by each component:
#a (the variation explained by distance alone)
#b (the variation explained by either distance or the environment)
#c (the variation explained by the environment alone)
#d (the unexplained variation).
Zeta.varpart(Zeta.msgdm(data.spec=gbif_pa,data.env=env_sampled,xy=location_sampled,order=zeta_low))
Zeta.varpart(Zeta.msgdm(data.spec=gbif_pa,data.env=env_sampled,xy=location_sampled,order=zeta_high))


#Estimate richness against WUI.  The WUI categories correspond to:
#0  non-WUI / background / no classification
#1	Forest/Shrubland/Wetland-dominated Intermix WUI — areas where buildings and wildland vegetation are intermixed and the wildland component is forest/shrub/wetland. 
#2	Forest/Shrubland/Wetland-dominated Interface WUI — buildings are near large patches of forest/shrub/wetland vegetation (but not mixed within). 
#3	Grassland-dominated Intermix WUI — buildings intermingle with grasslands. 
#4	Grassland-dominated Interface WUI — buildings near large grassland areas. 
#5	Non-WUI: Forest/Shrubland/Wetland-dominated — wildland vegetation area not classified as WUI. 
#6	Non-WUI: Grassland-dominated — grassland not in WUI. 
#7	Non-WUI: Urban — built/urban land not in a WUI context. 
#8	Non-WUI: Other — other land types (e.g., water, bare land) outside WUI definitions.
