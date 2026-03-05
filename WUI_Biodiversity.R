rm(list=ls())
# Load relevant packages

require(terra)
require(raster)
require(sf)
require(tidyverse)
require(data.table)
require(iNEXT)
require(zetadiv)

# Load users working directory
source('rvar/var.R')

# Set random number string
seed_n <- 1
set.seed(seed_n)

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
  
#Get number of records per county
Species_by_county <- st_join(gbif_spatial,Study_counties %>% select('NAME'), left=FALSE) %>% pull(NAME) %>% table
print(Species_by_county)
# LA: 85, O: 35, SBer: 5387, SD: 1933, SBar:5886, V:42

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

#WUI Data
if(file.exists('data/WUI_SA_1.tif')){
  #Get a list of the split WUI files
  WUI_list <- list.files(path='data/',pattern='WUI_SA_.*\\.tif$',full.names=TRUE)
  #Load the files into a list of spat rasters
  WUI_Data <- vrt(WUI_list,'data/WUI_Full.vrt',overwrite=TRUE)
} else {
  #Download Wildland Urban Interface raster of North America from: https://geoserver.silvis.forest.wisc.edu/geodata/globalwui/NA.zip
  #Gather a list of all directories in the WUI folder that begin with X: (Coordinate System )
  WUI_tiles <- paste(list.dirs('.fullDatasets/NA')[grepl('.fullDatasets/NA/X',list.dirs('.fullDatasets/NA'))],'/WUI.tif',sep='')
  
  #Get CRS of the rasters
  WUI_tile <- rast(WUI_tiles[1])
  WUI_crs <- crs(WUI_tile)
  
  #Project Study Counties data into WUI:Faster
  Study_counties_WUI_projection <- st_transform(Study_counties,WUI_crs)
  
  #set i
  i <- 1
  
  #Build an index of all the tiles, including a polygon of the bounds of each
  WUI_tiles_index <- map_dfr(WUI_tiles,function(file){
    
  #Print a progress update to soothe the soul
  if (i %% 50 == 0) {
    cat(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] Processing file ", i, " of ", length(WUI_tiles), "\n"))
  }
    
  #Get the header only of a raster file
  WUI_tile_header <- rast(file)
  #Get the extent of that header
  WUI_tile_extent <- ext(WUI_tile_header)
  
  #increment i
  i <<- i+1
  
  #Convert the extent into a simple polygon & attaches a tag based on the filename
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
  
  ##Split raster for GIThub
  #Select number of splits -> number of rows per split
  WUI_n_splits <- 5
  WUI_rows_per_chunk <- ceiling(nrow(WUI_Data) / WUI_n_splits)
  #Loop for each split
  for(i in 1:WUI_n_splits){
    #Define the row range for this particular chunk
    WUI_start_row <- ((i-1)*WUI_rows_per_chunk) + 1
    WUI_end_row <- min(i*WUI_rows_per_chunk,nrow(WUI_Data))
    
    #Crop the raster to those rows
    WUI_chunk <- WUI_Data[WUI_start_row:WUI_end_row,, drop=FALSE]
    
    #Write out each raster
    raster::writeRaster(WUI_chunk,paste0('data/WUI_SA_',i,'.tif'),overwrite=TRUE)
  }
  
  #Drop unnecessary vars
  rm(WUI_tiles,WUI_crs,Study_counties_WUI_projection,WUI_tiles_index,i,WUI_tile_header,WUI_tile_extent,WUI_n_splits,WUI_rows_per_chunk,WUI_start_row,WUI_end_row,WUI_chunk)
}

# Extract WHP Points
WHP_points <- terra::extract(WHP_Data,gbif_spatial)
gbif_spatial$WHP <- WHP_points$Band_1
rm(WHP_points)

# Add FHSZ Data
gbif_spatial <- st_join(gbif_spatial,FHSZ_Data, join=st_intersects)

# Add WUI Data
WUI_points <- terra::extract(WUI_Data,gbif_spatial)
gbif_spatial$WUI <- WUI_points$WUI_Full

#########
# Biodiversity Analysis: 'undersampled' thresholds of 20, 50 & 100
#########

diversity_n_thresholds <- c(20,50,100)

# Research Question: how much zeta diversity varies with WUI and distance for common (High zeta order) and rare species (Low zeta order).

for (n_threshold in diversity_n_thresholds) {
  # This errors out on windows PCs: TryCatch to handle the error gracefully
  
  
  # Group by a unique sample ID based on time & location
  # Remove undersampled Sites
  gbif_grouped <- gbif_spatial %>% 
    group_by(eventDate,as.character(geometry)) %>% 
    mutate(sampleid = cur_group_id()) %>% 
    filter(n() >= n_threshold) %>% ungroup()
  
  # Set FHSZ as a factor
  gbif_grouped$FHSZ <- as.factor(gbif_grouped$FHSZ)
  
  # Create presence/absence data
  gbif_taxa_count <- gbif_grouped[,c('sampleid','taxonKey')]
  #3,350 unique sample+species+location -> only presence noted, not number of times
  gbif_taxa_count <- gbif_taxa_count[!duplicated(gbif_taxa_count),]
  #Set to 1 i.e. species present in this sample.
  gbif_taxa_count <- gbif_taxa_count %>%
    dplyr::group_by(sampleid,taxonKey) %>%
    dplyr::mutate(taxonCount = n()) %>%
    ungroup()
  #Take the presence data, and add absence for every sampleid
  gbif_pres_abs <- gbif_taxa_count %>%
    pivot_wider(
      names_from = taxonKey,
      values_from = taxonCount,
      values_fill = 0
    )
  
  print(paste0("For n:",n_threshold,". ",ncol(gbif_pres_abs) - 2," species present"))
  
  #Create dataframe
  gbif_pres_abs <- as.data.frame(gbif_pres_abs)
  
  #Get values for WHP, FHSZ, and WUI for each sample
  gbif_data_lookup <- gbif_grouped %>%
    select(sampleid,WUI,WHP,FHSZ) %>%
    distinct() %>%
    st_drop_geometry()
  
  #Join the data to each sampleid in gbif_pres_abs
  gbif_pa_data <- gbif_pres_abs %>% left_join(gbif_data_lookup, by='sampleid')
  gbif_pa_data$WUI  <- as.factor(gbif_pa_data$WUI)
  gbif_pa_data$FHSZ <- as.factor(gbif_pa_data$FHSZ)
  rownames(gbif_pa_data) <- gbif_pa_data$sampleid
  gbif_pa_data$sampleid <- NULL
  
  
  zeta_orders <- c(2,10)
  
  Zeta_msgdm <- list()
  
  n_loop <- 1
  
  #Zeta.varpart Pairwise & n=10 similarity
  for (zeta_order in zeta_orders) {
    env_data <- c('WUI','WHP','FHSZ')
    for (n in 1:length(env_data)){
      env_combinations <- combn(env_data,n, simplify=FALSE)
      for (i in 1:length(env_combinations)){
        current_env_vars <- env_combinations[[i]]
        if (n == 1){
          current_env_data <- data.frame(st_drop_geometry(gbif_pa_data[,current_env_vars]))
          complete_env <- complete.cases(current_env_data)
          current_env_data <- data.frame(current_env_data[complete_env,])
          colnames(current_env_data) <- current_env_vars
          current_spec_data <- gbif_pa_data[,!(names(gbif_pa_data) %in% c('geometry','WUI','WHP','FHSZ'))]
          current_spec_data <- current_spec_data[complete_env,]
          current_xy_data <- gbif_pa_data[complete_env,]$geometry
          
        } else {
          current_env_data <- st_drop_geometry(gbif_pa_data[,current_env_vars])
          current_spec_data <- gbif_pa_data[,!(names(gbif_pa_data) %in% c('geometry','WUI','WHP','FHSZ'))]
          current_xy_data <- gbif_pa_data$geometry
        }
        #MS-GDM needs at least one numeric value. If NO column is numeric, designate the first column as numeric.
        if (!any(sapply(current_env_data, is.numeric))) {
          current_env_data[[1]] <- as.numeric(as.factor(current_env_data[[1]]))
          print(paste("Forced column", names(current_env_data)[1], "to numeric for MS-GDM."))
        }
        current_vars <- paste(current_env_vars,collapse=', ')
        print(paste("Current Var:",current_vars))
        Zeta_msgdm[[n_loop]] <- tryCatch({
          #Reset seed immediately as zetadiversity functions have their own montecarlo components.
          set.seed(seed_n)
          zeta_variation <- Zeta.varpart(
            Zeta.msgdm(
              data.spec=current_spec_data,
              data.env=current_env_data,
              xy=st_coordinates(current_xy_data),
              order=zeta_order
            ) 
          )
          print(paste0("For n:",n_threshold," & zeta order:",zeta_order))
          print(zeta_variation)
          print("a: envrionmental factors, b: distance, c: location, d:unexplained")
          # Build a data frame row
          data.frame(
            zeta_order = zeta_order,
            n_threshold = n_threshold, 
            env_vars = current_vars,
            status = "Success",
            error_msg = NA,
            # Note: Adjust extraction below based on how Zeta.varpart specifically names its output
            a_env = zeta_variation$`Adjusted Rsq`[4], 
            b_dist = zeta_variation$`Adjusted Rsq`[5],
            c_loc = zeta_variation$`Adjusted Rsq`[6],
            d_unexp = zeta_variation$`Adjusted Rsq`[7],
            stringsAsFactors = FALSE
          )
        },error=function(e){# Build a failed data frame row to capture the error safely
          data.frame(
            zeta_order = zeta_order,
            n_threshold = n_threshold,
            env_vars = current_vars,
            status = "Error",
            error_msg = conditionMessage(e),
            a_env = NA, b_dist = NA, c_loc = NA, d_unexp = NA,
            stringsAsFactors = FALSE
          )
        })
        n_loop <- n_loop + 1
      }
    }
  }
  
  all_Zeta_msgdm <- do.call(rbind, Zeta_msgdm)
  write.csv(all_Zeta_msgdm, paste("zeta_msgdm_results_n_",n_threshold,".csv", sep=""), row.names = FALSE)
  
  #Check decline
  zeta_dec <- Zeta.decline.ex(gbif_pa_data[, !(names(gbif_pa_data) %in% c('geometry','WUI', 'WHP', 'FHSZ'))], orders = 1:max(zeta_orders))
  Plot.zeta.decline(zeta_dec)
  #Decline is too high to run varpart on 3 or higher orders on windows
}


########
# Diversity measurements
########

#iNext ChaoRichness
# ChaoRichness(presabsdf,datatype='incidence_raw')
# iNEXT(presabsdf,q=0,datatype='incidence_raw')


########
# Correlation & Data Analysis:
########


#MAKE LEAFLET MAP OF ALL N 20,50,100 SAMPLES
#Plotting gbif data with CA Counties
#Plot GBIF Data & CA Counties
ggplot() +
  geom_sf(data=Study_counties, fill = 'lightgrey',color='black') +
  geom_sf(data=gbif_spatial, size=1,alpha=0.7, color='orange') +
  theme_minimal()

# Group by a unique sample ID based on time & location
# Remove undersampled Sites
gbif_grouped <- gbif_spatial %>% 
  group_by(eventDate,as.character(geometry)) %>% 
  mutate(sampleid = cur_group_id()) %>% 
  filter(n() >= 20) %>% ungroup()

# All data is from the Center for Biodiversity Genomics, in combination of samples containing at least 20 separate taxa, -> inferred to be eDNA samples
print(unique(gbif_grouped$institutionCode))

# Set FHSZ as a factor
gbif_grouped$FHSZ <- as.factor(gbif_grouped$FHSZ)

# Create presence/absence data
gbif_taxa_count <- gbif_grouped[,c('sampleid','taxonKey')]
#3,350 unique sample+species+location -> only presence noted, not number of times
gbif_taxa_count <- gbif_taxa_count[!duplicated(gbif_taxa_count),]
#Set to 1 i.e. species present in this sample.
gbif_taxa_count <- gbif_taxa_count %>%
  dplyr::group_by(sampleid,taxonKey) %>%
  dplyr::mutate(taxonCount = n()) %>%
  ungroup()
#Take the presence data, and add absence for every sampleid
gbif_pres_abs <- gbif_taxa_count %>%
  pivot_wider(
    names_from = taxonKey,
    values_from = taxonCount,
    values_fill = 0
  )


#Create dataframe
gbif_pres_abs <- as.data.frame(gbif_pres_abs)

#Get values for WHP, FHSZ, and WUI for each sample
gbif_data_lookup <- gbif_grouped %>%
  select(sampleid,WUI,WHP,FHSZ) %>%
  distinct() %>%
  st_drop_geometry()

#Join the data to each sampleid in gbif_pres_abs
gbif_pa_data <- gbif_pres_abs %>% left_join(gbif_data_lookup, by='sampleid')
gbif_pa_data$WUI  <- as.factor(gbif_pa_data$WUI)
gbif_pa_data$FHSZ <- as.factor(gbif_pa_data$FHSZ)
rownames(gbif_pa_data) <- gbif_pa_data$sampleid
gbif_pa_data$sampleid <- NULL


#Check for Correlations between FHSZ & WUI, WHP & WUI

#Mode function for most common categorical data
get_mode <- function(v) {
  #Remove any n/a
  uniq_v <- unique(na.omit(v))
  
  #return the most common value
  uniq_v[which.max(tabulate(match(v,uniq_v)))]
}

#Collapse gbif_grouped into mean (WHP) or mode (WUI,FHSZ) values for each unique sampleid
gbif_ANOVA <- gbif_grouped %>%
  st_drop_geometry() %>%
  group_by(sampleid) %>%
  summarise(
    #Continuous WHP: Mean is appropriate given single spatial coordinate
    WHP = mean(WHP,na.rm = TRUE),
    
    #Mode WUI & FHSZ: Mode is appropriate given categorical data & single spatial coordinate
    WUI = as.factor(get_mode(WUI)),
    FHSZ = as.factor(get_mode(FHSZ))
  ) %>%
  ungroup()

#Drop FHSZ (contain some n/a)
gbif_ANOVA_WUIWHP <- gbif_ANOVA %>% dplyr::select(WHP,WUI) %>% filter(!is.na(WHP) & !is.na(WUI))

#Drop n/a (FHSZ measures)
gbif_ANOVA_ALL <- gbif_ANOVA %>% filter(!is.na(WHP) & !is.na(WUI) & !is.na(FHSZ))

#ANOVA Analysis for All the data, How does WUI & FHSZ explain WHP
print(paste('WHP ~ WUI + FHSZ Analysis N:',nrow(gbif_ANOVA_ALL)))
combined_model <- aov(WHP ~ WUI + FHSZ, data=gbif_ANOVA_ALL)

print(summary(combined_model))
# Very Significant that WUI is a predictor of WHP
# Not significant that FHSZ is a predictor of WHP
print(TukeyHSD(combined_model))
# 2 & 3 not significantly different

#ANOVA Analysis for WHP & WUI
print(paste('WHP ~ WUI + Analysis N:',nrow(gbif_ANOVA_WUIWHP)))
partial_model <- aov(WHP ~ WUI, data=gbif_ANOVA_WUIWHP)

print(summary(partial_model))
# Very Significant that WUI is a predictor of WHP
print(TukeyHSD(partial_model))
# 5 & 1 are significantly different, and 5 is much riskier for WHP
# 8 & 1 are significantly different, and 1 is riskier for WHP
# 5 & 8 are significantly different, and 5 is much riskier for WHP
# Significant Risk 5->1->8 FSW Wild->FSW Intermix->Other(Water etc.)
# Wildlands more significant than Other for WHP

#Alpha Diversity comparison


#Combining FIRE Data: Do a FAMD to create an equivalent of PCA -> Single 'Fire Risk' statistic -> RF Model? -> Cluster

#1  Moderate. Lower level of wildfire hazard relative to other zones but still subject to wildfire behavior conditions.
#2	High. Elevated fire hazard due to fuels, terrain, and fire weather conditions.
#3	Very High. Highest level of wildfire hazard; areas most prone to wildfire spread and intensity.

#0  non-WUI / background / no classification
#1	Forest/Shrubland/Wetland-dominated Intermix WUI — areas where buildings and wildland vegetation are intermixed and the wildland component is forest/shrub/wetland. 
#2	Forest/Shrubland/Wetland-dominated Interface WUI — buildings are near large patches of forest/shrub/wetland vegetation (but not mixed within). 
#3	Grassland-dominated Intermix WUI — buildings intermingle with grasslands. 
#4	Grassland-dominated Interface WUI — buildings near large grassland areas. 
#5	Non-WUI: Forest/Shrubland/Wetland-dominated — wildland vegetation area not classified as WUI. 
#6	Non-WUI: Grassland-dominated — grassland not in WUI. 
#7	Non-WUI: Urban — built/urban land not in a WUI context. 
#8	Non-WUI: Other — other land types (e.g., water, bare land) outside WUI definitions.

#PLOTTING

#Plot WUI from gbif_grouped
ggplot() +
  geom_sf(data=Study_counties, fill = 'lightgrey',color='black') +
  geom_sf(data=gbif_grouped, size=2,alpha=0.7, color=gbif_grouped$WUI) +
  theme_minimal()

#Examine distribution of WUI indecies across gbif_grouped samples
WUI_gbif_grouped <- gbif_grouped %>% st_drop_geometry() %>% count(WUI,name='WUI_Count')

#Plot WUI Count from gbif_grouped
ggplot(WUI_gbif_grouped, aes(x=reorder(WUI,-WUI_Count),y=WUI_Count)) +
  geom_bar(stat='identity',fill='lightgrey', color='black') +
  geom_text(aes(label=WUI_Count),vjust = -0.5) +
  theme_minimal() +
  labs(
    title = 'Distribution of the WUI Types',
    x = 'WUI Index',
    y = 'Frequency'
  )
# 1-Intermix F/S/W, 5-Wildlands F/S/W, 8-Water/Bare, 6-Wildlands G
rm(WUI_gbif_grouped)

#Plot WHP from gbif_grouped
ggplot() +
  geom_sf(data=Study_counties, fill = 'lightgrey',color='black') +
  geom_sf(data=gbif_grouped, size=2,alpha=0.7, color=gbif_grouped$WHP) +
  theme_minimal()

#Examine distribution of WUI indecies across gbif_grouped samples (Binned into log)
WHP_gbif_grouped <- gbif_grouped %>% st_drop_geometry() %>% 
  mutate(
    WHP_Band = cut(WHP,breaks=2^(5:12),include.lowest = TRUE, dig.lab = 2)
  ) %>%
  count(WHP_Band,name='WHP_Count') %>%
  filter(!is.na(WHP_Band))

#Plot WUI Count from gbif_grouped
ggplot(WHP_gbif_grouped, aes(x=reorder(WHP_Band,-WHP_Count),y=WHP_Count)) +
  geom_bar(stat='identity',fill='lightgrey', color='black') +
  geom_text(aes(label=WHP_Count),vjust = -0.5) +
  theme_minimal() +
  labs(
    title = 'Distribution of the WHP Types',
    x = 'Log of WHP Index',
    y = 'Frequency'
  )
rm(WHP_gbif_grouped)

#Plot ANOVAs

ggplot(gbif_ANOVA_ALL, aes(x = WUI, y = WHP, fill=WUI)) +
  geom_boxplot(show.legend = TRUE) +
  labs(title = 'WHP by WUI',y='WHP',x='WUI Code') +
  theme_minimal()

ggplot(gbif_ANOVA_ALL, aes(x = FHSZ, y = WHP, fill=WUI)) +
  geom_boxplot(show.legend = TRUE) +
  labs(title = 'WHP by FHSZ',y='WHP',x='FHSZ Code') +
  theme_minimal()

ggplot(gbif_ANOVA_WUIWHP, aes(x = WUI, y = WHP, fill=WUI)) +
  geom_boxplot(show.legend = TRUE) +
  labs(title = 'WHP by WUI',y='WHP',x='WUI Code') +
  theme_minimal()