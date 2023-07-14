
library(tidyverse)
library(here)
library(raster)
library(sf)
library(fasterize)

# 2022 Imagery
img_2022_05_06a <- stack(here::here("..","..","..","planet","2023_floods","santa_clarita","santa_clarita_2022_05","files","20220506_174231_39_2421_3B_AnalyticMS_SR_8b_harmonized_clip.tif"))
img_2022_05_06b <- stack(here::here("..","..","..","planet","2023_floods","santa_clarita","santa_clarita_2022_05","files","20220506_181309_03_2475_3B_AnalyticMS_SR_8b_harmonized_clip.tif"))
img_2022_05_12a <- stack(here::here("..","..","..","planet","2023_floods","santa_clarita","santa_clarita_2022_05","files","20220512_174220_92_2440_3B_AnalyticMS_SR_8b_harmonized_clip.tif"))
img_2022_05_12b <- stack(here::here("..","..","..","planet","2023_floods","santa_clarita","santa_clarita_2022_05","files","20220512_174223_22_2440_3B_AnalyticMS_SR_8b_harmonized_clip.tif"))
img_2022_05_12c <- stack(here::here("..","..","..","planet","2023_floods","santa_clarita","santa_clarita_2022_05","files","20220512_181249_85_248b_3B_AnalyticMS_SR_8b_harmonized_clip.tif"))
img_2022_05_12d <- stack(here::here("..","..","..","planet","2023_floods","santa_clarita","santa_clarita_2022_05","files","20220512_181252_15_248b_3B_AnalyticMS_SR_8b_harmonized_clip.tif"))
# 2022 Composite
composite_2022 <- (img_2022_05_06a + 
             img_2022_05_12b +
             img_2022_05_12a +
             img_2022_05_12b +
             img_2022_05_12c +
             img_2022_05_12d) / 6
writeRaster(composite_2022, here::here("imagery","santa_clarita_2022_may_composite.tif"), overwrite=TRUE)

# 2023 Imagery
img_2023_05_08 <- stack(here::here("..","..","..","planet","2023_floods","santa_clarita","santa_clarita_2023_05","files","20230508_174407_33_24b5_3B_AnalyticMS_SR_8b_harmonized_clip.tif"))
img_2023_05_12a <- stack(here::here("..","..","..","planet","2023_floods","santa_clarita","santa_clarita_2023_05","files","20230512_174110_97_241b_3B_AnalyticMS_SR_8b_harmonized_clip.tif"))
img_2023_05_12b <- stack(here::here("..","..","..","planet","2023_floods","santa_clarita","santa_clarita_2023_05","files","20230512_182231_19_2254_3B_AnalyticMS_SR_8b_harmonized_clip.tif"))
# 2023 Composite
composite_2023 <- (img_2023_05_08 + 
                     img_2023_05_12a +
                     img_2023_05_12b) / 3
writeRaster(composite_2023, here::here("imagery","santa_clarita_2023_may_composite.tif"), overwrite=TRUE)

# Get the spectral angle between two raster images
spectralAngle <- function(img_1, img_2)
{
  num_bands <- dim(img_1)[[3]]
  if(num_bands != dim(img_2)[[3]])
    return(-1)
  
  sum_term <- img_1[[1]]*0
  magnitude_1 <- img_1[[1]]*0 + 1
  magnitude_2 <- img_2[[1]]*0 + 1
  
  for(band in 1:num_bands)
  {
    sum_term <- sum_term + img_1[[band]]*img_2[[band]]
    magnitude_1 <- magnitude_1 + img_1[[band]]^2
    magnitude_2 <- magnitude_2 + img_2[[band]]^2
  }
  magnitude_1 <- magnitude_1^0.5
  magnitude_2 <- magnitude_2^0.5
  
  return(acos(sum_term / magnitude_1 / magnitude_2))
}

spectral_angle_change <- spectralAngle(composite_2022, composite_2023)
writeRaster(spectral_angle_change, here::here("imagery","2023_2022_may_angle_change.tif"))

# Load field vegetation data 
veg_polygons <- st_read(here::here("..","..","vegetation_surveys","Veg_Polygons_SCR_2022_Stillwater","SCR_VegMap.shp"))
# Add numeric class for Arundo level
arundo_levels <- 1:8
names(arundo_levels) <- c("<1",  
                          "1-5",
                          "5-10",
                          "10-25",
                          "25-50", 
                          "50-75", 
                          "75-95",
                          ">95")
arundo_means <- c(0.5, 3, 7.5, 17.5, 37.5, 62.5, 85, 97.5) / 100
names(arundo_means) <- names(arundo_levels)
veg_polygons$ARDO_lvl <- arundo_levels[veg_polygons$ARDO_class]
veg_polygons$ARDO_mean <- arundo_means[veg_polygons$ARDO_class]
# Reproject data to raster CRS
veg_polygons <- st_transform(veg_polygons, st_crs(spectral_angle_change))
# Rasterize arundo level
arundo_lvl <- fasterize(veg_polygons, spectral_angle_change, field = "ARDO_lvl")
arundo_mean <- fasterize(veg_polygons, spectral_angle_change, field = "ARDO_mean")

# Generate NDVI imagery
ndvi_2022 <- (composite_2022[[8]] - composite_2022[[6]]) / (composite_2022[[8]] + composite_2022[[6]])
ndvi_2023 <- (composite_2023[[8]] - composite_2023[[6]]) / (composite_2023[[8]] + composite_2023[[6]])

writeRaster(ndvi_2022, here::here("imagery","santa_clarita_ndvi_2022.tif"), overwrite=TRUE)
writeRaster(ndvi_2023, here::here("imagery","santa_clarita_ndvi_2023.tif"), overwrite=TRUE)

# Get change in NDVI from 2022 to 2023
ndvi_change <- ndvi_2023 - ndvi_2022
ndvi_change[is.na(arundo_lvl)] <- NA
writeRaster(ndvi_change, here::here("imagery","ndvi_change_2022_2023.tif"), overwrite=TRUE)

# Output several products in one image for convenience
output_data_stack <- stack(spectral_angle_change, arundo_mean, ndvi_2022, ndvi_2023)
#   filter to remove areas NOT mapped in vegetation polygons (outside river channel)
output_data_stack[is.na(output_data_stack[[2]])] <- NA
names(output_data_stack) <- c("spectral_change",
                              "arundo_lvl",
                              "ndvi_2022",
                              "ndvi_2023")
writeRaster(output_data_stack, here::here("imagery","santa_clarita_change.tif"), overwrite=TRUE)

# Extract change statistics for each polygon
veg_polygons$UID <- as.numeric(veg_polygons$UID)
polygon_label_raster <- fasterize(veg_polygons, ndvi_change, field = "UID")
# Get label and raster data in one stack
all_polygon_data_stack <- stack(polygon_label_raster, 
                                output_data_stack)
names(all_polygon_data_stack) <- c("UID",
                                   names(output_data_stack))
all_polygon_data <- as.data.frame(all_polygon_data_stack) %>% 
  drop_na()
# Summarize spatial data in each polygon
polygon_data_summaries <- all_polygon_data %>% 
  group_by(UID) %>% 
  summarize(min_NDVI_2022 = min(ndvi_2022),
            mean_NDVI_2022 = mean(ndvi_2022),
            max_NDVI_2022 = max(ndvi_2022),
            min_NDVI_2023 = min(ndvi_2023),
            mean_NDVI_2023 = mean(ndvi_2023),
            max_NDVI_2023 = max(ndvi_2023),
            mean_ndvi_change = mean(ndvi_2023-ndvi_2022),
            min_sac = min(spectral_change),
            mean_sac = mean(spectral_change),
            max_sac = max(spectral_change))
# Add raster data to polygons after filtering for the polygons included in our subset
veg_polygons_output <- veg_polygons[veg_polygons$UID %in% unique(polygon_data_summaries$UID),]
veg_polygons_output <- merge(veg_polygons_output, polygon_data_summaries)
st_write(veg_polygons_output, here::here("imagery","santa_clarita_vegetation.gpkg"))

# Initial results, including 2 and 3-endmember models
#mesma_2022 <- stack(here::here("..","..","..","planet","2023_floods","santa_clarita","2022_may_composite.hdr_SMA_20230516T14H36M54S_shade_normalized.tif"))
#mesma_2023 <- stack(here::here("..","..","..","planet","2023_floods","santa_clarita","2023_may_composite.hdr_SMA_20230516T14H37M33S_shade_normalized.tif"))
# Final results, including only 3-endmember models
mesma_2022 <- stack(here::here("..","..","..","planet","2023_floods","santa_clarita","2022_may_composite.hdr_SMA_20230524T09H27M26S_shade_normalized.tif"))
mesma_2023 <- stack(here::here("..","..","..","planet","2023_floods","santa_clarita","2023_may_composite.hdr_SMA_20230524T10H09M36S_shade_normalized.tif"))

mesma_2022 <- mask(mesma_2022, veg_polygons_output)
mesma_2023 <- mask(mesma_2023, veg_polygons_output)

writeRaster(mesma_2022, here::here("..","..","..","planet","2023_floods","santa_clarita","2022_may_sma_masked_v2.tif"))
writeRaster(mesma_2023, here::here("..","..","..","planet","2023_floods","santa_clarita","2023_may_sma_masked_v2.tif"))

# Also do some masking for shade-normalized SMA data from other dates
#  Santa Clarita scenes
sc_mesma_dec <- stack(here::here("..","..","..","planet","2023_floods","santa_clarita","santa_clarita_2022_12_20/files/santa_clarita_2022_12_20.hdr_SMA_20230531T09H18M56S_shade_normalized.tif"))
sc_mesma_feb <- stack(here::here("..","..","..","planet","2023_floods","santa_clarita","santa_clarita_2023_02_19/files/santa_clarita_2023_02_19.hdr_SMA_20230531T09H16M48S_shade_normalized.tif"))
sc_mesma_mar <- stack(here::here("..","..","..","planet","2023_floods","santa_clarita","santa_claria_2023_03_07/files/santa_clarita_2023_03_07.hdr_SMA_20230531T09H20M14S_shade_normalized.tif"))
#  Saticoy scenes
sy_mesma_dec <- stack(here::here("..","..","..","planet","2023_floods","saticoy","saticoy_2022_12_23/files/saticoy_2022_12_23.hdr_SMA_20230531T10H46M18S_shade_normalized.tif"))
sy_mesma_feb <- stack(here::here("..","..","..","planet","2023_floods","saticoy","saticoy_2023_02_07/files/saticoy_2023_02_07.hdr_SMA_20230531T10H51M09S_shade_normalized.tif"))
sy_mesma_mar <- stack(here::here("..","..","..","planet","2023_floods","saticoy","saticoy_2023_03_07_psscene_analytic_8b_sr_udm2/files/saticoy_2023_03_07.hdr_SMA_20230531T10H56M11S_shade_normalized.tif"))

sc_mesma_dec <- mask(sc_mesma_dec, veg_polygons)
sc_mesma_feb <- mask(sc_mesma_feb, veg_polygons)
sc_mesma_mar <- mask(sc_mesma_mar, veg_polygons)
sy_mesma_dec <- mask(sy_mesma_dec, veg_polygons)
sy_mesma_feb <- mask(sy_mesma_feb, veg_polygons)
sy_mesma_mar <- mask(sy_mesma_mar, veg_polygons)

writeRaster(sc_mesma_dec, here::here("..","..","..","planet","2023_floods","santa_clarita","santa_clarita_2022_12_20/files/santa_clarita_2022_12_20_sma.tif"), overwrite=TRUE)
writeRaster(sc_mesma_feb, here::here("..","..","..","planet","2023_floods","santa_clarita","santa_clarita_2023_02_19/files/santa_clarita_2023_02_19_sma.tif"), overwrite=TRUE)
writeRaster(sc_mesma_mar, here::here("..","..","..","planet","2023_floods","santa_clarita","santa_claria_2023_03_07/files/santa_clarita_2023_03_07_sma.tif"), overwrite=TRUE)
writeRaster(sy_mesma_dec, here::here("..","..","..","planet","2023_floods","saticoy","saticoy_2022_12_23/files/saticoy_2022_12_23_sma.tif"), overwrite=TRUE)
writeRaster(sy_mesma_feb, here::here("..","..","..","planet","2023_floods","saticoy","saticoy_2023_02_07/files/saticoy_2023_02_07_sma.tif"), overwrite=TRUE)
writeRaster(sy_mesma_mar, here::here("..","..","..","planet","2023_floods","saticoy","saticoy_2023_03_07_psscene_analytic_8b_sr_udm2/files/santa_clarita_2023_03_07_sma.tif"), overwrite=TRUE)



# Change in fractional cover
gv_change <- mesma_2023[[1]] - mesma_2022[[1]]
npv_change <- mesma_2023[[2]] - mesma_2022[[2]]
soil_change <-mesma_2023[[3]] - mesma_2022[[3]]
fraction_change <- stack(gv_change, npv_change, soil_change)
writeRaster(fraction_change, here::here("..","..","..","planet","2023_floods","santa_clarita","2023_fraction_change_v2.tif"))

scoured <- (gv_change < 0) * (soil_change > 0) * (mesma_2023[[1]] < 0.05)
writeRaster(scoured, here::here("..","..","..","planet","2023_floods","santa_clarita","scoured.tif"), overwrite=TRUE)
killed <- (gv_change < 0) * (npv_change > 0) * (mesma_2023[[1]] < 0.05) * (scoured == 0)
writeRaster(killed, here::here("..","..","..","planet","2023_floods","santa_clarita","killed.tif"), overwrite=TRUE)
veg_increased <- (gv_change > 0) * (scoured == 0) * (killed == 0)
writeRaster(veg_increased, here::here("..","..","..","planet","2023_floods","santa_clarita","veg_increased.tif"), overwrite=TRUE)



# Get label and raster data in one stack
all_polygon_data_stack <- stack(polygon_label_raster, 
                                output_data_stack)
names(all_polygon_data_stack) <- c("UID",
                                   names(output_data_stack))
all_polygon_data <- as.data.frame(all_polygon_data_stack) %>% 
  drop_na()
# Summarize spatial data in each polygon
polygon_data_summaries <- all_polygon_data %>% 
  group_by(UID) %>% 
  summarize(min_NDVI_2022 = min(ndvi_2022),
            mean_NDVI_2022 = mean(ndvi_2022),
            max_NDVI_2022 = max(ndvi_2022),
            min_NDVI_2023 = min(ndvi_2023),
            mean_NDVI_2023 = mean(ndvi_2023),
            max_NDVI_2023 = max(ndvi_2023),
            mean_ndvi_change = mean(ndvi_2023-ndvi_2022),
            min_sac = min(spectral_change),
            mean_sac = mean(spectral_change),
            max_sac = max(spectral_change))



# Load in manually updated polygons (split based on field surveys)
veg_updated <- st_read(here::here("imagery","santa_clarita_veg_update.shp"))
# Update arundo percentage based on initial and updated values
veg_updated[!is.na(veg_updated$ARDO_lvl_2),]$ARDO_mean <- veg_updated[!is.na(veg_updated$ARDO_lvl_2),]$ARDO_lvl_2 / 100
veg_updated[is.na(veg_updated$scoured),]$scoured <- 0
veg_updated$polygon_id <- 1:nrow(veg_updated)
veg_updated$overall_area <- as.numeric(st_area(veg_updated))/4046.86
veg_updated$arundo_area <- veg_updated$overall_area * veg_updated$ARDO_mean
ardo_types <- c("Standing",
                "Knocked Down")
veg_updated$ardo_condition <- ardo_types[is.na(veg_updated$ardo_type) + 1]

veg_types <- c("Scrub",
               "Woodland",
               "Scrub",
               "Scrub",
               "Scrub",
               "Scrub",
               "Scrub",
               "Developed",
               "Herbaceous",
               "Scrub",
               "Scrub",
               "Riverwash",
               "Herbaceous",
               "Wetland",
               "Scrub",
               "Arundo",
               "Developed",
               "Scrub",
               "Scrub",
               "Riverwash",
               "Woodland",
               "Herbaceous",
               "Scrub",
               "Scrub",
               "Scrub",
               "Scrub",
               "Herbaceous")
names(veg_types) <- c( "Lepidospartum squamatum Shrubland Alliance",                                              
                       "Populus fremontii Forest Alliance",                                                       
                       "Artemisia californica Shrubland Alliance",                                                
                       "Artemisia californica - Salvia mellifera Shrubland Alliance",                             
                       "Baccharis pilularis Shrubland Alliance",                                                  
                       "Artemisia tridentata Shrubland Alliance",                                                 
                       "Atriplex canescens Shrubland Alliance",                                                   
                       "Developed",                                                                               
                       "Bromus (diandrus, hordeaceus) - Brachypodium distachyon Herbaceous Semi-Natural Alliance",
                       "Baccharis salicifolia Shrubland Alliance",                                                
                       "Eriogonum fasciculatum Shrubland Alliance",                                               
                       "Riverwash",                                                                               
                       "Non-native Grass and Forb Mapping Unit",                                                  
                       "Typha (angustifolia, domingensis, latifolia) Herbaceous Alliance",                        
                       "Sambucus nigra Shrubland Alliance",                                                       
                       "Phragmites australis - Arundo donax Herbaceous Semi-Natural Alliance",                    
                       "Disturbed",                                                                               
                       "Buckwheat",                                                                               
                       "Tamarix spp. Shrubland Semi-Natural Alliance",                                            
                       "Barren",                                                                                  
                       "Willow",                                                                                  
                       "Annual Grasses and Forbs",                                                                
                       "Baccharis (Riparian)",                                                                    
                       "Riparian Mixed Shrub",                                                                    
                       "Scalebroom",                                                                              
                       "Riversidean Alluvial Scrub",                                                              
                       "Brachypodium distachyon Herbaceous Semi-Natural Alliance")
veg_updated$VegClass <- veg_types[veg_updated$VegType]
veg_updated[veg_updated$scoured == 1,]$VegClass = "Riverwash"
st_write(veg_updated,
         here::here("imagery","final_veg_santa_clarita.shp"))

write_csv(as.data.frame(veg_updated) %>% 
            arrange(VegClass, -ARDO_mean, polygon_id) %>%
            dplyr::select(VegClass, ARDO_mean, polygon_id, overall_area, arundo_area, ardo_condition),
          here::here("arundo_table.csv"))

write_csv(as.data.frame(veg_updated) %>% 
            group_by(VegClass) %>%
            summarize(ARDO_mean = sum(arundo_area)/sum(overall_area),
                      overall_area = sum(overall_area),
                      arundo_area = sum(arundo_area)) %>%
            arrange(VegClass) %>%
            dplyr::select(VegClass, ARDO_mean, overall_area, arundo_area),
          here::here("arundo_table_class_avg.csv"))





