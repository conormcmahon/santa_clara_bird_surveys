
library(raster)
library(sf)
library(tidyverse)

before_floods <- stack(here::here("..","..","..","planet","2023_floods","saticoy","saticoy_2022_12_23","files","20221223_175927_32_2262_3B_AnalyticMS_SR_8b_harmonized_clip.tif"))
between_floods <- stack(here::here("..","..","..","planet","2023_floods","saticoy","saticoy_2023_02_07","files","20230208_183436_93_2402_3B_AnalyticMS_SR_8b_harmonized_clip.tif"))
after_floods <- stack(here::here("..","..","..","planet","2023_floods","saticoy","saticoy_2023_03_06","files","20230306_181619_70_2488_3B_AnalyticMS_SR_8b_harmonized_clip.tif"))

before_ndvi <- (before_floods[[8]]-before_floods[[6]]) / (before_floods[[8]]+before_floods[[6]])
between_ndvi <- (between_floods[[8]]-between_floods[[6]]) / (between_floods[[8]]+between_floods[[6]])
after_ndvi <- (after_floods[[8]]-after_floods[[6]]) / (after_floods[[8]]+after_floods[[6]])
writeRaster(before_ndvi, here::here("..","..","..","planet","2023_floods","saticoy","before_ndvi.tif"), overwrite=TRUE)
writeRaster(between_ndvi, here::here("..","..","..","planet","2023_floods","saticoy","between_ndvi.tif"), overwrite=TRUE)
writeRaster(after_ndvi, here::here("..","..","..","planet","2023_floods","saticoy","after_ndvi.tif"), overwrite=TRUE)

before_ndwi <- (before_floods[[4]]-before_floods[[8]]) / (before_floods[[4]]+before_floods[[8]])
between_ndwi <- (between_floods[[4]]-between_floods[[8]]) / (between_floods[[4]]+between_floods[[8]])
after_ndwi <- (after_floods[[4]]-after_floods[[8]]) / (after_floods[[4]]+after_floods[[8]])
writeRaster(before_ndwi, here::here("..","..","..","planet","2023_floods","saticoy","before_ndwi.tif"), overwrite=TRUE)
writeRaster(between_ndwi, here::here("..","..","..","planet","2023_floods","saticoy","between_ndwi.tif"), overwrite=TRUE)
writeRaster(after_ndwi, here::here("..","..","..","planet","2023_floods","saticoy","after_ndwi.tif"), overwrite=TRUE)

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

january_sma_change <- spectralAngle(before_floods[[2:8]], between_floods[[2:8]])
february_sma_change <- spectralAngle(between_floods[[2:8]], after_floods[[2:8]])
winter_sma_change <- spectralAngle(before_floods[[2:8]], after_floods[[2:8]])
writeRaster(january_sma_change, here::here("..","..","..","planet","2023_floods","saticoy","jan_sma_change.tif"), overwrite=TRUE)
writeRaster(february_sma_change, here::here("..","..","..","planet","2023_floods","saticoy","feb_sma_change.tif"), overwrite=TRUE)
writeRaster(winter_sma_change, here::here("..","..","..","planet","2023_floods","saticoy","total_sma_change.tif"), overwrite=TRUE)

# Load polygons from vegetation surveys
veg_polygons <- st_read(here::here("vegetation_surveys","Veg_Polygons_All","Veg_Polygons_All.shp"))
veg_polygons <- st_transform(veg_polygons, st_crs(january_sma_change))

# Add change information to polygons
veg_polygons$jan_change <- raster::extract(x=january_sma_change, y=veg_polygons[,1], fun=mean, na.rm=TRUE)
veg_polygons$feb_change <- raster::extract(x=february_sma_change, y=veg_polygons[,1], fun=mean, na.rm=TRUE)
veg_polygons$total_change <- raster::extract(x=winter_sma_change, y=veg_polygons[,1], fun=mean, na.rm=TRUE)
veg_polygons <- veg_polygons %>%
  rowwise() %>%
  mutate(max_change = max(jan_change, feb_change, total_change),
         flood_changed = max_change > 0.2)
st_write(veg_polygons %>% filter(max_change > 0.2), here::here("vegetation_surveys","veg_flood_change","veg_altered.shp"), delete_dsn=TRUE)
st_write(veg_polygons %>% filter(max_change <= 0.2), here::here("vegetation_surveys","veg_flood_change","veg_unchanged.shp"), delete_dsn=TRUE)
st_write(veg_polygons, here::here("vegetation_surveys","veg_flood_change","veg_flood_impacts.shp"), delete_dsn=TRUE)

veg_type_colors <- c("goldenrod1","red3","chartreuse","bisque4","cyan3","blue","green4")
veg_types <- sort(unique(veg_polygons$Type))

# Visualize Changes by Veg Type
ggplot(veg_polygons) + 
  geom_histogram(aes(x=max_change, fill=Type)) + 
  facet_wrap(~Type) + 
  scale_fill_manual(values=veg_type_colors,
                     labels=veg_types) + 
  theme(legend.position="none") + 
  xlab("Spectral Angle Change Post-flood") + 
  ylab("Frequency") + 
  geom_vline(xintercept=0.2, col="black", linetype="dashed")
# I used a threshold of 0.2 SMA for changed polygons based on a visual check - that seemed to be an important breakpoint 

