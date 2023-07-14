
# Load libraries
library(terra)
library(sf)
library(here)
library(tidyverse)

# Load source data
planet_2022 <- terra::rast(here::here("..","..","..","planet","2023_floods","saticoy","saticoy_2022_07_07","PSScene","20220707_181608_12_2499_3B_AnalyticMS_SR_8b_harmonized_clip.tif"))
planet_2023 <- terra::rast(here::here("..","..","..","planet","2023_floods","saticoy","saticoy_2023_07_13","PSScene","20230706_174627_23_24bc_3B_AnalyticMS_SR_8b_harmonized_clip.tif"))
saticoy_veg <- st_read(here::here("..","..","vegetation_surveys","Veg_Polygons_All","veg_polygons_all.shp"))
saticoy_veg <- st_transform(saticoy_veg, st_crs(planet_2023))

# Generate NDVI greenness imagery
ndvi_2022 <- (planet_2022[[8]]-planet_2022[[6]]) / (planet_2022[[8]]+planet_2022[[6]])
ndvi_2023 <- (planet_2023[[8]]-planet_2023[[6]]) / (planet_2023[[8]]+planet_2023[[6]])

# Reduction in NDVI after flood
ndvi_loss <- (ndvi_2023 - ndvi_2022)

# Areas with a lot of initial veg and a lot of loss
substantial_veg_loss <- (ndvi_loss < -0.15) * (ndvi_2022 > 0.5)

# Get average post-flood loss of NDVI in a given polygon 
getNDVILoss <- function(ind)
{
  polygon <- saticoy_veg[ind,]
  ndvi_change_values <- terra::extract(ndvi_loss, polygon)
  polygon$mean_loss <- mean(ndvi_change_values$nir, na.rm=TRUE)
  polygon$max_loss <- max(ndvi_change_values$nir, na.rm=TRUE)
  polygon$min_loss <- min(ndvi_change_values$nir, na.rm=TRUE)
  return(polygon)
}

new_polygons <- lapply(1:nrow(saticoy_veg), getNDVILoss)
new_polygons <- bind_rows(new_polygons)

# Get area, in acres
new_polygons$area <- as.numeric(st_area(new_polygons)) * 0.000247105

getArundoFracBins <- function(data)
{
  bin_low_ends <- c(0, 1, 5, 20)
  bin_high_ends <- c(1, 5, 20, 101)
  
  fillArundoBin <- function(index)
  {
    data_in_bin <- data %>% 
      filter(ArundoPerc >= bin_low_ends[[index]],
             ArundoPerc < bin_high_ends[[index]])
    if(nrow(data_in_bin) == 0)
    {
      return(data.frame(num_polygons = 0,
                        total_area = 0,
                        bin_low_end = bin_low_ends[[index]],
                        bin_high_end = bin_high_ends[[index]]))
    }
    return(data_in_bin %>% 
             summarize(num_polygons = n(),
                       total_area = sum(area),
                       bin_low_end = bin_low_ends[[index]],
                       bin_high_end = bin_high_ends[[index]]))
  }
  
  new_data <- lapply(1:length(bin_low_ends),
                               fillArundoBin)
  new_data <- bind_rows(new_data)
  
  return(new_data)
}

# Pre-flooding Arundo Area
print("")
print(paste("Pre-flood total area of Arundo is: ", sum(new_polygons$area), sep=""))
plot(new_polygons[,"ArundoPerc"])
print("")
print("Area overall, by Arundo density:")
getArundoFracBins(new_polygons)

# Flood-impacted Arundo Area
print("")
print("Arundo area by veg type, in flood-affected areas:")
getArundoFracBins(new_polygons %>% filter(mean_loss < -0.15))
print("")
print("Arundo area by veg type, in unaffected areas:")
getArundoFracBins(new_polygons %>% filter(mean_loss >= -0.15))






