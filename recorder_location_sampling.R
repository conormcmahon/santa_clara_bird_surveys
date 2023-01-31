

library(sf)
library(raster)
library(tidyverse)

veg_polygons <- st_read(here::here("vegetation_surveys","Veg_Polygons_All","Veg_Polygons_All.shp"))
veg_types <- unique(veg_polygons$Type)

set.seed(1)

# Function to randomly sample points among polygons of a given type 
generatePointsByType <- function(veg_type, veg_edge_buffer, inter_unit_distance, min_arundo, max_arundo) 
{
  # Get stats on polygons in a certain veg type
  print("")
  polygon_subset <- veg_polygons %>% filter(Type == veg_type,
                                            ArundoPerc > min_arundo,
                                            ArundoPerc < max_arundo) 
  print(paste("Working on polygons of type                ", veg_type,sep=""))
  print(paste("Number of Polygons:                        ", nrow(polygon_subset),sep=""))
  print(paste("Total area:                                ", round(sum(polygon_subset$Acres)*4046.86,1)," m^2", sep=""))
  print(paste("Average Arundo Pct:                        ", mean(polygon_subset$ArundoPerc * polygon_subset$Acres) / sum(polygon_subset$ArundoPerc),sep=""))
  print(paste("Min Arundo Pct:                            ", min(polygon_subset$ArundoPerc),sep=""))
  print(paste("Max Arundo Pct:                            ", max(polygon_subset$ArundoPerc),sep=""))
  print("")
  
  # Merge all polygons into one big polygon for that veg type
  one_big_polygon <- st_sf(st_union(polygon_subset))
  print(paste("Total area:                                ", round(st_area(one_big_polygon),1)," m^2", sep=""))
  
  # Buffer that polygon inwards, to prevent sampling of other habitat types
  one_big_polygon_subset <- st_buffer(one_big_polygon, -veg_edge_buffer)
  print(paste("After buffer, total sample-able area is:   ", round(st_area(one_big_polygon_subset),1)," m^2", sep=""))
  
  points <- st_sample(one_big_polygon_subset, 1)
  subset_polygon <- st_difference(one_big_polygon_subset, st_buffer(points, inter_unit_distance))
  # Repeatedly draw out random point locations, removing data to prevent overlap in sampling area
  while(nrow(subset_polygon) > 0)
  {
    new_point <- st_sample(subset_polygon, 1)
    points <- c(points, new_point)
    subset_polygon <- st_difference(subset_polygon, st_buffer(new_point, inter_unit_distance))
  }
  print(paste("Retrieved a total of ", length(points), " points within this class.", sep=""))
  return(points)
}


# Generate random points in each veg type, stratified at least *inter_unit_distance* m apart, and each at least *veg_edge_buffer* m from a veg-type boundary
generateAllPointsForParameterSet <- function(veg_edge_buffer, inter_unit_distance)
{
  # Get all points for each veg type
  woodland_high_arundo <- generatePointsByType("Woodland", veg_edge_buffer, inter_unit_distance, 30, 100)
  woodland_low_arundo <- generatePointsByType("Woodland", veg_edge_buffer, inter_unit_distance, 0, 30)
  riparian_scrub <- generatePointsByType("Riparian scrub", veg_edge_buffer, inter_unit_distance, 0, 100)
  typha_wetland <- generatePointsByType("Typha wetland", 20, inter_unit_distance, 0, 100)
  sage_scrub <- generatePointsByType("Sage scrub", veg_edge_buffer, inter_unit_distance, 0, 100)
  river_wash <- generatePointsByType("Riverwash", veg_edge_buffer, inter_unit_distance, 0, 100)
  grasses_and_forbs <- generatePointsByType("Grasses and forbes", veg_edge_buffer, inter_unit_distance, 0, 100)
  non_native <- generatePointsByType("Non-native", veg_edge_buffer, inter_unit_distance, 0, 100)
  
  all_points <- st_sf(c(woodland_high_arundo,
                        woodland_low_arundo,
                        riparian_scrub,
                        typha_wetland,
                        sage_scrub,
                        river_wash,
                        grasses_and_forbs,
                        non_native))
  
  return(st_join(all_points, veg_polygons))
}

# Print some summary data about a set of points 
printPointSummary <- function(points, name)
{
  print(as.data.frame(points %>% 
          group_by(Type) %>%
          summarize(count = n(),
                    min_arundo = min(ArundoPerc),
                    max_arundo = max(ArundoPerc),
                    mean_arundo = mean(ArundoPerc),
                    frac_with_surf_water = sum(SurfaceWat=="Y")/n())))
  st_write(points, here::here("point_examples",paste(name,".shp",sep="")))
}

# First, ensure each point is embedded within 50 m of contiguous habitat, and at least 50 m from each other point
example_1 <- generateAllPointsForParameterSet(50, 50)
printPointSummary(example_1)
# next, ensure each point is embedded within 50 m of contiguous habitat, and at least 100 m from each other point
example_2 <- generateAllPointsForParameterSet(50, 100)
printPointSummary(example_2)
# next, ensure each point is embedded within 50 m of contiguous habitat, and at least 100 m from each other point
example_3 <- generateAllPointsForParameterSet(50, 150)
printPointSummary(example_3)
# next, ensure each point is embedded within 50 m of contiguous habitat, and at least 100 m from each other point
example_4 <- generateAllPointsForParameterSet(100, 100)
printPointSummary(example_4)
# next, ensure each point is embedded within 50 m of contiguous habitat, and at least 100 m from each other point
example_5 <- generateAllPointsForParameterSet(150, 150)
printPointSummary(example_5)


