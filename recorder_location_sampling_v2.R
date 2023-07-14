

library(sf)
library(raster)
library(tidyverse)

veg_polygons <- st_read(here::here("vegetation_surveys","veg_flood_change","veg_flood_impacts.shp"))
veg_types <- unique(veg_polygons$Type)

set.seed(1)

# Function to randomly sample points among polygons of a given vegetation type 
#   Veg_Type - string for vegetation type (e.g. Woodland, Riverwash, etc.)
#   Veg_Edge_Buffer - min distance (m) a point can be from the boundary of a vegetation type
#   Inter_Unit_Distance - minimum distance between two recorders 
#   Min_Arundo - minimum Arundo percentage in polygons used 
#   Max_Arundo - maximum Arundo percentage in polygons used
generatePointsByType <- function(veg_type, veg_edge_buffer, inter_unit_distance, min_arundo, max_arundo, change_min, change_max) 
{
  # Print summary stats on polygons in a certain veg type
  print("")
  polygon_subset <- veg_polygons %>% filter(Type %in% veg_type,
                                            ArundoPerc > min_arundo,
                                            ArundoPerc <= max_arundo, 
                                            max_change < change_max,
                                            max_change >= change_min) 
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


# Apply above function across all vegetation types
# Returns a sf dataframe with all possible points randomly sampled in each veg type, given constraints
generateAllPointsForParameterSet <- function(veg_edge_buffer, inter_unit_distance, arundo_threshold, change_threshold)
{
  # Get all points for each veg type
  woodland_high_arundo_changed <- generatePointsByType(c("Woodland","Non-native"), veg_edge_buffer, inter_unit_distance, arundo_threshold, 100, change_threshold, 100)
  woodland_low_arundo_changed <- generatePointsByType(c("Woodland","Non-native"), veg_edge_buffer, inter_unit_distance, 0, arundo_threshold, change_threshold, 100)
  woodland_high_arundo_unchanged <- generatePointsByType(c("Woodland","Non-native"), veg_edge_buffer, inter_unit_distance, arundo_threshold, 100, 0, change_threshold)
  woodland_low_arundo_unchanged <- generatePointsByType(c("Woodland","Non-native"), veg_edge_buffer, inter_unit_distance, 0, arundo_threshold, 0, change_threshold)
  scrub_unchanged <- generatePointsByType(c("Sage scrub","Riparian scrub","Type wetland"), veg_edge_buffer, inter_unit_distance, 0, arundo_threshold, 0, change_threshold)
  scrub_changed <- generatePointsByType(c("Sage scrub","Riparian scrub","Type wetland"), veg_edge_buffer, inter_unit_distance, 0, arundo_threshold, 0, 0.2)
  
  # Combine point datasets, convert to sf dataframe
  all_points <- st_sf(c(woodland_high_arundo_changed,
                        woodland_low_arundo_changed,
                        woodland_high_arundo_unchanged,
                        woodland_low_arundo_unchanged,
                        scrub_unchanged,
                        scrub_changed))
  all_points_data <- st_join(all_points, veg_polygons)
  
  all_points_data$class_final <- as.character(all_points_data$Type)
  all_points_data[all_points_data$Type == "Woodland",]$class_final <- paste(c("Woodland","Scrub")[all_points_data[all_points_data$Type == "Woodland",]$class_final %in% c("Woodland","Non-native")],
                                                                            c(" High"," Low")[(all_points_data[all_points_data$Type == "Woodland",]$ArundoPerc > arundo_threshold)+1],
                                                                            c(" Unchanged"," Changed")[(all_points_data[all_points_data$Type == "Woodland",]$max_change > change_threshold)+1],
                                                                            sep="")
  return(all_points_data)
}


# Print some summary data about a set of points 
printPointSummary <- function(points, name)
{
  print(as.data.frame(points %>% 
                        group_by(class_final) %>%
                        summarize(count = n(),
                                  min_arundo = min(ArundoPerc),
                                  max_arundo = max(ArundoPerc),
                                  mean_arundo = mean(ArundoPerc),
                                  min_change = min(max_change),
                                  max_change = min(max_change),
                                  mean_change = min(max_change),
                                  frac_with_surf_water = sum(SurfaceWat=="Y")/n())))
  st_write(points, here::here("point_examples_v2",paste(name,".shp",sep="")), delete_dsn=TRUE)
}

# First, ensure each point is embedded within 50 m of contiguous habitat, and at least 50 m from each other point
example_1 <- generateAllPointsForParameterSet(50, 50, 30)
printPointSummary(example_1, "example_1")
# next, ensure each point is embedded within 50 m of contiguous habitat, and at least 100 m from each other point
example_2 <- generateAllPointsForParameterSet(50, 100, 30)
printPointSummary(example_2, "example_2")
# next, ensure each point is embedded within 50 m of contiguous habitat, and at least 150 m from each other point
example_3 <- generateAllPointsForParameterSet(50, 150, 30)
printPointSummary(example_3, "example_3")
# First, ensure each point is embedded within 50 m of contiguous habitat, and at least 50 m from each other point
example_1 <- generateAllPointsForParameterSet(25, 50, 30)
printPointSummary(example_1, "example_4")
# next, ensure each point is embedded within 50 m of contiguous habitat, and at least 100 m from each other point
example_2 <- generateAllPointsForParameterSet(25, 100, 30)
printPointSummary(example_2, "example_5")
# next, ensure each point is embedded within 50 m of contiguous habitat, and at least 150 m from each other point
example_3 <- generateAllPointsForParameterSet(25, 150, 30)
printPointSummary(example_3, "example_6")
# next, ensure each point is embedded within 50 m of contiguous habitat, and at least 100 m from each other point
# example_4 <- generateAllPointsForParameterSet(100, 100, 30)
# printPointSummary(example_4, "example_4")
# # next, ensure each point is embedded within 50 m of contiguous habitat, and at least 100 m from each other point
# example_5 <- generateAllPointsForParameterSet(150, 150, 30)
# printPointSummary(example_5, "example_5")


