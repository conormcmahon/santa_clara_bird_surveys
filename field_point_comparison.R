
library(raster)
library(sf)
library(tidyverse)
library(here)

# Open Questions
#   is there overlap between cover across species? total can be > 100
#   are there dates for the points? 

# Load data
saticoy_2022 <- st_read(here::here("..","..","vegetation_surveys","Veg_Points_Saticoy_2022","SurveyPoints.shp"))
# Add overall cover (Arundo + Percent1 + Percent2 + Percent3)
#   note - this may underestimate cover in very complex stands with no dominant vegetation species
#   also, may overestimate vegetation if some plants are overlying each other 
saticoy_2022$overall_cover <- (saticoy_2022$ArundoPerc+saticoy_2022$Percent1+saticoy_2022$Percent2+saticoy_2022$Percent3)

# Print a species list
sort(unique(c(saticoy_2022$Species1, saticoy_2022$Species2, saticoy_2022$Species3)))
woodland_species <- c("Euc glo", "Euc spp", "Jug cal", "Pla rac", "Plat rac", "Pop fre", "Pop spp", "Plat syc", "Pop tri", "Que agr", "Sal lae", "Sal las", "Salix spp", "Tam ram", "Tam spp")
drw_species <- c("Jug cal", "Pla rac", "Plat rac", "Pop fre", "Pop spp", "Plat syc", "Pop tri", "Sal lae", "Sal las", "Salix spp", "Tam ram", "Tam spp")
shrubby_riparian <- c("Bac Sal", "Bac sal", "Sal exi", "Salix exi")
shrub_species <- c("Bac pil", "Het arb", "Sam nig")
obligate_riparian_woody <- c(drw_species, shrubby_riparian)

# Get cover of various guilds 
saticoy_2022$woodland_cover <- (saticoy_2022$Species1 %in% woodland_species) * saticoy_2022$Percent1 + 
  (saticoy_2022$Species2 %in% woodland_species) * saticoy_2022$Percent2 + 
  (saticoy_2022$Species3 %in% woodland_species) * saticoy_2022$Percent3
saticoy_2022$drw_species <- (saticoy_2022$Species1 %in% drw_species) * saticoy_2022$Percent1 + 
  (saticoy_2022$Species2 %in% drw_species) * saticoy_2022$Percent2 + 
  (saticoy_2022$Species3 %in% drw_species) * saticoy_2022$Percent3
saticoy_2022$shrubby_riparian <- (saticoy_2022$Species1 %in% shrubby_riparian) * saticoy_2022$Percent1 + 
  (saticoy_2022$Species2 %in% shrubby_riparian) * saticoy_2022$Percent2 + 
  (saticoy_2022$Species3 %in% shrubby_riparian) * saticoy_2022$Percent3
saticoy_2022$obligate_riparian_woody <- (saticoy_2022$Species1 %in% obligate_riparian_woody) * saticoy_2022$Percent1 + 
  (saticoy_2022$Species2 %in% obligate_riparian_woody) * saticoy_2022$Percent2 + 
  (saticoy_2022$Species3 %in% obligate_riparian_woody) * saticoy_2022$Percent3

# Visualize variation in overall cover
ggplot(saticoy_2022) + 
  geom_histogram(aes(x=overall_cover), alpha=0.5) + 
  ggtitle("Overall Vegetation Cover")
# Visualize variation in Arundo cover
ggplot(saticoy_2022) + 
  geom_histogram(aes(x=ArundoPerc), alpha=0.5) + 
  ggtitle("Overall Vegetation Cover")
# Visualize variation in Woodland cover
ggplot(saticoy_2022) + 
  geom_histogram(aes(x=woodland_cover), alpha=0.5) + 
  ggtitle("Woodland Cover")
# Visualize variation in Woody Riparian Obligates
ggplot(saticoy_2022) + 
  geom_histogram(aes(x=obligate_riparian_woody), alpha=0.5) + 
  ggtitle("Cover of Woody Obligately Riparian Species")

# Output modified data 
st_write(saticoy_2022, here::here("..","..","vegetation_surveys","Veg_Points_Saticoy_2022","SurveyPointsModified.gpkg"))