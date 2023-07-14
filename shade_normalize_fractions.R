
library(raster)
library(tidyverse)

directory <- here::here("..","..","..","planet","2023_floods")
mesma_fraction_files <- c("santa_clarita/2023_may_composite.hdr_SMA_20230531T08H56M16S",
                          "santa_clarita/2022_may_composite.hdr_SMA_20230531T08H59M02S",
                          "santa_clarita/santa_clarita_2022_12_20/files/santa_clarita_2022_12_20.hdr_SMA_20230531T09H18M56S",
                          "santa_clarita/santa_clarita_2023_02_19/files/santa_clarita_2023_02_19.hdr_SMA_20230531T09H16M48S",
                          "santa_clarita/santa_claria_2023_03_07/files/santa_clarita_2023_03_07.hdr_SMA_20230531T09H20M14S",
                          "saticoy/saticoy_2022_12_23/files/saticoy_2022_12_23.hdr_SMA_20230531T10H46M18S",
                          "saticoy/saticoy_2023_02_07/files/saticoy_2023_02_07.hdr_SMA_20230531T10H51M09S",
                          "saticoy/saticoy_2023_03_07_psscene_analytic_8b_sr_udm2/files/saticoy_2023_03_07.hdr_SMA_20230531T10H56M11S")

normalizeFractions <- function(filename)
{
  print(paste("Working on file ", paste(directory, "/", filename, sep="")))
  input_img <- stack(paste(directory, "/", filename, sep=""))
  
  normalized_img <- input_img[[1:3]]
  
  shade_frac <- input_img[[4]]
  non_shade_frac <- 1-shade_frac
  
  normalized_img <- normalized_img / non_shade_frac
  
  writeRaster(normalized_img, paste(directory, "/", filename, "_shade_normalized.tif", sep=""), overwrite=TRUE)
}

lapply(mesma_fraction_files, normalizeFractions)