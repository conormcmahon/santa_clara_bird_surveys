
library(here)
library(RStoolbox)
library(tidyverse)

# Load input library as a dataframe
#   first column is wavelength
#   each subsequent column is a spectrum
input_library <- readSLI(here::here("..","..","..","planet","mesma","010627_westusa_all.sli"))
# Separate wavelength vector
wavelengths_in <- input_library$wavelength
# Reformat data into a matrix of spectra (without wavelength)
input_library_mat <- as.matrix(input_library[,2:ncol(input_library)])

planet_band_centers <- c(443, 490, 531, 565, 610, 665, 705, 865)
planet_fwhm <- c(20, 50, 36, 36, 20, 31, 15, 40)

getResponseAtBand <- function(band)
{
  within_lower <- wavelengths_in > planet_band_centers[[band]] - planet_fwhm[[band]]/2
  within_higher <- wavelengths_in < planet_band_centers[[band]] + planet_fwhm[[band]]/2
  bands_within_range <- within_lower * within_higher
  
  colMeans(input_library_mat[which(as.logical(bands_within_range)),])
}

# Get mean spectral values within each target band
output_df_no_wavelength <- bind_rows(lapply(1:length(planet_band_centers), getResponseAtBand))
# Apply a scale factor shift to make library match intended imagery 
scale_factor_shift <- 10
output_df_no_wavelength <- output_df_no_wavelength * scale_factor_shift
# Add wavelength back into library
output_library <- data.frame(wavelength = planet_band_centers)
output_library <- cbind(output_library, output_df_no_wavelength)
writeSLI(output_library, 
         here::here("..","..","..","planet","mesma","west_usa_planet.sli"),
         wavl.units = "Nanometers")
# Copy over metadata
metadata_in <- read_csv(here::here("..","..","..","planet","mesma","010627_westusa_all.csv"))
metadata_in$Name = names(output_library[,2:ncol(output_library)])
write_csv(metadata_in, here::here("..","..","..","planet","mesma","west_usa_planet.csv"))






