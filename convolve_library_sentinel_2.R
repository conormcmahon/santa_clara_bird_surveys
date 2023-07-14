
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

s2_band_centers <- c(492.4, 559.8, 664.6, 704.1, 740.5, 782.8, 832.8, 864.7, 1613.7, 2202.4)
s2_fwhm <- c(66, 36, 31, 15, 15, 20, 106, 21, 91, 175)

getResponseAtBand <- function(band)
{
  within_lower <- wavelengths_in > s2_band_centers[[band]] - s2_fwhm[[band]]/2
  within_higher <- wavelengths_in < s2_band_centers[[band]] + s2_fwhm[[band]]/2
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
write_csv(metadata_in, here::here("..","..","..","sentinel_2","mesma","west_usa_s2.csv"))






