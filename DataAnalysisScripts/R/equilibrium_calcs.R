## Atmospheric equilibrium calculations
## Author: Abby Lewis

## First, load data
# Atmos CO2
atmos_co2 <- read.delim("https://gml.noaa.gov/webdata/ccgg/trends/co2/co2_mm_mlo.txt",
  comment.char = "#", header = F, sep = "",
  col.names = c(
    "Year", "Month", "Time", "CO2_Concentration",
    "Interpolated", "Trend", "Days", "Uncertainty"
  ),
  na.strings = c("-99.99")
) %>%
  group_by(Year, Month) %>%
  summarize(mean = mean(Interpolated, na.rm = T))

# Atmos CH4
atmos_ch4 <- read.delim("https://gml.noaa.gov/webdata/ccgg/trends/ch4/ch4_mm_gl.txt",
  comment.char = "#", header = F, sep = "",
  col.names = c(
    "Year", "Month", "Decimal_year", "Average",
    "Average_unc", "Trend", "Trend_unc"
  ),
  na.strings = c("-9.99")
) %>%
  group_by(Year, Month) %>%
  summarize(mean = mean(Average, na.rm = T))

# Function to calculate CO2 equilibrium concentration in µmol/L
# From Weiss 1974
calc_CO2_equil <- function(temp_C, year, month) {
  CO2_ppm <- atmos_co2$mean[match(
    paste(year, month),
    paste(atmos_co2$Year, atmos_co2$Month)
  )]
  conc_umol_L <- marelac::gas_satconc(S = 0, t = temp_C, species = "CO2", atm = CO2_ppm / 1000 / 1000)
  return(conc_umol_L)
}

# Function to calculate CH4 equilibrium concentration in µmol/L
# Based on Yamamoto et al. (1976)
calc_CH4_equil <- function(temp_C, year, month) {
  CH4_ppb <- atmos_ch4$mean[match(
    paste(year, month),
    paste(atmos_ch4$Year, atmos_ch4$Month)
  )]
  conc_umol_L <- marelac::gas_satconc(S = 0, t = temp_C, species = "CH4", atm = CH4_ppb / 1000 / 1000 / 1000)
  return(conc_umol_L)
}

# Function to calculate DO equilibrium concentration in µmol/L
calc_DO_equil <- function(temp_C) {
  conc_umol_L <- marelac::gas_satconc(S = 0, t = temp_C, species = "O2")
  return(conc_umol_L)
}
