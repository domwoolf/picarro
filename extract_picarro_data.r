# Program:
# extract_picarro.r
#
# Version 0.4
# Revised Jan 5 2018, 
# Author Dominic Woolf
#

#############################################################################################################
#  Load Packages 
#############################################################################################################
library(data.table)
library(ggplot2)
library(cowplot)
library(tcltk)
if (!require(regSmooth, quietly=T)) { # grab package regSmooth from github if it is not installed already
  library(devtools)
  install_github('domwoolf/regSmooth')
  library(regSmooth) # provides advanced data smoothing by Tikhonov regularization
}


extract_picarro = function(data_path = NA, lambda = 1e-4) {
  # This function returns a data.table with respired CO2 concentration and d13C for each sample at each samping time
  #
  # Required inputs:
  #   data_path = directory containing all raw data from picarro analyzer within nested subdirectories
  #               if data_path is not provided in the function call, it can be selected interactively
  #   In addition to raw Picarro data, the data_path must also contain the following 2 files:
  #   fractional_volume.txt  = provides fraction of head space gas divided by total volume including residal gas in the analyzer and pipework
  #   logfile (filename must start with the characters "logfile"), which contains epoch time, sample ID, and step ID:
  #           Step 2 = sample analysis
  #           Step 3 = end of step 2
  #           Step 5 = purge analysis
  #           Step 6 = end of step 5
  #
  #############################################################################################################
  #  Load Data Files 
  #############################################################################################################
  if (is.na(data_path)) data_path = tk_choose.dir(default = "", caption = "Select data directory")
  logfile = list.files(data_path, pattern = 'logfile*', full.names = TRUE)
  if (length(logfile) != 1L) stop('Need exactly one logfile in the data path')
  log_data = fread(logfile)
  picarro_files = list.files(data_path, pattern='.*[.]dat', recursive = T, full.names = TRUE)
  picarro_dates = strsplit( picarro_files, '-') # extract date from file names to sort chronologically
  picarro_dates = sapply(picarro_dates, `[`, 2) # date is the 2nd part of the file name
  picarro_files = picarro_files[order(picarro_dates)] # sort files chronologically
  pic.data = fread(picarro_files[1]) # read first picarro data file
  for(i in 2:length(picarro_files)) { # then append all the other files into the same data.table
    pic.data <- rbind(pic.data, fread(picarro_files[i]), fill=TRUE) 
  }
  pic.data = pic.data[EPOCH_TIME > log_data$epoch[1]] # discard any picarro data from before the first entry in the log file
  pic.data[, log_previous := findInterval(EPOCH_TIME, log_data$epoch, left.open  = TRUE)] # find which line of the logfile is the last one before current data epoch time
  pic.data[, step := log_data[log_previous, step]]   # read (from logfile) which step we are on at current time
  pic.data[, sample := log_data[log_previous, sample]] # read (from logfile) which sample we are on at current time
  fractional_volume  <- fread(paste0(data_path,"/fractional_volume.txt"))
  
  #############################################################################################################
  #  Data Preprocessing  
  #############################################################################################################
  columns_to_keep = c("EPOCH_TIME", "HP_12CH4_dry", "HP_Delta_iCH4_Raw", "12CO2_dry", "Delta_Raw_iCO2", 'log_previous',
                      'step', 'sample')
  drop.cols = setdiff(names(pic.data), columns_to_keep)
  pic.data[, (drop.cols) := NULL]   # drop unused columns
  pic.data = pic.data[step %in% c('step2', 'step5')] # only keep sample and purge analysis data
  pic.data = pic.data[complete.cases(pic.data)] # remove any rows with NA
  pic.data[step == 'step2', step := 'respiration'] # give the steps meaningful names
  pic.data[step == 'step5', step := 'purge']
  cycle.length = which(log_data[-1, paste(step,sample)] == log_data[1, paste(step,sample)])[1] # how many lines in log file before we return to the same sample and step
  
  log_data[, cycle := ceiling((1:.N) / cycle.length)] # number the data collection cycles. Each cycle represents how long it takes to go through all samples and the start over
  pic.data[, cycle := log_data[log_previous, cycle]]
  
  pic.data[step=="purge", cycle := cycle + 1L]     # Associate purge data with analysis data that follows (rather than precedes) it
  pic.data[, combined_sample_cycle := factor(paste(formatC(sample, digits = 4, flag="0"), formatC(cycle, digits = 4, flag="0"), sep=","))]
  # smooth the CO2 data before finding maximum... to protect against noise spikes
  pic.data[, CO2_12_smooth := regSmooth(EPOCH_TIME, `12CO2_dry`, lambda = 1e-4)$yhat, by = .(combined_sample_cycle, step)]
  
  #############################################################################################################
  #  Reduce data to one row per sample measurement  
  #############################################################################################################
  long.data = pic.data[, .(CO2 = max(CO2_12_smooth, na.rm = T),# Use maximum CO2 concentration for respired CO2 in the headspace
                           d13C = Delta_Raw_iCO2[which.max(CO2_12_smooth)], # d13C at same time as maxCO2
                           last.co2 = CO2_12_smooth[.N], # Use final CO2 concentration for purge CO2 in the headspace
                           last.d13C = Delta_Raw_iCO2[.N],
                           dilute.co2 = CO2_12_smooth[1],  # dilute.co2 is 1st value from each sample/cycle (i.e. concentration in pipework before headspace gas reaches the analyser)
                           dilute.d13C = Delta_Raw_iCO2[1]
  ), by=.(sample, cycle, combined_sample_cycle, step)] 
  long.data[step == 'purge',  `:=`(CO2 = last.co2, d13C = last.d13C) ]
  
  # convert from long to wide format with separate columns for purge and respiration data
  short.data = dcast(long.data, combined_sample_cycle ~ step, 
                     value.var = c('sample','cycle','CO2', 'd13C', 'dilute.co2', 'dilute.d13C'))
  short.data = short.data[complete.cases(short.data)] # There's no purge data to go with the first respiration, nor respiration data to go with the final purge - so we remove these lines
  short.data[, (c('sample_purge', 'cycle_purge')) := NULL] # the dcast gave us duplicates of the sample and cycle columns, so lets get rid of the extras
  setnames(short.data, c('sample_respiration', 'cycle_respiration'), c('sample', 'cycle'))
  
  #############################################################################################################
  #  Correct CO2 concentrations for residual co2 in pipework & analyzer
  #############################################################################################################
  # fractional_volume from the file fractional_volume.txt... 
  # ...provides headspace divided by total volume (incl. analyzer and pipework)
  short.data = merge(short.data, fractional_volume, all.x = T, by='sample')
  
  # headspace CO2 is measured CO2, adjusted for the fraction of CO2 from pipework dilutant
  # NB note there were missing parentheses in previous version of code!
  short.data[, headspace.co2  := (CO2_respiration - (1-fractional_volume) * dilute.co2_respiration) / fractional_volume ]
  short.data[, headspace.co2_purge  := (CO2_purge - (1-fractional_volume) * dilute.co2_purge) / fractional_volume ]
  
  # Need to adjust d13C by quantity of dilutant co2 relative to quantity respired (not concentrations!!)
  short.data[, fractional_mass  := fractional_volume * CO2_respiration/headspace.co2] #like fractional volume, but on a mass of carbon basis
  short.data[, headspace.d13c := (d13C_respiration - (1-fractional_mass)*dilute.d13C_respiration)/fractional_mass]
  
  short.data[, respired.co2 := headspace.co2 - headspace.co2_purge] # subtract ppm co2 that was in the jar to begin with
  return(short.data)
}

# Example
short.data = extract_picarro()
ggplot(short.data[-1], aes(cycle, respired.co2)) +
  geom_line() +
  facet_wrap(~ sample)
