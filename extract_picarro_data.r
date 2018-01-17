# Program: extract_picarro_data.r
# Version 0.6
# Author Dominic Woolf
# Revised Jan 17 2018, 


#############################################################################################################
#  Load Packages 
#############################################################################################################
ensure_library = function (lib.name){
  # Tries to load package, installing from CRAN if not already installed.
  # Returns logical value for whether successful
  x = require(lib.name, quietly = TRUE, character.only = TRUE)
  if (!x) {
    install.packages(lib.name, dependencies = TRUE, quiet = TRUE)
    x = require(lib.name, quietly = TRUE, character.only = TRUE)
  }
  x
}
invisible({
ensure_library('data.table')
ensure_library('ggplot2')
ensure_library('cowplot')
ensure_library('tcltk')
})

if (!require(regSmooth, quietly=T)) { # grab package regSmooth from github if it is not installed already
  library(devtools)
  install_github('domwoolf/regSmooth')
  library(regSmooth) # provides advanced data smoothing by Tikhonov regularization
}

#############################################################################################################
#  Some helper functions to select a directory interactively in a platform-independednt and fail-safe way 
#############################################################################################################
select_directory_method = function() {
  # Tries out a sequence of potential methods for selecting a directory to find one that works 
  # The fallback default method if nothing else works is to get user input from the console
  if (exists('utils::choose.dir')) {
    .dir.method = 'choose.dir'
  } else if (rstudioapi::isAvailable() & rstudioapi::getVersion() > '1.1.287') {
    .dir.method = 'RStudioAPI'
    ensure_library('rstudioapi')
  } else if(ensure_library('tcltk') & 
            class(try({tt  <- tktoplevel(); tkdestroy(tt)}, silent = TRUE)) != "try-error") {
    .dir.method = 'tcltk'
  } else if (ensure_library('gWidgets2') & ensure_library('RGtk2')) {
    .dir.method = 'gWidgets2RGtk2'
  } else if (ensure_library('rJava') & ensure_library('rChoiceDialogs')) {
    .dir.method = 'rChoiceDialogs'
  } else {
    .dir.method = 'console'
  }
  .dir.method
}

choose_directory = function(method = .dir.method, title = 'Select data directory') {
  # Interactively select a directory using predetermined method
  # If no method has already been established, calls select_directory_method to find one that works 
  if (missing(method)) {
    if (!exists('.dir.method')) assign('.dir.method', select_directory_method(), envir = .GlobalEnv)
    method = .dir.method
  }
  switch (method,
          'choose.dir' = choose.dir(caption = title),
          'RStudioAPI' = selectDirectory(caption = title),
          'tcltk' = tk_choose.dir(caption = title),
          'rChoiceDialogs' = rchoose.dir(caption = title),
          'gWidgets2RGtk2' = gfile(type = 'selectdir', text = title),
          readline('Please enter directory path: ')
  )
}


#############################################################################################################
# Data extraction begins here
#############################################################################################################
extract_picarro = function(data_path = NA, lambda = 1e-4) {
  # This function returns a data.table with respired CO2 concentration and d13C for each sample at each samping time
  #
  # Required inputs:
  #   data_path = directory containing all raw data from picarro analyzer within nested subdirectories
  #               if data_path is not provided in the function call, it can be selected interactively
  #   NB picarro files are assumed to end in .dat -- No other files in the data path should have this ending
  #   In addition to raw Picarro data, the data_path must also contain the following 2 files:
  #   fractional_volume.txt  = provides fraction of head space gas divided by total volume including residal gas in the analyzer and pipework
  #   logfile (filename must start with the characters "logfile"), which contains epoch time, sample ID, and step ID:
  #           Step 2 = sample analysis
  #           Step 3 = end of step 2
  #           Step 5 = purge analysis
  #           Step 6 = end of step 5
  #
  # Optional parameters:
  #   lambda is the smoothing parameter - see help(regSmooth)
  #
  #############################################################################################################
  #  Load Data Files 
  #############################################################################################################
  if (is.na(data_path)) data_path = choose_directory(title = "Select data directory")
  logfile = list.files(data_path, pattern = 'logfile*', full.names = TRUE)
  if (length(logfile) != 1L) stop('Need exactly one logfile in the data path')
  log_data = fread(logfile)
  fractional_volume = fread(paste0(data_path,"/fractional_volume.txt"))
  
  picarro_files = list.files(data_path, pattern='.*[.]dat$', recursive = T, full.names = TRUE)
  if (!(length(picarro_files) >= 1)) stop('No Picarro data files found in selected directory')
  columns_to_keep = c("EPOCH_TIME", "HP_12CH4_dry", "HP_Delta_iCH4_Raw", "12CO2_dry", "Delta_Raw_iCO2")
  fread.picarro = function(fname) {
    pic.data.i = fread(fname)
    if (length(setdiff(columns_to_keep, names(pic.data.i)))) {
      stop(paste('Columns', setdiff(columns_to_keep, names(pic.data.i)), 'missing from', fname))
    }
    pic.data.i
  }
  pic.data = fread.picarro(picarro_files[1]) # read first picarro data file
  for(i in 2:length(picarro_files)) { # then append all the other files into the same data.table
    pic.data <- rbind(pic.data, fread.picarro(picarro_files[i]), fill=TRUE)
  }
  setkey(pic.data, 'EPOCH_TIME') # sort data chronologically
  
  #############################################################################################################
  #  Data Preprocessing  
  #############################################################################################################
  drop.cols = setdiff(names(pic.data), columns_to_keep)
  pic.data[, (drop.cols) := NULL]   # drop unused columns
  pic.data = pic.data[EPOCH_TIME > log_data$epoch[1]] # discard any picarro data from before the first entry in the log file
  pic.data[, log_previous := findInterval(EPOCH_TIME, log_data$epoch, left.open  = TRUE)] # find which line of the logfile is the last one before current data epoch time
  pic.data[, step   := log_data[log_previous, step]]   # read (from logfile) which step we are on at current time
  pic.data[, sample := log_data[log_previous, sample]] # read (from logfile) which sample we are on at current time
  
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
  cat('Smoothing data. This may take a minute.\n')
  pic.data[, CO2_12_smooth := regSmooth(EPOCH_TIME, `12CO2_dry`, lambda = 1e-4)$yhat, by = .(combined_sample_cycle, step)]
  cat('Finished smoothing :-)')
  #############################################################################################################
  #  Reduce data to one row per sample measurement  
  #############################################################################################################
  long.data = pic.data[, .(epoch = (EPOCH_TIME[1] + EPOCH_TIME[.N])/2,  # Use middle of sampling step as representative time
                           CO2 = max(CO2_12_smooth, na.rm = T), # Use maximum CO2 concentration for respired CO2 in the headspace
                           d13C = Delta_Raw_iCO2[which.max(CO2_12_smooth)], # d13C at same time as maxCO2
                           last.co2 = CO2_12_smooth[.N],  # Use final CO2 concentration for purge CO2 in the headspace
                           last.d13C = Delta_Raw_iCO2[.N],
                           dilute.co2 = CO2_12_smooth[1],  # dilute.co2 is 1st value from each sample/cycle (i.e. concentration in pipework before headspace gas reaches the analyser)
                           dilute.d13C = Delta_Raw_iCO2[1],
                           CH4 = HP_12CH4_dry[which.max(CO2_12_smooth)], # CH4 at same time as maxCO2
                           d13CH4 = HP_Delta_iCH4_Raw[which.max(CO2_12_smooth)], # d13C at same time as maxCO2
                           last.ch4 = HP_12CH4_dry[.N],  # Use final CH4 concentration for purge CH4                           
                           last.d13CH4 = HP_Delta_iCH4_Raw[.N],
                           dilute.ch4 = HP_12CH4_dry[1],  # dilute.ch4 is 1st value from each sample/cycle (i.e. concentration in pipework before headspace gas reaches the analyser)
                           dilute.d13CH4 = HP_Delta_iCH4_Raw[1]
  ), by=.(sample, cycle, combined_sample_cycle, step)] 
  long.data[step == 'purge',  `:=`(CO2 = last.co2, 
                                   d13C = last.d13C,
                                   CH4 = last.ch4,
                                   d13CH4 = last.d13CH4)]
  
  # convert from long to wide format with separate columns for purge and respiration data
  short.data = dcast(long.data, 
                     combined_sample_cycle ~ step, 
                     value.var = c('epoch','sample','cycle',
                                   'CO2', 'd13C', 'dilute.co2', 'dilute.d13C',
                                   "CH4", "d13CH4", "dilute.ch4", "dilute.d13CH4"))
  short.data = short.data[complete.cases(short.data)] # There's no purge data to go with the first respiration, nor respiration data to go with the final purge - so we remove these lines
  short.data[, (c('sample_purge', 'cycle_purge', 'epoch_purge')) := NULL] # the dcast gave us duplicates of the sample and cycle columns, so lets get rid of the extras
  setnames(short.data, c('epoch_respiration', 'sample_respiration', 'cycle_respiration'), c('epoch','sample', 'cycle'))
  
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

  #############################################################################################################
  #  Apply the same correction factors to the CH4 data
  #############################################################################################################
  short.data[, headspace.ch4  := (CH4_respiration - (1-fractional_volume) * dilute.ch4_respiration) / fractional_volume ]
  short.data[, headspace.ch4_purge  := (CH4_purge - (1-fractional_volume) * dilute.ch4_purge) / fractional_volume ]
  short.data[, headspace.d13ch4 := (d13CH4_respiration - (1-fractional_mass)*dilute.d13CH4_respiration) / fractional_mass]
  short.data[, evolved.ch4 := headspace.ch4 - headspace.ch4_purge] # subtract ppm co2 that was in the jar to begin with

  #############################################################################################################
  #  Add a days since start column, and we're done!
  #############################################################################################################
  short.data[, day := (epoch - epoch[1]) / (60 * 60 * 24), by = sample]
  return(short.data)
}


#############################################################################################################
# Examples and vignettes
#############################################################################################################
short.data = extract_picarro()
ggplot(short.data[-1], aes(day, respired.co2)) +
  geom_line() +
  facet_wrap(~ sample)

ggplot(short.data[-1], aes(day, evolved.ch4)) +
  geom_line() +
  facet_wrap(~ sample)
