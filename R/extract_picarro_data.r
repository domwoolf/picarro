# Program: extract_picarro_data.r
# Author Dominic Woolf
# Revised Jan 19 2018, 
# Contains functions to extract raw picarro data and convert to respired CO2 & CH4 

co2.density_30C = 0.001778
co2.carbon.frac = 12.01 / 44.0095

#' Interactive dialogue to allow user to select a directory
#' 
#' Function to select a directory interactively in a platform-independent way
#' @details
#' This function allows directory selection in a way that should work 
#' on any (most?) platforms.
#' @param caption Optional text to print in title of selection dialogue
#' @export
#' @return Length one character vector containing path to selected directory
#' @author
#'   Dominic Woolf.
#'   d.woolf@cornell.edu
#' @examples
#' choose_directory()
choose_directory = function(caption = 'Select data directory') {
  if (exists('utils::choose.dir')) {
    choose.dir(caption = caption) 
    } else {
      tk_choose.dir(caption = caption)
    }
  }


#' Data extraction from Lehmann lab Picarro isotope analyzer
#'
#' Function to convert raw data from picarro analyzer to respiration data
#' @details
#' This function returns a data.table with respired CO2 and CH4 
#' concentrations and delta-13C values for each jar at each samping time
#'
#' The data path \strong{must} contain the following three things: 
#' \enumerate{
#'   \item All raw Picarro data files for this experiment. Can be in nested subdirectories.
#'   NB picarro files are assumed to end in .dat -- No other files in the data path should have this ending
#'   \item fractional_volume.txt  = provides fraction of head space gas divided by total volume including residal gas in the analyzer and pipework
#'   \item logfile (filename must start with the characters "logfile"), which contains epoch time, sample ID, and step ID 
#'           ("Step 2" = sample analysis;
#'            "Step 3" = end of step 2;
#'            "Step 5" = purge analysis;
#'            "Step 6" = end of step 5).
#' }
#' The returned data.table columns, and their descriptions are:
#' 
#' jar = jar number to which row corresponds
#' epoch = epoch time of measurement
#' cycle = number of sampling cycle (i.e. increments each time a specific jar is analyzed)
#' combined_jar_cycle = jar number and cycle number pasted together
#' CO2_purge = Picarro-measured CO2 concentration at the end of the purge step
#' CO2_respiration = Picarro-measured CO2 concentration at the peak (smoothed maximum) of the sampling step 
#' d13C_purge = Picarro-measured d13-CO2 at the end of the purge step
#' d13C_respiration = Picarro-measured d13-CO2 at the concentration-peak (smoothed maximum) of the sampling step 
#' <dilute.co2 is 1st value from each jar/cycle. I.e. concentration in pipework before headspace gas reaches the analyser>
#' dilute.co2_purge = dilute.co2 concentration in purge step
#' dilute.co2_respiration = dilute.co2 concentration in sampling step
#' dilute.d13C_purge = dilute.co2 d13C in purge step
#' dilute.d13C_respiration = dilute.co2 d13C in sampling step
#' fractional_volume = jar-specific correction factor for fraction of gas at picarro that originates from within the headspace rather than the pipework
#' headspace.co2 = CO2_respiration, corrected for fractional_volume
#' headspace.co2_purge = CO2_purge, corrected for fractional_volume
#' fractional_mass = fractional_volume converted to mass of CO2 basis
#' respired.co2 = CO2_respiration, corrected for fractional_volume
#' headspace.d13c = d13C_respiration, corrected for fractional_volume \strong{This correction provides highly suspect values.  Use d13C_respiration values instead}
#' day = epoch time converted into days since the first measurement on that sample

#' <Methane (CH4 has columns corresponding to the CO2 descriptions above):
#'    CH4_purge, CH4_respiration, d13CH4_purge, d13CH4_respiration, dilute.ch4_purge, 
#'    dilute.ch4_respiration, dilute.d13CH4_purge, dilute.d13CH4_respiration
#'    CH4_purge, headspace.ch4, headspace.ch4_purge, headspace.d13ch4, evolved.ch4
#'
#' @param data_path Directory containing all raw data from picarro analyzer within nested subdirectories. If not provided, user will be asked for path interactively.
#' @param lambda Optional smoothing parameter (positive real numeric). Higher values give greater smoothing. Smaller values follow data more closely. See help(regSmooth) for more details.
#' @import data.table 
#' @export
#' @return Provides a data.table with respired CO2 and CH4 concentrations and delta-13C values for each jar at each samping time
#' @author
#'   Dominic Woolf.
#'   d.woolf@cornell.edu
#' @examples
#' short.data = extract_picarro()
#' library(ggplot2)
#' ggplot(short.data[-1], aes(day, respired.co2)) +
#'   geom_line() +
#'   facet_wrap(~ jar)
#' 
#' ggplot(short.data[-1], aes(day, evolved.ch4)) +
#'   geom_line() +
#'   facet_wrap(~ jar)
extract_picarro = function(data_path = NA, lambda = 1e-4) {
  #
  #############################################################################################################
  #  Load Data Files 
  #############################################################################################################
  if (is.na(data_path)) data_path = choose_directory(caption = "Select data directory")
  logfile = list.files(data_path, pattern = 'logfile*', full.names = TRUE)
  if (length(logfile) != 1L) stop('Need exactly one logfile in the data path')
  log_data = fread(logfile)
  fractional_volume = fread(paste0(data_path,"/fractional_volume.txt"))
  
  picarro_files = list.files(data_path, pattern='.*[.]dat$', recursive = T, full.names = TRUE)
  if (!(length(picarro_files) >= 1)) stop('No Picarro data files found in selected directory')
  cat('I will be very happy to extract the data for you. Please be patient, I\'m working as fast as I can.\n')
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
  pic.data[, jar := log_data[log_previous, sample]] # read (from logfile) which sample we are on at current time
  
  pic.data = pic.data[step %in% c('step2', 'step5')] # only keep sample and purge analysis data
  pic.data = pic.data[complete.cases(pic.data)] # remove any rows with NA
  pic.data[step == 'step2', step := 'respiration'] # give the steps meaningful names
  pic.data[step == 'step5', step := 'purge']
  cycle.length = which(log_data[-1, paste(step,sample)] == log_data[1, paste(step,sample)])[1] # how many lines in log file before we return to the same jar and step
  
  log_data[, cycle := ceiling((1:.N) / cycle.length)] # number the data collection cycles. Each cycle represents how long it takes to go through all jars and the start over
  pic.data[, cycle := log_data[log_previous, cycle]]
  
  pic.data[step=="purge", cycle := cycle + 1L]     # Associate purge data with analysis data that follows (rather than precedes) it
  pic.data[, combined_jar_cycle := factor(paste(formatC(jar, digits = 4, flag="0"), formatC(cycle, digits = 4, flag="0"), sep=","))]
  # smooth the CO2 data before finding maximum... to protect against noise spikes
  cat('Smoothing data. This may take a minute.\n')
  pic.data[, CO2_12_smooth := regSmooth(EPOCH_TIME, `12CO2_dry`, lambda = 1e-4)$yhat, by = .(combined_jar_cycle, step)]
  cat('Finished smoothing.\n')
  #############################################################################################################
  #  Reduce data to one row per jar measurement  
  #############################################################################################################
  long.data = pic.data[, .(epoch = (EPOCH_TIME[1] + EPOCH_TIME[.N])/2,  # Use middle of sampling step as representative time
                           CO2 = max(CO2_12_smooth, na.rm = T), # Use maximum CO2 concentration for respired CO2 in the headspace
                           d13C = Delta_Raw_iCO2[which.max(CO2_12_smooth)], # d13C at same time as maxCO2
                           last.co2 = CO2_12_smooth[.N],  # Use final CO2 concentration for purge CO2 in the headspace
                           last.d13C = Delta_Raw_iCO2[.N],
                           dilute.co2 = CO2_12_smooth[1],  # dilute.co2 is 1st value from each jar/cycle (i.e. concentration in pipework before headspace gas reaches the analyser)
                           dilute.d13C = Delta_Raw_iCO2[1],
                           CH4 = HP_12CH4_dry[which.max(CO2_12_smooth)], # CH4 at same time as maxCO2
                           d13CH4 = HP_Delta_iCH4_Raw[which.max(CO2_12_smooth)], # d13C at same time as maxCO2
                           last.ch4 = HP_12CH4_dry[.N],  # Use final CH4 concentration for purge CH4                           
                           last.d13CH4 = HP_Delta_iCH4_Raw[.N],
                           dilute.ch4 = HP_12CH4_dry[1],  # dilute.ch4 is 1st value from each jar/cycle (i.e. concentration in pipework before headspace gas reaches the analyser)
                           dilute.d13CH4 = HP_Delta_iCH4_Raw[1]
  ), by=.(jar, cycle, combined_jar_cycle, step)] 
  long.data[step == 'purge',  `:=`(CO2 = last.co2, 
                                   d13C = last.d13C,
                                   CH4 = last.ch4,
                                   d13CH4 = last.d13CH4)]
  
  # convert from long to wide format with separate columns for purge and respiration data
  short.data = dcast(long.data, 
                     combined_jar_cycle ~ step, 
                     value.var = c('epoch','jar','cycle',
                                   'CO2', 'd13C', 'dilute.co2', 'dilute.d13C',
                                   "CH4", "d13CH4", "dilute.ch4", "dilute.d13CH4"))
  short.data = short.data[complete.cases(short.data)] # There's no purge data to go with the first respiration, nor respiration data to go with the final purge - so we remove these lines
  short.data[, (c('jar_purge', 'cycle_purge', 'epoch_purge')) := NULL] # the dcast gave us duplicates of the sample and cycle columns, so lets get rid of the extras
  setnames(short.data, c('epoch_respiration', 'jar_respiration', 'cycle_respiration'), c('epoch','jar', 'cycle'))
  
  #############################################################################################################
  #  Correct CO2 concentrations for residual co2 in pipework & analyzer
  #############################################################################################################
  # fractional_volume from the file fractional_volume.txt provides headspace divided by total volume (incl. analyzer and pipework)
  short.data = merge(short.data, fractional_volume, all.x = T, by='jar')
  
  # headspace CO2 is measured CO2, adjusted for the fraction of CO2 from pipework dilutant
  short.data[, headspace.co2  := (CO2_respiration - (1-fractional_volume) * dilute.co2_respiration) / fractional_volume ]
  short.data[, headspace.co2_purge  := (CO2_purge - (1-fractional_volume) * dilute.co2_purge) / fractional_volume ]
  
  # Need to adjust d13C by *quantity* of dilutant co2 relative to *quantity* respired (not concentrations!!)
  short.data[, fractional_mass  := fractional_volume * CO2_respiration/headspace.co2] #like fractional volume, but on a mass of carbon basis
  short.data[, headspace.d13c := (d13C_respiration - (1-fractional_mass)*dilute.d13C_respiration)/fractional_mass]
  
  # subtract ppm co2 that was in the jar to begin with
  short.data[, respired.co2 := headspace.co2 - headspace.co2_purge] 

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
  short.data[, day := (epoch - epoch[1]) / (60 * 60 * 24), by = jar]
  cat('All done. Have a nice day :-)\n')
  return(short.data)
}


#' Isotope partitioning of CO2 and CH4 data from Picarro analyzer
#'
#' Function to partition respiration data from picarro analyzer by source, using d13C of two substrates
#' @details
#' This function returns a data.table with respired CO2 and CH4 
#' concentrations and delta-13C values for each jar at each samping time
#'
#' The data path \strong{must} contain the following: 
#' \enumerate{
#'   \item jars.txt = file specifying the treatment in each jar. Also has columns for quantity of each substrate, and the headspace volume of the jar.
#'   \item d13c.txt = file specifying d13C of two substrates for each treatment. Treatment names must correspond with those in jars.txt
#' }
#'
#' @param respiration_data data.table containing extracted respiration data from Picarro analyzer. Typically created by function extract_picarro
#' @param data_path Directory containing jars.txt and treatments.txt.  If not provided, user will be asked for path interactively.
#' @import data.table 
#' @export
#' @return Provides a data.table with the original respired CO2 and CH4 concentrations and delta-13C values 
#' for each jar at each samping time, from the input data.table. 
#' Columns are added to this data.table for mass of CO2 and CH4 respired during ech sampling period,
#' cumulative values of the same - both total and partitioned per substrate. 
#' @author
#'   Dominic Woolf.
#'   d.woolf@cornell.edu
#' @examples
#' short.data = extract_picarro()
#' partition.data = partition_picarro(short.data)
partition_picarro = function(respiration_data = short.data, data_path = NA, remove_first = TRUE) {
  respiration_data = copy(respiration_data)
  # read in data files
  if (is.na(data_path)) data_path = choose_directory(caption = "Select data directory")
  jarfile = list.files(data_path, pattern = '^jars.txt$', full.names = TRUE)
  if (length(jarfile) != 1L) stop('Need exactly one jars.txt file in the data path')
  jar_data = fread(jarfile)
  treatmentsfile = list.files(data_path, pattern = '^treatments.txt$', full.names = TRUE)
  if (length(treatmentsfile) != 1L) stop('Need exactly one treatments.txt file in the data path')
  treatments_data = fread(treatmentsfile)
  
  #first reading is always off (possibly due to entrained air in the soil)
  # so, we remove the first measurement for each jar (if remove_first is true)
  if (remove_first) respiration_data = respiration_data[, .SD[-1], by = jar, .SDcols= names(respiration_data)]
  
  # merge jar and treament metadata into the repiration data
  respiration_data = respiration_data[jar_data, on='jar']
  respiration_data = respiration_data[treatments_data, on='treatment']
  
  # convert CO2 concentrations to mass of carbon, based on headspace volume
  respiration_data[, respired.co2_C.mg := respired.co2 * headspace.vol * co2.density_30C * co2.carbon.frac]

  #partition the CO2-carbon
  respiration_data[, amendment.co2_C.mg := respired.co2_C.mg * (amendment.d13c - d13C_respiration)/(amendment.d13c - soil.d13c)]
  respiration_data[, soil.co2_C.mg := respired.co2_C.mg - amendment.co2_C.mg]
  respiration_data[, cum.co2_C.mg := cumsum(respired.co2_C.mg), by = jar]
  respiration_data[, cum.amendment.co2_C.mg := cumsum(amendment.co2_C.mg), by = jar]
  respiration_data[, cum.soil.co2_C.mg := cumsum(soil.co2_C.mg), by = jar]
  return(invisible(respiration_data))
}
