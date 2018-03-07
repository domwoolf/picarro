# picarro
Contains functions for processing and analyzing picarro respiration data from Lehmann lab parallel incubation rig

Installation instructions in R:

  1. library(devtools)
  2. install_github('/domwoolf/picarro')

Required inputs:

  data_path = directory containing all raw data from picarro analyzer within nested subdirectories.
    If data_path is not provided in the function call, it can be selected interactively.
    **NB picarro files are assumed to end in .dat -- No other files in the data path should have this ending**.
    
  In addition to raw Picarro data, the data_path must also contain the following 2 files:
  
  1. fractional_volume.txt  = provides fraction of head space gas divided by total volume including residal gas in the analyzer and pipework
  2. logfile (filename must start with the characters "logfile"), which contains epoch time, sample ID, and step ID:
     Step 2 = sample analysis,
     Step 3 = end of step 2,
     Step 5 = purge analysis,
     Step 6 = end of step 5

