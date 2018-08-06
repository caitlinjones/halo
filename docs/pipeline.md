### final_pipeline.R

#### Main pipeline input is a single configuration file containing all parameters needed for all analyses to be run. 

#### Step 1:
* Study configuration and study meta data is read and parsed.
* If __annotations_file__ is not specified in study config, all files in __annotations_dirs__ will be parsed and stored.
* If __--markExclusions__ is used and __rawDataDir__ or __rawDataFiles__ is provided, an EXCLUDE column will be added to data indicating which cells to ignore during analysis. A cell is determined to be excluded based on information in:
    *  __\*FOVannotations.xlsx__ in meta data directory
    *  files in __driftDir__ or __driftFiles__
    *  Halo annotation boundary files in __annotationsDirs__ or __annotationsFile__
  Output files will be written to __dataDir__ and files will be named with '_Excl' appended to the associated raw data file name.
 
