#!/home/byrne/R/R-3.4.3/bin/R

options( java.parameters = c("-Xss2560k", "-Xmx8g") )
options('python_cmd'='/opt/common/CentOS_6-dev/python/python-3.5.1/bin/python')

library(argparse)

###
# Parse user input
###
parser <- ArgumentParser()

parser$add_argument("-r", "--rootDir", type="character", default=NULL,
                    help="root directory for current project")

parser$add_argument("-n", "--studyName", type="character", default=NULL,
                    help="name of current study")

parser$add_argument("-m", "--manifestFile", type="character", default=NULL,
                    help="file containing all project parameters")

parser$add_argument("-rf", "--rawDataFiles", type="character", default=NULL,
                    help="comma-separated (no spaces) list of *.rda files WITHOUT exclusions; specify only if running --markExclusions")

parser$add_argument("-rd", "--rawDataDir", type="character", default=NULL,
                    help="directory containing raw *.rda files, WITHOUT exclusions; specify only if running --markExclusions; NOTE: if rawDataFiles are given as well, value passed here will be ignored.")

parser$add_argument("-df", "--dataFiles", type="character", default=NULL,
                    help="comma-separated (no spaces) list of *.rda files")

parser$add_argument("-dd", "--dataDir", type="character", default=NULL,
                    help="directory containing *.rda files to be analyzed; NOTE: if dataFiles are given as well, value passed here will be ignored.")    
 
parser$add_argument("-mf", "--metaFiles", type="character", default=NULL,
                    help="comma-separated (no spaces) list of meta data files in *.xlsx format; see documentation for formats") 

parser$add_argument("-md", "--metaDir", type="character", default=NULL,
                    help="directory containing all *.xlsx meta data files in formats described in documentation; NOTE: if metaFiles is also given, value passed here will be ignored")

parser$add_argument("--driftDir", type="character", default=NULL,
                    help="directory containing drift summaries for each sample")

parser$add_argument("-p", "--pad", type="double", default=20, 
                    help="number in pixels to trim from each FOV")

parser$add_argument("--driftThreshold", type="double", default=0.1, 
                    help="maximum fraction of drift allowed in order to keep a cell")         

parser$add_argument("--maxDistanceFromInterface", type="double", default=360,
                    help="total distance inside of and outside of tumor interface to include in infiltration analyses") 

parser$add_argument("--bandWidth", type="double", default=10,
                    help="width of intervals within maxDistanceFromInterface to examine for spatial analyses")

parser$add_argument("--writeCSVFiles", action="store_true", default=TRUE, 
                    help="save all intermediate data in csv files")

parser$add_argument("-NW", "--noCSVfiles", dest="writeCSVFiles", action="store_false", 
                    help="do NOT write any intermediate results to file; default=FALSE (intermediate CSV files ARE written)")

parser$add_argument("-l", "--log", type="character", default=NULL,
                    help="log file")

parser$add_argument("--setDefaultDirectoryStructure", action="store_true", default=TRUE,
                    help="set up results directories in default structure")

parser$add_argument("-ND", "--noDefaultDirectoryStructure", dest="setDefaultDirectoryStructure", action="store_false",
                    help="turn OFF the default setting to set default directory structure")

parser$add_argument("-cd", "--configDir", type="character", default=NULL, 
                    help="directory of config files, default is rootDir/studyName/config")

parser$add_argument("-db", "--debugDir", type="character", default=NULL,
                   help="debug plots will be saved to this directory, default is rootDir/studyName/debug")

parser$add_argument("--infiltrationDir", type="character", default=NULL,
                   help="all infiltration analyses will go here, default is rootDir/studyName/infiltration")

parser$add_argument("--infiltrationDensityDir", type="character", default=NULL,
                    help="all infiltration density data and plots will go here, default is infiltrationDir/density")

parser$add_argument("--infiltrationDensityFile", type="character", default=NULL,
                    help="comma separated list of CSV files containing precomputed infiltration density data") 

parser$add_argument("--infiltrationAreaDir", type="character", default=NULL,
                    help="all infiltration area data will go here, default is infiltrationDir/area")

parser$add_argument("--infiltrationAreaFile", type="character", default=NULL,
                    help="comma separated list of CSV files containing pre computed infiltration area data")

parser$add_argument("--fovStatsDir", type="character", default=NULL,
                    help="all FOV based analyses will go here, default is rootDir/studyName/fov_data")

parser$add_argument("--fovDensityDir", type="character", default=NULL,
                    help="all FOV density data and plots will go here, default is fovStatsDir/density")

parser$add_argument("--fovDensityFile", type="character", default=NULL,
                    help="comma separated list of CSV files containing pre computed FOV density data")

parser$add_argument("--fovAreaDir", type="character", default=NULL,
                    help="all FOV area data will go here, default is fovStatsDir/area")

parser$add_argument("--fovAreaFile", type="character", default=NULL,
                    help="comma separated list of CSV files containing pre computed FOV area data")

parser$add_argument("--markerConfigFile", type="character", default=NULL,
                    help="YAML file containing marker configuration (see documentation for details)")

parser$add_argument("--plotConfigFile", type="character", default=NULL,
                    help="YAML file containing plot configuration (see documentation for details)")

parser$add_argument("--cellTypeConfigFile", type="character", default=NULL,
                    help="YAML file containing configuration of cell type markers & combinations (see documentation for details)")

parser$add_argument("--annotationsDirs", type="character", default=NULL,
                    help="directory where HALO XML annotations (coordinates) files live")

parser$add_argument("--boundaryColors", type="character", default=NULL,
                    help="comma-separated, key-value pairs containing color codes that exist in HALO annotations files and the boundary type that they represent; exapmle: '65280:Tum,65535:Exc,255:Epi'")

parser$add_argument("--boundaryReassignmentFile", type="character", default=NULL,
                    help="comma-separated file indicating HALO annotations that are incorrect (see documentation for details)")

parser$add_argument("--maxG", type="double", default=5.0, help="TO DO: ADD DESCRIPTION; SOMETHING TO DO WITH CALCULATING AREA")

parser$add_argument("--debug", action="store_true", default=FALSE,
                    help="print extra output for debugging")

parser$add_argument("--updateExistingConfigFiles", action="store_true", default=FALSE,
                    help="if provided, overwrite existing study, marker and plot config files with additional and updated parameters as pipeline progresses")

parser$add_argument("--noConfigOverwrite", dest="updateExistingConfigFiles", action="store_false", default=FALSE, 
                    help="do NOT overwrite existing config files, even to add parameters")


args <- parser$parse_args()
print(args)

library(halodev)

### run configuration
configureStudy(study_name                   = args$studyName,
               study_dir                    = args$rootDir,
               set_default_directory_struct   = args$setDefaultDirectoryStruct,
               config_dir                   = args$configDir,
               study_config_file             = args$manifestFile, 
               raw_data_dir                = args$rawDataDir,
               raw_data_files              = args$rawDataFiles,
               data_dir                    = args$dataDir,
               data_files                  = args$dataFiles,
               meta_dir                    = args$metaDir,
               meta_files                  = args$metaFiles,
               drift_dir                   = args$driftDir,
               drift_files                 = args$driftFiles,
               annotations_dirs            = args$annotationsDirs,
               annotations_files           = args$annotationsFiles,
               infiltration_dir            = args$infiltrationDir,
               infiltration_density_dir    = args$infiltrationDensityDir,
               infiltration_area_dir       = args$infiltrationAreaDir,
               infiltration_density_file   = args$infiltrationDensityFile,
               infiltration_area_file      = args$infiltrationAreaFile,
               fov_stats_dir               = args$fovStatsDir,
               fov_density_dir             = args$fovDensityDir,
               fov_area_dir                = args$fovAreaDir,
               fov_density_file            = args$fovDensityFile,
               fov_area_file               = args$fovAreaFile,
               marker_analysis_config_file = args$markerConfigFile,
               cell_type_config_file       = args$cellTypeConfigFile,
               plot_config_file            = args$plotConfigFile,
               log                         = args$log,
               drift_threshold             = args$driftThreshold,
               pad                         = args$pad,
               band_width                  = args$bandWidth,
               max_distance_from_interface = args$maxDistanceFromInterface,
               debug                       = args$debug,
               updateExistingConfigFiles   = args$updateExistingConfigFiles, 
               write_csv_files             = args$writeCSVFiles)
