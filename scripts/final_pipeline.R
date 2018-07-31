options( java.parameters = c("-Xss2560k", "-Xmx8g") )
options('python_cmd'='/opt/common/CentOS_6-dev/python/python-3.5.1/bin/python')

write("Loading libraries...",stdout())
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("halodev"))

source("/home/byrne/halo/dev/halodev/R/constants.R")
source("/home/byrne/halo/dev/halodev/R/marker_counts.R")
source("/home/byrne/halo/dev/halodev/R/process_meta.R")
source("/home/byrne/halo/dev/halodev/R/util.R")
source("/home/byrne/halo/dev/halodev/R/exclusions.R")
source("/home/byrne/halo/dev/halodev/R/melanoma_spatial.R")
source("/home/byrne/halo/dev/halodev/R/spatial_util.R")
source("/home/byrne/halo/dev/halodev/R/validation.R")
source("/home/byrne/halo/dev/halodev/R/marker_combo_table.R")
source("/home/byrne/halo/dev/halodev/R/plots.R")
source("/home/byrne/halo/dev/halodev/R/interface.R")

###
## Parse user input
###
parser <- ArgumentParser()

## args are preferrably all in manifest
parser$add_argument("-m", "--manifest", type="character", default=NULL,
                    help="file containing all project parameters; run ?initializeProject for details")
parser$add_argument("-e", "--markExclusions", action="store_true", default=FALSE,
                    help="raw *.rda files are being given; mark exclusions based on meta data")
parser$add_argument("-s", "--start", type="character", default="pre-rethreshold", 
                    help="Options: ['pre-rethreshold'|'post-rethreshold']; starting point for pipeline")
parser$add_argument("--debug", action="store_true", default=FALSE,
                    help="print extra output for debugging")

args <- parser$parse_args()
if(!args$start %in% c("pre-rethreshold","post-rethreshold")){ stop(paste0("Unrecognized start point: ",args$start)) }
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

updatedExclusions <- FALSE
updatedComboTable <- FALSE
updatedFOVarea <- FALSE
updatedFOVdensity <- FALSE
updatedInfiltrationArea <- FALSE
updatedInfiltrationDensity <- FALSE
updatedThresholds <- FALSE

################################################
####           INITIALIZE PROJECT           ####
################################################
pp <- NULL
if(!is.null(args$manifest)){
    print("Reading study manifest...")
    pp <- read_yaml(args$manifest)
    print(pp)
} else {
    usage()
}


#pp <- validateConfig(pp) ## TO DO

################################################
#########      GET ALL ANNOTATIONS     #########
################################################
if(!is.null(pp$marker_analysis_config_file)){
    print("Reading maker config...")
    mCfg <- tryCatch({
             read_yaml(pp$marker_analysis_config_file)
          }, err=function(e){
             print(e)
          })
}

if(!is.null(pp$celltype_config_file)){
    print("Reading cell type config file...")
    ctCfg <- tryCatch({
               read_yaml(pp$celltype_config_file)
           }, err=function(e){
               print(e)
           })
}

###
## check for meta data files
###
if(!is.null(pp$meta_files)){
    metaFiles <- pp$meta_files
    print("Meta files:")
    print(paste0("  ",metaFiles))
} else if(!is.null(pp$meta_dir)){
    metaFiles <- file.path(pp$meta_dir,dir(pp$meta_dir)[grep("\\.xlsx",dir(pp$meta_dir))])
    print("Meta files:")
    print(paste0("  ",metaFiles))
} else {
    stop("No meta data given.")
}
###
## remove files that someone currently has open, indicated by "~" in basename
###
if(length(grep("~",metaFiles)) > 0){
    metaFiles <- metaFiles[-grep("~",metaFiles)]
}

###
## read study annotations
###
print("Reading sample annotations...")
sampleAnnFile <- metaFiles[grep("SampleAnnotations.xlsx",metaFiles)]
if(length(sampleAnnFile) == 0){
    stop("No sample annotation found.")
}
print("Reading FOV annotations...")
fovAnnFile <- metaFiles[grep("FOVannotations.xlsx",metaFiles)]
if(length(fovAnnFile) == 0){
    stop("No FOV annotation found.")
}
sampAnn <- as.tibble(read.xlsx(sampleAnnFile,1))
fovAnnotations <- as.tibble(read.xlsx(fovAnnFile,1))


###
### get halo boundaries
###
print("Getting halo boundaries...")
lst <- cellDive.getAllBoundaryAnnotations(pp)
allHaloAnnotations <- lst$dat
updatedBoundaries <- lst$updated
pp <- lst$pp
cellDive.updateManifest(pp, args$manifest)


################################################
#########        START  ANALYSES       #########
################################################

###
### mark exclusions
###
if(args$markExclusions || updatedBoundaries){
    print("Marking exclusions...")
    pp <- cellDive.markExclusions(pp)
    updatedExclusions <- TRUE
    cellDive.updateManifest(pp, args$manifest)
}

if(args$start == 'pre-rethreshold' || updatedExclusions){
    print("Starting pre-rethreshold analysis")
    ###
    ### read in all data with exclusions (before rethresholding)
    ###
    ## if only potting and density has already been calculated, no need to load data
    allDat <- NULL
    if(updatedExclusions || updatedBoundaries || !"raw_marker_combo_table_file" %in% names(pp) || is.null(pp$raw_marker_combo_table_file) ||
        length(pp$raw_marker_combo_table_file) == 0 || !file.exists(pp$raw_marker_combo_table_file)){
        allDat <- cellDive.loadAllData(pp,"data")
    }
    ###
    ### make RAW combination table
    ###
    lst <- cellDive.writeMarkerComboTables(pp, allDat, "raw_marker_combo_table_file", ctCfg, updatedExclusions)
    allTbls <- lst$dat
    updatedComboTable <- lst$updated
    pp <- lst$pp
    cellDive.updateManifest(pp, args$manifest)

    ###
    ### total FOV areas
    ###
    print("Getting total FOV areas...")
    lst <- cellDive.calculateTotalFOVarea(allDat, pp, metaFiles, updatedExclusions, allHaloAnnotations=allHaloAnnotations)
    fovAreas <- lst$dat
    pp <- lst$pp
    updatedFOVarea <- lst$updated
    cellDive.updateManifest(pp, args$manifest)

    ###
    ### total FOV densities
    ###
    print("Getting total FOV densities...")
    if(is.null(allDat) &&  (!"fov_density_file" %in% names(pp) || length(pp$fov_density_file) == 0 || !file.exists(pp$fov_density_file) || updatedFOVarea)) {
        allDat <- cellDive.loadAllData(pp,"data")
    } 
    parsedMarkerConfig <- getMarkerConfig(pp$marker_analysis_config_file, pp$plot_config_file)
    lst <- cellDive.calculateTotalFOVdensity(allDat, parsedMarkerConfig, fovAreas, pp, updatedFOVarea, updatedComboTable)
    fovDensity <- lst$dat
    pp <- lst$pp
    updatedFOVdensity <- lst$updated
    cellDive.updateManifest(pp, args$manifest)

    ###
    ### plot total FOV densities
    ###
    if(updatedFOVdensity){
        print("Printing total FOV density plots...")
        printTotalDensityPlots(fovDensity, fovAreas, parsedMarkerConfig, yScaleConsistency="population", absoluteDensity=TRUE,
                              densityPercentage=FALSE, byFOV=FALSE, summarize=FALSE, stacked=FALSE,
                              sampleOrder=NULL, separateLegend=TRUE, outDir=pp$fov_density_dir)
    } else {
        print("No changes to total FOV densities. No plotting necessary.")
    }

    ###
    ### infiltration 
    ###
    if(is.null(allDat) && (!"infiltration_density_file" %in% names(pp) || length(pp$infiltration_density_file) == 0 ||
        !file.exists(pp$infiltration_density_file) || updatedInfiltrationArea)){
        print("Reading all halo data...")
        allDat <- cellDive.loadAllData(pp,"data")
    }

    ###
    ### infiltration areas 
    ###
    print("Getting infiltration areas...")
    lst <- cellDive.calculateInfiltrationArea(allDat, pp, haloAnnotations, updatedExclusions)
    infiltrationAreas <- lst$dat$area
    bandAssignments <- lst$dat$bandAssignments
    pp <- lst$pp
    updatedInfiltrationArea <- lst$updated
    cellDive.updateManifest(pp, args$manifest)

    ###
    ### infiltration densities
    ###
    print("Getting infiltration densities...")
    parsedMarkerConfig <- getMarkerConfig(pp$marker_analysis_config_file, pp$plot_config_file)
    lst <- cellDive.calculateInfiltrationDensity(parsedMarkerConfig, infiltrationAreas, bandAssignments, pp, updatedInfiltrationArea)
    infiltrationDensity <- lst$dat
    pp <- lst$pp
    updatedInfiltrationDensity <- lst$updated
    cellDive.updateManifest(pp, args$manifest)

printInfiltrationDensityPlots(infiltrationDensity, bandWidth=pp$band_width, parsedMarkerConfig, yScaleConsistency="population", absoluteDensity=TRUE, densityPercentage=TRUE, byFOV=TRUE, summarize=TRUE, stacked=TRUE, sampleOrder=NULL, separateLegend=TRUE, infiltrationAreas=infiltrationAreas, bandAssignments=bandAssignments, outDir=pp$infiltration_density_dir)

}

