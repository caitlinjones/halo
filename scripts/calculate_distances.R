#!/home/byrne/R/R-3.4.3/bin/Rscript

options(stringsAsFactors = FALSE)
options( java.parameters = c("-Xss2560k", "-Xmx8g") )
options('python_cmd'='/opt/common/CentOS_6-dev/python/python-3.5.1/bin/python')

write("Loading libraries...",stdout())
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("halodev"))

#source("/home/byrne/halo/dev/halodev/R/contstants.R")
#source("/home/byrne/halo/dev/halodev/R/melanoma_spatial.R")
#source("/home/byrne/halo/dev/halodev/R/exclusions.R")
#source("/home/byrne/halo/dev/halodev/R/process_meta.R")
#source("/home/byrne/halo/dev/halodev/R/interface.R")

pp <- NULL

if(interactive()){
    pp <- list(#manifest="/home/byrne/halo/data/Melanoma_IL2__Final/CombinedCohorts/study_config.yaml",
               data_dir="/home/byrne/halo/data/Melanoma_IL2__Final/CombinedCohorts/objectAnalysisData",
               meta_dir="/ifs/tcga/socci/Multiomyx/HaloData/Melanoma_IL2__Final/CombinedData",
               meta_files=NULL,
               band_width=10,
               max_distance_from_interface=360,
               writeCSVfiles=TRUE,
               infiltration_dir="/home/byrne/halo/data/Melanoma_IL2__Final/CombinedCohorts/infiltration",
               annotations_file="/home/byrne/halo/data/Melanoma_IL2__Final/CombinedCohorts/allHaloAnnotations.rda"
               )
} else {
    parser <- ArgumentParser()
    
    ## args are preferably all in manifest 
    parser$add_argument("-m", "--manifest", type="character", default=NULL,
                        help="file containing all project parameters; run ?initializeProject for details")
    args <- parser$parse_args()

    pp <- read_yaml(args$manifest)
}

### check for meta data files
if(!is.null(pp$meta_files)){
    metaFiles <- pp$meta_files
} else if(!is.null(pp$meta_dir)){
    metaFiles <- file.path(pp$meta_dir,dir(pp$meta_dir)[grep("\\.xlsx",dir(pp$meta_dir))])
} else {
    stop("No meta data given.")
}
## remove files that someone currently has open, indicated by "~" in basename
if(length(grep("~",metaFiles)) > 0){ metaFiles <- metaFiles[-grep("~",metaFiles)] }

### read study annotations
allStudyAnn <- cellDive.loadStudyAnnotations(metaFiles)

### get halo boundaries
allHaloAnnotations <- cellDive.getAllBoundaryAnnotations(pp, sampAnn=allStudyAnn)$dat

### load halo data
allDat <- cellDive.loadAllData(pp, "data")

### add distances to data
interfaceBins <- ((-pp$max_distance_from_interface/pp$band_width):(pp$max_distance_from_interface/pp$band_width))*10 
distances <- calcDistancesFromTumorBoundary(allDat, haloAnnotations=allHaloAnnotations, 
                                writeCSVfiles=pp$writeCSVfiles, outDir=pp$infiltration_dir, interfaceBins=interfaceBins)

### remove marker info?
#distances <- distances %>% select(Sample,SPOT,UUID,X,Y,Distance,Band)
