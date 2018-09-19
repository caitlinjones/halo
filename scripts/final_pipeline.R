
options(stringsAsFactors = FALSE)
options( java.parameters = c("-Xss2560k", "-Xmx8g") )
options('python_cmd'='/opt/common/CentOS_6-dev/python/python-3.5.1/bin/python')

write("Loading libraries...",stdout())
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("halodev"))

#source("/home/byrne/halo/dev/halodev/R/constants.R")
#source("/home/byrne/halo/dev/halodev/R/process_meta.R")
#source("/home/byrne/halo/dev/halodev/R/util.R")
#source("/home/byrne/halo/dev/halodev/R/exclusions.R")
#source("/home/byrne/halo/dev/halodev/R/melanoma_spatial.R")
#source("/home/byrne/halo/dev/halodev/R/spatial_util.R")
#source("/home/byrne/halo/dev/halodev/R/validation.R")
#source("/home/byrne/halo/dev/halodev/R/marker_combo_table.R")
#source("/home/byrne/halo/dev/halodev/R/plots.R")
#source("/home/byrne/halo/dev/halodev/R/interface.R")

if(!interactive()){
    ###
    ## Parse user input
    ###
    parser <- ArgumentParser()

    ## args are preferrably all in manifest
    parser$add_argument("-m", "--manifest", type="character", default=NULL,
                        help="file containing all project parameters; run ?initializeProject for details")
    parser$add_argument("-e", "--markExclusions", action="store_true", default=FALSE,
                        help="raw *.rda files are being given; mark exclusions based on meta data")
    parser$add_argument("-f", "--fovStats", action="store_true", default=FALSE,
                        help="run stats on total-FOV level")
    parser$add_argument("-i", "--infiltrationStats", action="store_true", default=FALSE,
                        help="run stats on interface level (infiltration)")
    parser$add_argument("-fp", "--forceAllPlotting", action="store_true", default=FALSE,
                        help="rerun all plotting regardless of whether data has changed")
    parser$add_argument("--noFOVplots", action="store_true", default=FALSE,
                        help="turn off plotting of total FOV-based data")
    parser$add_argument("--noInfiltrationPlots", action="store_true", default=FALSE,
                        help="turn off plotting of interface-based data")
    parser$add_argument("--debug", action="store_true", default=FALSE,
                        help="print extra output for debugging")

    args <- parser$parse_args()
} else {
    ### for most plotting tests
    args <- list(manifest="config/study_config.yaml", 
                 fovStats=FALSE, 
                 noFOVplots=FALSE, 
                 infiltrationStats=TRUE, 
                 noInfiltrationPlots=TRUE, 
                 markExclusions=FALSE)
    ### for macrophage pies
    #args <- list(manifest="config/macrophage_analysis_config.yaml",
    #             fovStats=TRUE,
    #             noFOVplots=TRUE,
    #             infiltrationStats=FALSE,
    #             markExclusions=FALSE)

    #### for COMBINED cohorts
    args <- list(manifest="config/study_config.yaml",
                 fovStats=TRUE,
                 infiltrationStats=TRUE,
                 markExclusions=FALSE)
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

updatedBoundaries <- FALSE
updatedExclusions <- FALSE
updatedComboTable <- FALSE
updatedFOVarea <- FALSE
updatedFOVdensity <- FALSE
updatedInfiltrationArea <- FALSE
updatedInfiltrationDensity <- FALSE

################################################
####           INITIALIZE PROJECT           ####
################################################
pp <- NULL
if(!is.null(args$manifest)){
    flog.info("Reading study manifest...")
    pp <- read_yaml(args$manifest)
    flog.debug(paste(names(pp), as.vector(pp), sep=": ", collapse="\n"))
} else {
    usage()
}


#pp <- validateConfig(pp) ## TO DO

################################################
#########      GET ALL ANNOTATIONS     #########
################################################
if(!is.null(pp$cell_type_config_file)){
    flog.info("Reading cell type config file...")
    ctCfg <- tryCatch({
               read_yaml(pp$cell_type_config_file)
           }, err=function(e){
               print(e)
           })
}

###
## check for meta data files
###
if(!is.null(pp$meta_files)){
    metaFiles <- pp$meta_files
    flog.info("Meta files:")
    flog.info(paste0("  ",metaFiles))
} else if(!is.null(pp$meta_dir)){
    metaFiles <- file.path(pp$meta_dir,dir(pp$meta_dir)[grep("\\.xlsx",dir(pp$meta_dir))])
    flog.info("Meta files:")
    flog.info(paste0("  ",metaFiles))
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
flog.info("Reading sample annotations...")
sampleAnnFile <- metaFiles[grep("SampleAnnotations.xlsx",metaFiles)]
if(length(sampleAnnFile) == 0){
    stop("No sample annotation found.")
}
flog.info("Reading FOV annotations...")
fovAnnFile <- metaFiles[grep("FOVannotations.xlsx",metaFiles)]
if(length(fovAnnFile) == 0){
    stop("No FOV annotation found.")
}
sampAnn <- as.tibble(read.xlsx(sampleAnnFile,1))
fovAnnotations <- as.tibble(read.xlsx(fovAnnFile,1))

### right now Cohort1 Sample is equal to sampAnn$Sample_name, not cell_dive_id
if(pp$studyName == "Cohort1"){
    sampAnn$CELL_DIVE_ID <- c("PR","CR","Untreated") 
}

###
### get halo boundaries
###
flog.info("Getting halo boundaries...")
lst <- cellDive.getAllBoundaryAnnotations(pp, sampAnn=sampAnn)
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
    flog.info("Marking exclusions...")
    pp <- cellDive.markExclusions(pp)
    updatedExclusions <- TRUE
    cellDive.updateManifest(pp, args$manifest)
}

flog.info("Starting analysis")
###
### read in all data with exclusions (before rethresholding)
###
## if only plotting and density has already been calculated, no need to load data
allDat <- NULL
if(updatedExclusions || updatedBoundaries || 
   is.null(pp$raw_marker_combo_table_file) ||
   length(pp$raw_marker_combo_table_file) == 0 || !file.exists(pp$raw_marker_combo_table_file) ||
   (args$fovStats && 
       (is.null(pp$fov_area_file) || length(pp$fov_area_file) == 0 || !file.exists(pp$fov_area_file) || 
        is.null(pp$fov_density_file) || length(pp$fov_density_file) == 0) || !file.exists(pp$fov_density_file)) ||
   (args$infiltrationStats && 
       (is.null(pp$infiltration_area_file) || length(pp$infiltration_area_file) == 0 || !file.exists(pp$infiltration_area_file) ||
        is.null(pp$infiltration_density_file) || length(pp$infiltration_density_file) == 0 || !file.exists(pp$infiltration_density_file))
   )
  ){
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


#if(pp$studyName == "Cohort1"){
#    sampleOrder <- c("Untreated","PR","CR")
#}

if(args$fovStats){

    sampLabelMapping <- sampAnn %>% select(Sample=CELL_DIVE_ID, Sample_name, Response, Patient) %>%
                        arrange(desc(Response), Patient) %>%
                        mutate(SampleLabel=paste(Sample,Sample_name,sep=","))
    sampleOrder <- unique(sampLabelMapping$SampleLabel)

    ###
    ### total FOV areas
    ###
    flog.info("Getting total FOV areas...")
    lst <- cellDive.calculateTotalFOVarea(allDat, pp, metaFiles, updatedExclusions, allHaloAnnotations=allHaloAnnotations)
    fovAreas <- lst$dat
    pp <- lst$pp
    updatedFOVarea <- lst$updated
    cellDive.updateManifest(pp, args$manifest)

    ###
    ### total FOV densities
    ###
    flog.info("Getting total FOV densities...")
    if(is.null(allDat) &&  (!"fov_density_file" %in% names(pp) || length(pp$fov_density_file) == 0 || !file.exists(pp$fov_density_file) || updatedFOVarea)) {
        allDat <- cellDive.loadAllData(pp,"data")
    }     
    parsedMarkerConfig <- getMarkerConfig(pp$marker_analysis_config_file, pp$plot_config_file)
    lst <- cellDive.calculateTotalFOVdensity(allDat, parsedMarkerConfig, fovAreas, pp, updatedFOVarea, updatedComboTable)
    fovDensity <- lst$dat
    pp <- lst$pp
    updatedFOVdensity <- lst$updated
    cellDive.updateManifest(pp, args$manifest)

    fovDensity <- fovDensity %>% left_join(sampLabelMapping, by=c("Sample"))
    fovDensity$Sample <- fovDensity$SampleLabel

    fovAreas <- fovAreas %>% left_join(sampLabelMapping, by=c("Sample"))
    fovAreas$Sample <- fovAreas$SampleLabel


    ###
    ### plot total FOV densities
    ###
    if(!args$noFOVplots && (updatedFOVdensity || args$forceAllPlotting)){
        flog.info("Printing total FOV density plots...")
        for(ms in unique(parsedMarkerConfig$MarkerSet)){
            msConfig <- parsedMarkerConfig %>% filter(MarkerSet == ms)
            den <- fovDensity %>% filter(CellType %in% unique(msConfig$CellType))        
            printTotalFOVDensityPlots(den, fovAreas, msConfig, yScaleConsistency="population", 
                               absoluteDensity=TRUE, densityPercentage=TRUE, byFOV=TRUE, summarize=TRUE, 
                               stacked=TRUE, sampleOrder=sampleOrder, outDir=pp$fov_density_dir)
        }

        ### plot cell identity densities
        if("cell_type_config_file" %in% names(pp)){
            parsedCellTypeConfig <- getMarkerConfig(pp$cell_type_config_file, pp$plot_config_file)
        }

        cellIdDen <- summarizeFOVDataByCellTypeDefinition(fovDensity, parsedCellTypeConfig, fovAreas, pp$plot_config_file)
        den <- cellIdDen$den
        den$Sample <- factor(den$Sample, levels=sampleOrder)

        printTotalFOVDensityPlots(den, fovAreas, cellIdDen$config, yScaleConsistency="population",
                               absoluteDensity=TRUE, densityPercentage=TRUE, byFOV=TRUE, summarize=TRUE,
                               stacked=TRUE, sampleOrder=sampleOrder, outDir=pp$fov_density_dir, forceStack=TRUE)


    } else {
        flog.info("No changes to total FOV densities. No plotting necessary.")
    }
}

if(args$infiltrationStats){
    ###
    ### infiltration 
    ###
    if(is.null(allDat) && (!"infiltration_area_file" %in% names(pp) || length(pp$infiltration_area_file) == 0 ||
        !file.exists(pp$infiltration_area_file))){
        flog.info("Reading all halo data...")
        allDat <- cellDive.loadAllData(pp,"data")
    }

    ###
    ### infiltration areas 
    ###
    flog.info("Getting infiltration areas...")
    lst <- cellDive.calculateInfiltrationArea(allDat, pp, allHaloAnnotations, updatedExclusions)
    infiltrationAreas <- lst$dat$area
    bandAssignments <- lst$dat$bandAssignments
    pp <- lst$pp
    updatedInfiltrationArea <- lst$updated
    cellDive.updateManifest(pp, args$manifest)

    ###
    ### infiltration densities
    ###
    flog.info("Getting infiltration densities...")
    parsedMarkerConfig <- getMarkerConfig(pp$marker_analysis_config_file, pp$plot_config_file)
    lst <- cellDive.calculateInfiltrationDensity(parsedMarkerConfig, infiltrationAreas, bandAssignments, sampAnn, pp, updatedInfiltrationArea)
    infiltrationDensity <- lst$dat
    pp <- lst$pp
    updatedInfiltrationDensity <- lst$updated
    cellDive.updateManifest(pp, args$manifest)


    ###
    ### plot infiltration densities
    ###
    if(!args$noInfiltrationPlots && (updatedInfiltrationDensity || args$forceAllPlotting)){

        #### PREP LABELS FOR PLOTS; CURRENTLY SPECIFIC TO MELANOMA COHORT2 DATA
        if(pp$studyName == "Cohort2"){
            smps <- unlist(lapply(as.vector(unique(infiltrationDensity$Sample)), function(x){ unlist(strsplit(x,","))[1] }))
            smpOrdr <- order(unlist(lapply(unique(smps), function(x){ as.numeric(unlist(strsplit(x,"_"))[2]) })))

            sampleOrder <- unique(infiltrationDensity$Sample)[smpOrdr]
            sampleLabelOrder <- unique(infiltrationDensity$SampleLabel)[smpOrdr]
            infiltrationDensity$Sample <- factor(infiltrationDensity$SampleLabel, levels=sampleLabelOrder)

            infiltrationDensity$Response[which(infiltrationDensity$Response == "Untreated")] <- "UT"

            sampLabelMapping <- sampAnn %>%
                                select(CELL_DIVE_ID, Sample_name) %>%
                                mutate(Sample=CELL_DIVE_ID, SampleLabel=paste(CELL_DIVE_ID, Sample_name, sep=",")) %>%
                                select(-c(CELL_DIVE_ID, Sample_name))

            infiltrationAreas <- infiltrationAreas %>%
                                 left_join(sampLabelMapping, by=c("Sample")) %>%
                                 mutate(Sample=SampleLabel) %>%
                                 select(-(SampleLabel))
        }
        ########


        #### PREP LABELS for COHORT 1 plots
        ##### THIS IS TERRIBLE AND TEMPORARY AND NEEDS TO BE RESOLVED IN META DATA FILES
        if(pp$studyName == "Cohort1"){
            infiltrationDensity$Response[which(infiltrationDensity$Response == "NA")] <- "UT"
            infiltrationDensity$Response[is.na(infiltrationDensity$Response)] <- "UT"
            infiltrationDensity$Sample_name <- paste(paste0("Pt",infiltrationDensity$Patient),infiltrationDensity$Response,sep="_")
            infiltrationDensity$SampleLabel <- paste(infiltrationDensity$Sample,infiltrationDensity$Sample_name, sep=",")

            smpOrdr <- c(2,1)

            sampLabelMapping <- infiltrationDensity %>% select(Sample,SampleLabel) %>% unique()
            infiltrationDensity$Sample <- factor(infiltrationDensity$SampleLabel, levels=unique(infiltrationDensity$SampleLabel)[smpOrdr])
            sampleOrder <- unique(as.vector(infiltrationDensity$Sample))[smpOrdr]


            infiltrationAreas <- infiltrationAreas %>%
                                 left_join(sampLabelMapping, by=c("Sample")) %>%
                                 mutate(Sample=SampleLabel) %>%
                                 select(-(SampleLabel))

        }

        if(!"BandArea" %in% names(infiltrationAreas)){
            infiltrationAreas <- infiltrationAreas %>% gather(3:ncol(infiltrationAreas), key="Band", value="BandArea")
        }   

        ###
        ### plot infiltration densities
        ###

        byFOV <- TRUE
        indivPopulations <- TRUE
        ptsPerPage <- NULL
        sumCols <- 3
        sumRows <- 3

        ### NOTE: it is not necessary to do this for each marker set (i.e., you can pass the entire
        ###       parsedMarkerConfig to the printInfiltrationDensityPlots() function, but this way, 
        ###       plots will be printed for each marker set before generating plots for the next set
        for(ms in unique(parsedMarkerConfig$MarkerSet)){
            msConfig <- parsedMarkerConfig %>% filter(MarkerSet == ms)
            den <- infiltrationDensity %>% filter(CellType %in% unique(msConfig$CellType))
            printInfiltrationDensityPlots(den, bandWidth=pp$band_width, msConfig, 
                   yScaleConsistency="population", absoluteDensity=TRUE, densityPercentage=TRUE, byFOV=byFOV, 
                   summarize=TRUE, stacked=TRUE, sampleOrder=sampleOrder, infiltrationAreas=infiltrationAreas, 
                   outDir=pp$infiltration_density_dir, indivPopulations=indivPopulations, sampleSummaryPtsPerPage=ptsPerPage,
                   sampleSummaryCols=sumCols,sampleSummaryRows=sumRows)
        }


        ### plot cell identity densities
        if("cell_type_config_file" %in% names(pp)){
            parsedCellTypeConfig <- getMarkerConfig(pp$cell_type_config_file, pp$plot_config_file)
        }

        sampIDs <- sampAnn %>% select(Sample=CELL_DIVE_ID, Sample_name, Patient, Response)
        sampIDs$Response[is.na(sampIDs$Response)] <- "UT"
        sampIDs$Response[which(sampIDs$Response == "NA")] <- "UT"

        if(pp$studyName == "Cohort1"){
            # for Cohort 1
            sampIDs$Sample_name <- paste(paste0("Pt",sampIDs$Patient), sampIDs$Response,sep="_")
        }

        ba <- bandAssignments %>% left_join(sampIDs, by=c("Sample")) %>%
              mutate(SampleLabel=paste0(Sample,",",Sample_name))
        ba$Sample <- ba$SampleLabel

        cellIdDen <- summarizeByCellTypeDefinition(ba, infiltrationAreas, parsedCellTypeConfig, pp$plot_config_file)
        den <- cellIdDen$den
        den$Sample <- factor(den$Sample, levels=sampleOrder)
        tmp <- infiltrationDensity %>% 
               filter(CellType %in% parsedCellTypeConfig$CellType) %>% 
               select(Sample,Sample_name,Patient,Response) %>%
               unique()
        den <- den %>% left_join(tmp, by=c("Sample"))


        printInfiltrationDensityPlots(den, bandWidth=pp$band_width, cellIdDen$config, 
                   yScaleConsistency="population", absoluteDensity=TRUE, densityPercentage=TRUE, byFOV=byFOV, 
                   summarize=TRUE, stacked=TRUE, sampleOrder=sampleOrder, infiltrationAreas=infiltrationAreas, 
                   outDir=pp$infiltration_density_dir, forceStack=TRUE, indivPopulations=indivPopulations, 
                   sampleSummaryPtsPerPage=ptsPerPage,sampleSummaryCols=sumCols,sampleSummaryRows=sumRows)    
    }
}

