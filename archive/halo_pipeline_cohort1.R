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

###
# Parse user input
###
parser <- ArgumentParser()

## args are preferrably all in manifest
parser$add_argument("-m", "--manifest", type="character", default=NULL,
                    help="file containing all project parameters; run ?initializeProject for details")
parser$add_argument("-e", "--markExclusions", action="store_true", default=FALSE,
                    help="raw *.rda files are being given; mark exclusions based on meta data")
parser$add_argument("--debug", action="store_true", default=FALSE,
                    help="print extra output for debugging")

args <- parser$parse_args()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

updatedExclusions <- FALSE
updatedComboTable <- FALSE
updatedFOVarea <- FALSE
updatedFOVdensity <- FALSE
updatedInfiltrationArea <- FALSE
updatedInfiltrationDensity <- TRUE

################################################
####           INITIALIZE PROJECT           ####
################################################
pp <- NULL
print(args$manifest)
if(!is.null(args$manifest)){
    pp <- initializeProject(args$manifest,type="counts")
    print(pp)
} else {
    usage()
}

#pp <- validateConfig(pp)

################################################
#########      GET ALL ANNOTATIONS     #########
################################################
if(!is.null(pp$marker_config_file)){
    mCfg <- tryCatch({
             read_yaml(pp$marker_config_file)
          }, err=function(e){
             print(e)
          })
}

## check for meta data files
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
## remove files that someone currently has open, indicated by "~" in basename
if(length(grep("~",metaFiles)) > 0){
    metaFiles <- metaFiles[-grep("~",metaFiles)]
}

## read study annotations
sampleAnnFile <- metaFiles[grep("SampleAnnotations.xlsx",metaFiles)]
if(length(sampleAnnFile) == 0){
    stop("No sample annotation found.")
}
sampAnn <- as.tibble(read.xlsx(sampleAnnFile,1))

fovAnnFile <- metaFiles[grep("FOVannotations.xlsx",metaFiles)]
if(length(fovAnnFile) == 0){
    stop("No FOV annotation found.")
}
fovAnnotations <- as.tibble(read.xlsx(fovAnnFile,1))


## read halo boundary annotations
allHaloAnnotations <- NULL
if(!is.null(pp$annotations_file)){
    allHaloAnnotations <- readRDS(pp$annotations_file)
} else if(!is.null(pp$annotations_dirs)){
    if(length(pp$annotations_dirs) == 1 && basename(pp$annotations_dirs) == "HaloCoordinates"){
        pp$annotations_dirs = file.path(pp$annotations_dirs, dir(pp$annotations_dirs))
    }
    print("Getting all HALO boundary annotations")
    for(s in sampAnn$Sample_name){
        sampHaloAnn <- getAllHaloAnnotations(s, pp$annotations_dirs, boundaryColors=pp$boundary_colors)
        allHaloAnnotations[[s]] <- sampHaloAnn
    }
    saveRDS(allHaloAnnotations, file="allHaloAnnotations.rda")
}
if(is.null(allHaloAnnotations)){
    stop("No boundary annotations found. Can not run exclusions OR infiltration analysis.")
}


################################################
####             MARK EXCLUSIONS            ####
################################################               
if(args$markExclusions){
    ## check for data
    if(!is.null(pp$raw_data_files)){
        dataFiles <- pp$raw_data_files
    } else if(!is.null(pp$raw_data_dir)){
        dataFiles <- file.path(pp$raw_data_dir,dir(pp$raw_data_dir)[grep("\\.rda",dir(pp$raw_data_dir))])
    } else {
        stop("No data given.")
    }
    print("Data files:")
    print(paste0("  ",dataFiles))

    ## check for drift/loss summaries
    if(!is.null(pp$drift_files)){
        driftFiles <- pp$drift_files
        print("Drift files:")
        print(paste0("  ",driftFiles))
    } else if(!is.null(pp$drift_dir)){
        driftFiles <- file.path(pp$drift_dir,dir(pp$drift_dir)[grep(".*summary.*txt",dir(pp$drift_dir))])
        print("Drift files:")
        print(paste0("  ",driftFiles))
    } else {
        driftFiles <- NULL
        warning("No drift summaries given.")
    }

    ## figure out where to write rda files with exclusion column
    outDir <- getwd()
    if("data_dir" %in% names(pp) && !is.null(pp$data_dir)){
        outDir <- pp$data_dir
    } else if("out_dir" %in% names(pp) && !is.null(pp$out_dir)){
        outDir <- pp$out_dir
    }
    if(!file.exists(outDir)){
        dir.create(outDir, recursive=T, showWarnings=FALSE)
    }

    for(x in seq(dataFiles)){
        df <- dataFiles[x]

        print(paste0("Reading file ",df))
        dat <- readRDS(df)
        dat$Sample <- gsub("_.*","",dat$Sample) #### FIGURE OUT A BETTER WAY

        ## to map FOV annotations to samples
        samp <- unique(dat$Sample)
        cellDiveId <- sampAnn %>% filter(Sample_name == samp) %>% pull(CELL_DIVE_ID) %>% as.character()
        #cellDiveId <- samp   #### RIGHT NOW THIS IS INCONSISTENT BETWEEN COHORTS

        drift <- NULL
        ## figure out which drift file to use
        if(!is.null(driftFiles)){
            dlFile <- driftFiles[grep(paste0(cellDiveId,"_|",samp,"_"), driftFiles, ignore.case=TRUE)]
            if(length(dlFile) == 0){
                warning("No drift file found. Skipping drift exclusions")
            } else {
                print(paste0("Reading drift file",dlFile))
                drift <- as.tibble(read.delim(dlFile, header=T, sep="\t"))
            }
        }

        sampHaloAnnotations <- allHaloAnnotations[[samp]]        

        print(paste0("Marking exclusions for ",samp))
        dat <- markExclusions(samp, dat, drift, fovAnnnotations, cellDiveId, haloAnn=sampHaloAnnotations,
                              printPlots=TRUE, debugDir=pp$debug_dir, boundaryColors=pp$boundary_colors,
                              boundaryReassignmentFile=pp$boundaryReassignmentFile)

        saveRDS(dat, file=file.path(outDir,gsub("\\.rda","_Excl.rda",basename(df))))

        ## update project param list to include new data files, which will be used
        ## from here on instead of the raw files 
        pp$data_files <- c(pp$data_files, file.path(outDir,gsub("\\.rda","_Excl.rda",basename(df))))   
    }

    updatedExclusions <- TRUE
}

################################################
####               READ DATA                ####
################################################
dataFiles <- NULL
allDat <- NULL
if(is.null(pp$data_files)){
    if(is.null(pp$data_dir)){
        stop("No data files provided. Please modify config file and rerun.")
    }
    dataFiles <- file.path(pp$data_dir, dir(pp$data_dir)[grep(".rda",dir(pp$data_dir))])
} else {
    dataFiles <- pp$data_files
}

for(df in dataFiles){
    dat <- readRDS(df)
    dat <- dat %>% 
           filter(EXCLUDE=="", ValueType=="Positive") %>%
           select(Sample, SubSample, SPOT, UUID, Marker, matches("XM|YM"), Value)
    allDat <- allDat %>%
              bind_rows(dat)
}

################################################
####         MARKER COMBINATION TABLE       ####
################################################
allTbls <- NULL
if(!is.null(pp$raw_marker_combo_table_file) && !updatedExclusions){
    wb <- loadWorkbook(pp$raw_marker_combo_table_file)
    sheets <- getSheets(wb)
    for(sheet in names(sheets)){
        allTbls[[sheet]] <- as.tibble(as.data.frame(sheets[[sheet]]))
    }
} else {
    if(is.null(mCfg)){
        stop("Marker configuration is NULL. Can not generate marker table.")
    }
    indivMarkers <- unique(mCfg$all_cell_type_markers)
    allCountTbls <- markerComboCounts(dat=allDat, indivMarkers)

    for(x in 1:length(allCountTbls)){
        s <- names(allCountTbls)[x]
        tbl <- allCountTbls[[s]]

        flatMeta <- flattenMetaData(metaDir=pp$meta_dir, metaFiles=pp$meta_files)
        cellTypes <- flatMeta$MarkerCombinationAnnotation
        iTbl <- interpretMarkerCombos(tbl, cellTypes)

        allCountTbls[[s]] <- iTbl %>%
                unite("Interpretation", cellType, cellTypeNum, comboString, sep = " ") %>%
                dplyr::select(Count, PCT, CUM.PCT, Interpretation, everything())

        totalsByType <- as.tibble(iTbl) %>%
                        filter(CUM.PCT <= 0.95) %>%
                        group_by(cellType) %>%
                        summarize(Count=sum(Count),PCT=sum(PCT))

        rep_na <- as.list(rep("",ncol(iTbl)))
        names(rep_na) <- names(iTbl)
        allCountTbls[[paste0(s,"_by_cell_type")]] <- iTbl %>%
                        bind_rows(totalsByType) %>%
                        replace_na(rep_na) %>%
                        arrange(cellType) %>%
                        unite("Interpretation", cellType, cellTypeNum, comboString, sep = " ") %>%
                        dplyr::select(Count, PCT, CUM.PCT, Interpretation, everything())

    }
    uniqCT <- sort(unique(cellTypes$CellTypes$Cell_type))
    mcw <- markerComboWorkbook(allCountTbls, uniqCT)
    xlsxFile <- file.path(pp$study_dir, paste0("raw_marker_combo_table_",format(Sys.time(), "%Y%m%d"),".xlsx"))
    saveWorkbook(mcw, xlsxFile)

    pp$raw_marker_combo_table_file <- xlsxFile

    updatedComboTable <- TRUE
}     

write(as.yaml(pp, indent=4, indent.mapping.sequence=TRUE), file=args$manifest)

################################################
####         CALCULATE TOTAL FOV AREA       ####  
################################################

if(is.null(pp$fov_area_file) || !file.exists(pp$fov_area_file) || updatedExclusions){
    if(is.null(allHaloAnnotations) && is.null(pp$annotations_dirs)){
        flatMeta <- flattenMetaData(metaDir=pp$meta_dir, metaFiles=pp$meta_files)
        annotationsDirs <- convertBoundaryFilePaths(unique(flatMeta$SlideAnnotation$Path_to_boundary_files),"")
        annotationsDirs <- file.path(annotationsDirs,dir(annotationsDirs))
        pp$annotations_dirs <- annotationsDirs
    } else {
        annotationsDirs <- NULL
    }
    ## TO DO: ADD THESE CHECKS TO VALIDATION FUNCTION
    if(!"write_csv_files" %in% names(pp) || is.null(pp$write_csv_files)){ pp$write_csv_files = TRUE }
    if(!"max_g" %in% names(pp) || is.null(pp$max_g)){ pp$max_g = 5 } 
    #if(is.null(mCfg)){ # error }

    fovAreas <- calculateAreaTotalFOV(allDat, metaFiles, pp$fov_area_dir, haloAnnotations=allHaloAnnotations,
                                      annotationsDirs=NULL, writeCSVfiles=pp$write_csv_files, maxG=pp$max_g)
    if(pp$write_csv_files){
        pp$fov_area_file <- file.path(pp$fov_area_dir,"All_samples_area.csv")
    }
    updatedFOVarea <- TRUE
} else {
    fovAreas <- as.tibble(read.delim(pp$fov_area_file, header=T, sep=","))
}

if(is.null(pp$fov_density_file) || !file.exists(pp$fov_density_file) || updatedFOVarea){
    fovDensity <- calculateMarkerDensityTotalFOV(allDat, fovAreas, mCfg$cell_type_marker_combinations, 
                                                 writeCSVfiles=pp$write_csv_files, densityDir=pp$fov_density_dir)
    if(pp$write_csv_files){
        pp$fov_density_file <- file.path(pp$fov_density_dir,"All_samples_density.csv")
    }
    updatedFOVdensity <- TRUE
} else {
    fovDensity <- as.tibble(read.delim(pp$fov_density_file, header=T, sep=","))
}

write(as.yaml(pp, indent=4, indent.mapping.sequence=TRUE), file=args$manifest)

# 4) infiltration area/density
#if(is.null(pp$infiltration_area_file) || !file.exists(pp$infiltration_area_file || updatedExclusions){
        
#    infArea <- calculateInterfaceArea(iFiles, allDat, writeCSVfiles=pp$write_csv_files, 
#                                      outDir=pp$infiltration_area_dir)
#} else {

#}


# 5) total FOV density plots

# 6) infiltration density plots


# 7) if updateConfigFile, overwrite config file used with final config values
