##### interface.R

#' Process and load any existing project data 
#' 
#' Process and load any existing project data in order to start analysis
#' where previous analyses left off
#' 
#' @param   pp  list of project parameters and values
#' @return list of all data loaded, including area, density, project configs and status
#' @export
cellDive.initCombinedProject <- function(args){

    allDat <- NULL
    updated <- list(Boundaries=FALSE,
                    Exclusions=FALSE, 
                    ComboTable=FALSE, 
                    FOVarea=FALSE, 
                    FOVdensity=FALSE, 
                    InfiltrationArea=FALSE, 
                    InfiltrationDensity=FALSE)
    dat <- list()
    pp <- NULL

    ### get project parameters
    if(!is.null(args$manifest)){
        flog.info("Reading study manifest...")
        pp <- read_yaml(args$manifest)
        flog.debug(paste(names(pp), as.vector(pp), sep=": ", collapse="\n"))
    } else {
        usage()
    }

    #pp <- validateConfig(pp) ## TO DO

    ### GET ALL ANNOTATIONS ###
    dat$ctCfg <- NULL
    if(!is.null(pp$cell_type_config_file)){
        flog.info("Reading cell type config file...")
        dat$ctCfg <- tryCatch({
                   read_yaml(pp$cell_type_config_file)
               }, err=function(e){
                   print(e)
               })
    }

    ###
    ### check for meta data files
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
    ### read study annotations
    ###
    allStudyAnn <- cellDive.loadStudyAnnotations(metaFiles)
    sampleOrder <- allStudyAnn %>% arrange(desc(`Lesion Response`), Patient) %>% pull(SampleLabel) %>% unique()
    dat$allStudyAnn <- allStudyAnn
    dat$sampleOrder <- sampleOrder


    ###
    ### get halo boundaries
    ###
    flog.info("Getting halo boundaries...")
    lst <- cellDive.getAllBoundaryAnnotations(pp, sampAnn=allStudyAnn)
    allHaloAnnotations <- lst$dat
    updated$Boundaries <- lst$updated
    pp <- lst$pp
    cellDive.updateManifest(pp, args$manifest)
    dat$allHaloAnnotations <- allHaloAnnotations


    ### 
    ### get marker/plot configs    
    ###
    flog.info("Parsing marker configurations...")
    flog.info("  (1) config for main density analysis")
    parsedMrkrCfg <- getMarkerConfig(pp$marker_analysis_config_file, pp$plot_config_file)

    if("cell_type_config_file" %in% names(pp)){
        flog.info("  (2) config for cell type analysis")
        parsedCTCfg <- getMarkerConfig(pp$cell_type_config_file, pp$plot_config_file)
        flog.info("      * combining cell type and main configs")
        completeMrkrCfg <- bind_rows(parsedMrkrCfg, parsedCTCfg)

        parsedCTCfg$Population <- gsub(",PCK26-","",parsedCTCfg$Population)
        parsedCTCfg$CellType <- gsub(",PCK26-","",parsedCTCfg$CellType)
        parsedCTCfg <- parsedCTCfg %>% unique()
    } else {
        completeMrkrCfg <- parsedMrkrCfg
    }
    completeMrkrCfg <- completeMrkrCfg %>% unique()

    parsedMrkrCfg$Population <- gsub(",PCK26-","",parsedMrkrCfg$Population)
    parsedMrkrCfg$CellType <- gsub(",PCK26-","",parsedMrkrCfg$CellType)
    parsedMrkrCfg <- parsedMrkrCfg %>% unique()

    dat$parsedMarkerConfig <- parsedMrkrCfg
    dat$completeMarkerConfig <- completeMrkrCfg


    if(args$fovStats){
        if(is.null(allDat) && (!"fov_area_file" %in% names(pp) || length(pp$fov_area_file) == 0 ||
            !file.exists(pp$fov_area_file) || !"fov_density_file" %in% names(pp) || length(pp$fov_density_file) == 0 ||
            !file.exists(pp$fov_density_file))){
            flog.info("Reading all halo data...")
            allDat <- cellDive.loadAllData(pp,"data")
        }


        ###
        ### total FOV areas
        ###
        flog.info("Getting total FOV areas...")
        lst <- cellDive.calculateTotalFOVarea(allDat, pp, metaFiles, updated$Exclusions, 
                                              allHaloAnnotations=allHaloAnnotations)
        fovAreas <- lst$dat %>% unique() %>% filter(Area != "Area") # *specific to combined data to ensure proper
                                                                    # combining of area files (for this one I did
                                                                    # that manually; if combining cohorts is going
                                                                    # to be a thing, eventually write code to do it
                                                                    # properly
        fovAreas$SPOT <- as.numeric(fovAreas$SPOT)
        fovAreas$Area <- as.numeric(fovAreas$Area)
        fovAreas <- fovAreas %>% left_join(allStudyAnn, by=c("Sample","SPOT"))
        pp <- lst$pp
        updated$FOVarea <- lst$updated
        cellDive.updateManifest(pp, args$manifest)

        ###
        ### total FOV densities
        ###
        flog.info("Getting total FOV densities...")
        lst <- cellDive.calculateTotalFOVdensity(allDat, completeMrkrCfg, fovAreas, pp, 
                                                 updated$FOVarea, updated$ComboTable)
        fovDensity <- lst$dat %>% filter(Sample != "Sample") # *see note above about manual combining of area files; 
                                                             # same applies here
        #fovDensity$SPOT <- as.numeric(fovDensity$SPOT)
        pp <- lst$pp
        updated$FOVdensity <- lst$updated
        cellDive.updateManifest(pp, args$manifest)

        ### ** SPECIFIC TO COMBINED MELANOMA DATA ** ###
        fovDensity <- fovDensity %>% left_join(allStudyAnn, by=c("Sample","SPOT"))
        fovDensity$Sample <- fovDensity$SampleLabel
        fovDensity$CellType <- gsub(",PCK26-","",fovDensity$CellType) # *specific to combined melanoma data; if combining
                                                                      # cohorts with different markers is going to be a 
                                                                      # thing, write code to do this properly
        ## this must be done AFTER density 
        fovAreas$Sample <- fovAreas$SampleLabel

        dat$fovAreas <- fovAreas
        dat$fovDensity <- fovDensity
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
        lst <- cellDive.calculateInfiltrationArea(allDat, pp, allHaloAnnotations, updated$Exclusions)
        infiltrationAreas <- lst$dat$area %>% filter(Sample != "Sample")
        infiltrationAreas$SPOT <- as.numeric(infiltrationAreas$SPOT)
        infiltrationAreas$BandArea <- as.numeric(infiltrationAreas$BandArea)
        infiltrationAreas <- infiltrationAreas %>% left_join(allStudyAnn, by=c("Sample","SPOT"))
        infiltrationAreas$Sample <- infiltrationAreas$SampleLabel
        infiltrationAreas <- infiltrationAreas %>% select(Sample, SPOT, Band, BandArea)

        bandAssignments <- lst$dat$bandAssignments %>% left_join(allStudyAnn, by=c("Sample","SPOT"))
        bandAssignments$Sample <- bandAssignments$SampleLabel

        pp <- lst$pp
        updated$InfiltrationArea <- lst$updated
        cellDive.updateManifest(pp, args$manifest)

        ###
        ### infiltration densities
        ###
        flog.info("Getting infiltration densities...")
        lst <- cellDive.calculateInfiltrationDensity(parsedMrkrCfg, infiltrationAreas, bandAssignments, 
                                                     allStudyAnn, pp, updated$InfiltrationArea)
        infiltrationDensity <- lst$dat
        infiltrationDensity$Sample <- infiltrationDensity$SampleLabel
        infiltrationDensity$CellType <- gsub(",PCK26-","",infiltrationDensity$CellType)

        pp <- lst$pp
        updated$InfiltrationDensity <- lst$updated
        cellDive.updateManifest(pp, args$manifest)

        dat$infiltrationDensity <- infiltrationDensity
        dat$infiltrationAreas <- infiltrationAreas
        dat$bandAssignments <- bandAssignments
    }

    return(list(dat=dat, pp=pp, updated=updated))
}





#' Overwrite manifest file with updated values
#'
#' Overwrite manifest file with most current values
#'
#' @param projectParams   key-value pairs of project parameters in list format
#' @param manifest        file to which parameters should be written in YAML format
#' @return nothing
#' @export
cellDive.updateManifest <- function(projectParams, manifest){
    write(as.yaml(projectParams, indent=4, indent.mapping.sequence=TRUE), file=manifest)
}

#' Load data from each rda file into one table
#' 
#' Read each data file, select the minimum columns required for downstream
#' analyses and store as one table in memory. Columns include: 
#'     Sample
#'     SubSample
#'     SPOT
#'     UUID
#'     Marker
#'     XMin/Max
#'     YMin/Max
#'     Value
#' Data is filtered for cells whose EXCLUDE value == "" and ValueType == "Positive"
#' @param projectParams   key-value pairs of project parameters in list format
#' @param whichFileSet    name of key (minus "_dir" and "_files" suffixes) of files to 
#'                        be loaded (e.g., 'raw_data' or 'data' or 'rethresholded_data') 
#' @return tibble containing all data required for downstream analysis
#' @export
cellDive.loadAllData <- function(projectParams, whichFileSet){
    dir_key <- paste0(whichFileSet,"_dir")
    files_key <- paste0(whichFileSet,"_files") 
    pp <- projectParams
    allDat <- NULL
    flog.info("Loading all pre-processed halo data...")
    if(is.null(pp[[files_key]])){
        if(is.null(pp[[dir_key]])){
            stop("No data files provided. Please modify config file and rerun.")
        }
        dataFiles <- file.path(pp[[dir_key]], dir(pp[[dir_key]])[grep(".rda",dir(pp[[dir_key]]))])
    } else {
        dataFiles <- pp[[files_key]]
    }

    for(df in dataFiles){
        flog.info(paste0("Reading file ",df))
        dat <- readRDS(df)
        dat <- dat %>%
               filter(EXCLUDE=="", ValueType=="Positive") %>%
               select(Sample, SubSample, SPOT, UUID, Marker, matches("XM|YM"), Value)
        ## remove cells that are NOT DAPI pos
        dapiNeg <- dat$UUID[which(dat$Marker == "DAPI" && dat$Value == 0)]
        if(length(dapiNeg) > 0){
            dat <- dat %>% filter(!UUID %in% dapiNeg)
        }
        allDat <- allDat %>%
                   bind_rows(dat)
    }
    return(allDat)
}

#' Load and parse all Sample+FOV annotations
#' 
#' Load and parse all Sample+FOV annotations
#' 
#' @param metaFiles  vector of meta data files, including at least
#'                   *SampleAnnotations.xlsx and *FOVannotations.xlsx
#' @return tibble of joined Sample+FOV annotations
#' @export
cellDive.loadStudyAnnotations <- function(metaFiles){

    flog.info("Reading FOV annotations...")
    fovAnnFile <- metaFiles[grep("FOVannotations.xlsx",metaFiles)]
    if(length(fovAnnFile) == 0){
        stop("No FOV annotation found.")
    }
    fovAnnotations <- as.tibble(read.xlsx(fovAnnFile,1,check.names=FALSE))

    flog.info("Reading Sample annotations...")
    sampleAnnFile <- metaFiles[grep("SampleAnnotations.xlsx",metaFiles)]
    if(length(sampleAnnFile) == 0){
        stop("No sample annotation found.")
    }
    sampAnn <- as.tibble(read.xlsx(sampleAnnFile,1,check.names=FALSE))
    sampAnn$rda_sample <- sampAnn$CELL_DIVE_ID
    sampAnn$rda_sample[which(sampAnn$rda_sample == "melanoma_untreated")] <- "Untreated"
    sampAnn$rda_sample[which(sampAnn$rda_sample == "melanoma_PR")] <- "PR"
    sampAnn$rda_sample[which(sampAnn$rda_sample == "melanoma_CR")] <- "CR"

    allStudyAnn <- full_join(sampAnn, fovAnnotations, by=c("CELL_DIVE_ID")) %>%
                   select(Sample=rda_sample, Patient, Sample_number, FOV_number, FOV_type,
                          `Lesion Response`, Treatment, `Prior Rx (systemic)`, `Num. IL2 injections`,
                          `IL2 interval (days)`, `Lesion Response` ) %>%
                   mutate(SampleLabel=paste(paste0("Pt_",Patient),paste0("s_",Sample_number),sep=".")) %>%
                   rename(SPOT=FOV_number)

    return(allStudyAnn)
}

#' Mark cells in *.rda file to be excluded from analysis
#'
#' Add 'EXCLUDE' column to each sample tibble indicating reason(s) why a cell should
#' be excluded from analysis. Possible reasons include: 
#'    (1) cell falls outside padding border
#'    (2) lab determined FOV and/or Marker should be excluded for technical reasons
#'    (3) drift/loss percentage is above maximum allowed value
#'    (4) cell was determined by Halo to fall inside an area to be excluded
#' Write updated tibble to *.rda file with "_Excl.rda" appended to original file name
#' @param projectParams  a list of parameters for entire project including meta directory/files, drift directory, halo annotations directories, and raw data directory
#' @return updated projectParams including data directory where new *.rda files are saved
#' @export
cellDive.markExclusions <- function(projectParams){
    pp <- projectParams
    ## check for data
    if(!is.null(pp$raw_data_files)){
        dataFiles <- pp$raw_data_files
    } else if(!is.null(pp$raw_data_dir)){
        dataFiles <- file.path(pp$raw_data_dir,dir(pp$raw_data_dir)[grep("\\.rda",dir(pp$raw_data_dir))])
    } else {
        stop("No data given.")
    }
    flog.debug("Data files:")
    flog.debug(paste0("  ",dataFiles,"\n"))

    ## check for drift/loss summaries
    if(!is.null(pp$drift_files)){
        driftFiles <- pp$drift_files
        flog.debug("Drift files:")
        flog.debug(paste0("  ",driftFiles))
    } else if(!is.null(pp$drift_dir)){
        driftFiles <- file.path(pp$drift_dir,dir(pp$drift_dir)[grep(".*summary.*txt",dir(pp$drift_dir))])
        flog.debug("Drift files:")
        flog.debug(paste0("  ",driftFiles))
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

        flog.info(paste0("Reading file ",df))
        dat <- readRDS(df)
        #dat$Sample <- gsub("_.*","",dat$Sample) #### FIGURE OUT A BETTER WAY

        ## to map FOV annotations to samples
        samp <- unique(dat$Sample)
        #cellDiveId <- sampAnn %>% filter(Sample_name == samp) %>% pull(CELL_DIVE_ID) %>% as.character()
        cellDiveId <- samp   #### RIGHT NOW THIS IS INCONSISTENT BETWEEN COHORTS

        drift <- NULL
        ## figure out which drift file to use
        if(!is.null(driftFiles)){
            dlFile <- driftFiles[grep(paste0(cellDiveId,"_|",samp,"_"), driftFiles, ignore.case=TRUE)]
            if(length(dlFile) == 0){
                warning("No drift file found. Skipping drift exclusions")
            } else {
                flog.info(paste0("Reading drift file",dlFile))
                drift <- as.tibble(read.delim(dlFile, header=T, sep="\t"))
            }
        }

        sampHaloAnnotations <- allHaloAnnotations[[samp]]

        flog.info(paste0("Marking exclusions for ",samp))
        dat <- markExclusions(samp, dat, drift, fovAnnnotations, cellDiveId, haloAnn=sampHaloAnnotations,
                              printPlots=TRUE, debugDir=pp$debug_dir, boundaryColors=pp$boundary_colors,
                              boundaryReassignmentFile=pp$boundaryReassignmentFile)

        saveRDS(dat, file=file.path(outDir,gsub("\\.rda","_Excl.rda",basename(df))))

        ## update project param list to include new data files, which will be used
        ## from here on instead of the raw files 
        pp$data_files <- c(pp$data_files, file.path(outDir,gsub("\\.rda","_Excl.rda",basename(df))))
    }

    return(pp)
}


#' Load and/or parse all Halo boundary annotations
#' 
#' If boundary annotations have already been parsed and saved to *.rda file, load file. Otherwise, parse all annotations, save to *.rda file and return parsed annotations.
#'
#' @param pp        list of all project parameters, from which it will be determined whether annotations have been previously parsed. List must include annotations_dirs and/or annotations_file parameters.
#' @param sampAnn   parsed sample annotations; only required if pp$annotations_file is NULL or the file does not exist/is empty
#' @return list containing three values: 
#'           (1) dat = list of all Halo boundary annotations organized by sample and FOV
#'           (2) updated = logical indicating whether annotations were just parsed (TRUE) or previously parsed (FALSE)
#'           (3) pp = updated list of project parameters
#' @export
cellDive.getAllBoundaryAnnotations <- function(pp,sampAnn=NULL){
    allHaloAnnotations <- NULL
    updated <- FALSE
    if("annotations_file" %in% names(pp) && (!is.null(pp$annotations_file) || !file.exists(pp$annotations_file) || length(pp$annotations_file == 0))){
        if(grepl(".rda$",pp$annotations_file)){
            flog.info(paste0("Loading halo boundary annotations from RDA file: ",pp$annotations_file))
            allHaloAnnotations <- readRDS(pp$annotations_file)
        } else {
            stop(paste0("Annotations file has unrecognized format: ",pp$annotations_file))
        }
    } else if("annotations_dirs" %in% names(pp) && !is.null(pp$annotations_dirs)){
        if(basename(pp$annotations_dirs) == "HaloCoordinates"){
            pp$annotations_dirs <- file.path(pp$annotations_dirs, dir(pp$annotations_dirs))
        }
        flog.info("Getting all HALO boundary annotations")
        for(s in sampAnn$CELL_DIVE_ID){
            sampHaloAnn <- getAllHaloAnnotations(s, pp$annotations_dirs, boundaryColors=pp$boundary_colors)
            allHaloAnnotations[[s]] <- sampHaloAnn
        }
        if(!is.null(pp$annotations_file) && length(pp$annotations_file) == 1){
            saveRDS(allHaloAnnotations, file=pp$annotations_file)
        } else {
            saveRDS(allHaloAnnotations, file="allHaloAnnotations.rda")
            pp$annotations_file <- "allHaloAnnotations.rda"
        }
        updated=TRUE
    } else {
        flog.info("No boundary annotations found. Can not run exclusions OR infiltration analysis.")
    }

    return(list(dat=allHaloAnnotations, updated=updated, pp=pp))
}

#' Write XLSX table of counts for all marker combinations in all samples
#' 
#' Write XLSX table of counts for all marker combinations in all samples
#' 
#' @param pp                list of all project parameters
#' @param allDat            tibble of all Halo data
#' @param whichTableFile    name of parameter in pp to look at when checking to see if table already exists
#' @param markerConfig      parsed marker configuration (value returned from getMarkerConfig())
#' @param updatedExclusions logical indicating whether exclusions have been updated
#' @return list of three items:  
#'           (1) dat = list containing one item, the count table generated
#'           (2) updated = logical indicating whether a new count table was generated (TRUE) or an old one was loaded (FALSE) 
#'           (3) pp = updated list of project parameters
#' @export
cellDive.writeMarkerComboTables <- function(pp, allDat, whichTableFile, markerConfig, updatedExclusions){
    allCountTbls <- NULL
    updated <- FALSE
    mCfg <- markerConfig
    if(!is.null(pp[[whichTableFile]]) && file.exists(pp[[whichTableFile]]) && !updatedExclusions){
        flog.info(paste0("Marker combo table: ",pp[[whichTableFile]]," already exists. No update needed"))
    } else {

        indivMarkers <- unique(mCfg$all_cell_type_markers)
        allCountTbls <- markerComboCounts(dat=allDat, indivMarkers, oneSheet=T)
        flatMeta <- flattenMetaData(metaDir=pp$meta_dir, metaFiles=pp$meta_files)
        cellTypes <- flatMeta$MarkerCombinationAnnotation

        tbl <- allCountTbls[[1]]
        iTbl <- interpretMarkerCombos(tbl, cellTypes)
        allCountTbls[["AllCounts"]] <- iTbl %>%
                                 select(-c(comboString,cellTypeNum,Cell_Type,Cell_Subtype)) %>%
                                 select(Markers, Cell_Type=cellType, Cell_Subtype=subtype, Total, everything()) 

        if(is.null(pp[[whichTableFile]])){
            xlsxFile <- file.path(pp$study_dir, paste0(gsub("_file","",whichTableFile),"_",format(Sys.time(), "%Y%m%d"),".xlsx"))
        } else {
            xlsxFile <- pp[[whichTableFile]]
        }
        writeXLSXsheet(as.data.frame(allCountTbls[[1]]), xlsxFile, sheetName="AllCountTbls", append = FALSE)
        pp[[whichTableFile]] <- xlsxFile
        updated <- TRUE
    }

    return(list(dat=allCountTbls, updated=updated, pp=pp))
}

#' Calculate area of each entire FOV
#' 
#' Generate a table of three columns: Sample, SPOT, Area
#' 
#' @param allDat             tibble of all Halo data
#' @param pp                 list of all project parameters
#' @param metaFiles          vector of meta files
#' @param updatedExclusions  logical indicating whether exclusions were updated during current run of pipeline
#' @param allHaloAnnotations list of parsed Halo annotations organized by sample and then FOV
#' @param fovAreaPrefix      prefix added to 'fov_area_file' and/or 'fov_area_dir' in names of pp list (leave "" if no prefix)
#' @return list of three items:  
#'           (1) dat = tibble of three columns: Sample, SPOT, Area 
#'           (2) updated = logical indicating whether a new area table was generated (TRUE) or an old one was loaded (FALSE) 
#'           (3) pp = updated list of project parameters
#' @export
cellDive.calculateTotalFOVarea <- function(allDat, pp, metaFiles, updatedExclusions, allHaloAnnotations=NULL, fovAreaPrefix=""){
    fovAreas <- NULL
    updated <- FALSE
    file_key <- "fov_area_file"
    dir_key <- "fov_area_dir"
    if(!fovAreaPrefix == ""){
        file_key <- paste0(fovAreaPrefix,"_",file_key)
        dir_key <- paste0(fovAreaPrefix,"_",dir_key)
    }
    if(is.null(pp[[file_key]]) || length(pp[[file_key]])==0 || !file.exists(pp[[file_key]]) || updatedExclusions){
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
                                          annotationsDirs=annotationsDirs, writeCSVfiles=pp$write_csv_files, maxG=pp$max_g)
        if(pp$write_csv_files){
            fileName <- "All_samples_area.csv"
            if(!fovAreaPrefix==""){
                fileName <- paste0(fovAreaPrefix,"_",fileName)
            }
            pp[[file_key]] <- file.path(pp$fov_area_dir,fileName)
        }
        updated <- TRUE
    } else {
        fovAreas <- as.tibble(read.delim(pp[[file_key]], header=T, sep=","))
    }

    return(list(dat=fovAreas, updated=updated, pp=pp))
}

#' Calculate density of markers in each full FOV
#' 
#' Generate table of four columns: Sample, SPOT, Counts, Density
#' 
#' @param allDat          all Halo data
#' @param markerConfig    tibble of parsed marker configuration (value returned from getMarkerConfig)
#' @param fovAreas        tibble of areas for each FOV in each Sample (value returned from cellDive.calculateTotalFOVarea)
#' @param pp              list of all project parameters
#' @param updatedFOVareas logical indicating whether FOV areas were updated during current pipeline run
#' @param updatedComboTable  logical indicating whether combo table was updated during current pipeline run
#' @param fovDensityPrefix  prefix added to prefix added to 'fov_density_file' and/or 'fov_density_dir' in names of pp list (leave "" if no prefix) 
#' @return list of three items:  
#'           (1) dat = tibble of three columns: Sample, SPOT, Count, Density 
#'           (2) updated = logical indicating whether a new density table was generated (TRUE) or an old one was loaded (FALSE) 
#'           (3) pp = updated list of project parameters
#' @export
cellDive.calculateTotalFOVdensity <- function(allDat, markerConfig, fovAreas, pp, updatedFOVarea, updatedComboTable, fovDensityPrefix=""){
    mCfg <- markerConfig
    fovDensity <- NULL
    updated <- FALSE
    file_key <- "fov_density_file"
    dir_key <- "fov_density_dir"
    if(!fovDensityPrefix == ""){
        file_key <- paste0(fovDensityPrefix,"_",file_key)
        dir_key <- paste0(fovDensityPrefix,"_",dir_key)
    }

    markerCombos <- NULL
    if("cell_type_marker_combinations" %in% names(mCfg) && !is.null(mCfg$cell_type_marker_combinations)){
        markerCombos <- mCfg$cell_type_marker_combinations
    } else if("CellType" %in% names(mCfg) && !is.null(mCfg$CellType)){
        markerCombos <- unique(mCfg$CellType)
    } else {
        warning("Could not parse marker config. No densities calculated.")
        return(list(dat=NULL, pp=pp, updated=FALSE))
    }

    if(is.null(pp[[file_key]]) || length(pp[[file_key]])==0 || !file.exists(pp[[file_key]]) || updatedFOVarea || updatedComboTable){
        fovDensity <- calculateMarkerDensityTotalFOV(allDat, fovAreas, markerCombos,
                                                     writeCSVfiles=pp$write_csv_files, densityDir=pp$fov_density_dir)    
        if(pp$write_csv_files){
            fileName <- "all_samples_density.csv"
            if(!fovDensityPrefix==""){
                fileName <- paste0(fovDensityPrefix,"_",fileName)
            }
            pp[[file_key]] <- file.path(pp$fov_density_dir,fileName)
        }
        updated <- TRUE
    } else {
        fovDensity <- as.tibble(read.delim(pp[[file_key]], header=T, sep=","))
    }
    
    return(list(dat=fovDensity, pp=pp, updated=updated))
}

#' Calculate infiltration area
#' 
#' Generate table with columns: Sample, SPOT, Band, Area
#' 
#' @param allDat           tibble of all halo data
#' @param pp               list of all project parameters
#' @param haloAnnotations  list of parsed halo annotations organized by sample and then by FOV
#' @param updatedExclusions  logical indicating whether exclusions (haloAnnotations) were updated during current run of the pipeline
#' @return list of three items:
#'           (1) dat = list where area = table of Sample,SPOT,Area and bandAssignments = table of 0/1 for each marker in each cell, plus a column "Band" indicating roughly how far the cell lies from a tumor interface
#'           (2) updated = logical indicating whether table being returned was generated (TRUE) or loaded from existing file (FALSE)
#'           (3) pp = updated list of project parameters
#' @export 
cellDive.calculateInfiltrationArea <- function(allDat, pp, haloAnnotations, updatedExclusions){
    updated <- FALSE

    if(is.null(pp$infiltration_area_file) || length(pp$infiltration_area_file) == 0 || !file.exists(pp$infiltration_area_file) ||
       is.null(pp$infiltration_band_assignments_file) || length(pp$infiltration_band_assignments_file) == 0 || 
       !file.exists(pp$infiltration_band_assignments_file) || updatedExclusions){

        interfaceBins <- (-(pp$max_distance_from_interface/pp$band_width):(pp$max_distance_from_interface/pp$band_width))*pp$band_width

        samp_ia <- calculateInterfaceArea(allDat, haloAnnotations=haloAnnotations, writeCSVfiles=pp$write_csv_files,
                                          maxG=pp$max_g, outDir=pp$infiltration_area_dir, interfaceBins=interfaceBins)

        ia <- samp_ia$area %>% gather(3:ncol(samp_ia$area), key="Band", value="BandArea")
        ba <- samp_ia$bandAssignments

        pp$infiltration_area_file <- file.path(pp$infiltration_area_dir, "all_samples_interface_area.csv")
        pp$infiltration_band_assignments_file <- file.path(pp$infiltration_area_dir, "all_samples_band_assignments.csv")
        updated <- TRUE
    } else {
        ba <- as.tibble(read.delim(pp$infiltration_band_assignments_file,header=T,sep=",",check.names=FALSE))
        ia <- as.tibble(read.delim(pp$infiltration_area_file,header=T,sep=",",check.names=FALSE)) 
        if(!"BandArea" %in% names(ia)){
            ia <- ia %>% gather(3:ncol(ia), key="Band", value="BandArea")
        }
    }

    return(list(dat=list(area=ia, bandAssignments=ba), pp=pp, updated=updated))
}

#' Add columns Cell_Type and Cell_Subtype to marker counts file
#' 
#' Given a table where the first three columns = Markers, Cell_Type, and Cell_Subtype and Cell_Type/Subtype are blank, and the remaining columns contain total counts for each marker combination in a single sample, fill in blank columns using meta data files
#'
#' @param xlfile            XLSX file containing table described above; default=NULL; if not specified, markerComboTable must be given
#' @param markerComboTable  marker combination table returned from cellDive.getMarkerCombinationTable(); if not specified, must provide xlfile
#' @param outDir            directory where xlsx file should be saved
#' @param allCellTypes      parsed cell types as returned by flattenMetaData()
#' @return marker combo table with columns Cell_Type and Cell_Subtype filled in
#' @export
cellDive.fillInMarkerComboInterpretations <- function(xlfile=NULL, markerComboTable=NULL, outDir=NULL, allCellTypes=NULL){
    if(is.null(outDir)){ outDir <- getwd() }
    if(is.null(markerComboTable)){
        if(is.null(xlfile)){ stop("Please provide either xlfile or markerComboTable") }
        markerComboTable <- as.tibble(read.xlsx(xlfile,1))
    } 
    markerComboTbl2 <- interpretMarkerCombos(markerComboTable, allCellTypes)
    markerComboTbl2$Cell_Type <- markerComboTbl2$cellType
    markerComboTbl2$Cell_Subtype <- markerComboTbl2$subtype

    markerComboTbl3 <- markerComboTbl2 %>% select(-c(cellType,subtype,cellTypeNum,comboString))
    markerComboTbl3[is.na(markerComboTbl3)] <- ""
    writeXLSXsheet(as.data.frame(markerComboTbl3), file.path(outDir,gsub(".xlsx","_withInterpretations.xlsx",xlfile)), sheetName="Sheet1", append=FALSE)
    return(markerComboTbl3)
}

#' Calculate infiltration density
#'
#' Generate table with five columns: Sample,SPOT,Band,Counts,Density
#' 
#' @param markerConfig             parsed marker configuration (as returned by getMarkerConfig())
#' @param infiltrationAreas        area table containing Sample,SPOT,Band,Area
#' @param bandAssignments          table containing Band column where each value is an assignment for a single cell to a specific distance band around a tumor interface
#' @param pp                       list of all project parameters
#' @param sampAnn                  table of all sample annotations
#' @param updatedInfiltrationAreas logical indicating whether infiltration areas were updated during the current pipeline run
#' @param infiltrationDensityPrefix  prefix added to infiltration_density_file/dir in list of project params; leave equal to "" to indicate no prefix added
#' @return list of three items:
#'           (1) dat = table including Sample,SPOT,Band,Counts,Density 
#'           (2) updated = logical indicating whether table being returned was generated (TRUE) or loaded from existing file (FALSE)
#'           (3) pp = updated list of project parameters
#' @export 
cellDive.calculateInfiltrationDensity <- function(markerConfig, infiltrationAreas, bandAssignments, 
                                                  sampAnn, pp, updatedInfiltrationAreas, 
                                                  infiltrationDensityPrefix=""){
    mCfg <- markerConfig
    infiltrationDensity <- NULL
    updated <- FALSE
    file_key <- "infiltration_density_file"
    dir_key <- "infiltration_density_dir"
    if(!infiltrationDensityPrefix == ""){
        file_key <- paste0(infiltrationDensityPrefix,"_",file_key)
        dir_key <- paste0(infiltrationDensityPrefix,"_",dir_key)
    }

    markerCombos <- NULL
    if("cell_type_marker_combinations" %in% names(mCfg) && !is.null(mCfg$cell_type_marker_combinations)){
        markerCombos <- mCfg$cell_type_marker_combinations
    } else if("CellType" %in% names(mCfg) && !is.null(mCfg$CellType)){
        markerCombos <- unique(mCfg$CellType)
    } else {
        warning("Could not parse marker config. No densities calculated.")
        return(list(dat=NULL, pp=pp, updated=FALSE))
    }

    if(is.null(pp[[file_key]]) || length(pp[[file_key]])==0 || !file.exists(pp[[file_key]]) || updatedInfiltrationAreas){
        fileName <- file.path(pp$infiltration_density_dir,"all_samples_density.csv")
        if(!infiltrationDensityPrefix==""){
            fileName <- file.path(pp$infiltration_density_dir,paste0(infiltrationDensityPrefix,"_",fileName))
        }

        flog.info("Calculating infiltration density for all samples")

        infiltrationDensity <- calculateInfiltrationDensity(infiltrationAreas, bandAssignments, 
                                                            markerCombos, fileName, writeCSVfiles=FALSE) 

        if(pp$write_csv_files){
            pp[[file_key]] <- fileName 
            write_csv(as.data.frame(infiltrationDensity), fileName) 
        }
        updated <- TRUE
    } else {
        infiltrationDensity <- as.tibble(read.delim(pp[[file_key]], header=T, sep=",")) %>%
                               filter(Sample != "Sample")
    }

    if(!is.null(sampAnn)){
        infiltrationDensity <- infiltrationDensity %>% 
                               left_join(sampAnn, by=c("Sample","SPOT"))
    }
    return(list(dat=infiltrationDensity, pp=pp, updated=updated))
}

