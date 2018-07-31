##### interface.R

#' Overwrite manifest file with updated values
#'
#' Overwrite manifest file with most current values
#'
#' @param projectParams   key-value pairs of project parameters in list format
#' @param manifest        file to which parameters should be written in YAML format
#' @return nothing
#' @export
cellDive.updateManifest <- function(projectParams, manifest){
    write(as.yaml(pp, indent=4, indent.mapping.sequence=TRUE), file=manifest)
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
    print("Loading all pre-processed halo data...")
    if(is.null(pp[[files_key]])){
        if(is.null(pp[[dir_key]])){
            stop("No data files provided. Please modify config file and rerun.")
        }
        dataFiles <- file.path(pp[[dir_key]], dir(pp[[dir_key]])[grep(".rda",dir(pp[[dir_key]]))])
    } else {
        dataFiles <- pp[[files_key]]
    }

    for(df in dataFiles){
        print(paste0("Reading file ",df))
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

    return(pp)
}


cellDive.getAllBoundaryAnnotations <- function(pp){
    allHaloAnnotations <- NULL
    updated <- FALSE
    if(!is.null(pp$annotations_file)){
        if(grepl(".txt$",pp$annotations_file)){
            print(paste0("Loading halo boundary annotations from TXT file: ",pp$annotations_file))
            haloAnnotationsTxt <- read.delim(pp$annotations_file,header=T,sep="\t")
            haloAnnotationsTxt <- split(haloAnnotationsTxt, haloAnnotationsTxt$Sample)
            for(s in names(haloAnnotationsTxt)){
                fovAnn <- split(haloAnnotationsTxt[[s]], haloAnnotationsTxt[[s]]$FOV)
                for(fov in unique(haloAnnotationsTxt)){
                    btAnn <- split(fovAnn[[fov]], fovAnn[[fov]]$RegionCode)
                    allHaloAnnotations[[s]][[fov]] <- btAnn
                }
            }
        } else if(grepl(".rda$",pp$annotations_file)){
            print(paste0("Loading halo boundary annotations from RDA file: ",pp$annotations_file))
            allHaloAnnotations <- readRDS(pp$annotations_file)
        } else {
            stop(paste0("Annotations file has unrecognized format: ",pp$annotations_file))
        }
    } else if(!is.null(pp$annotations_dirs)){
        if(basename(pp$annotations_dirs) == "HaloCoordinates"){
            pp$annotations_dirs <- dir(pp$annotations_dirs)
        }
        print("Getting all HALO boundary annotations")
        allHaloAnnTXT <- NULL
        #for(s in sampAnn$Sample_name){
        for(s in sampAnn$CELL_DIVE_ID){
            sampHaloAnn <- getAllHaloAnnotations(s, pp$annotations_dirs, boundaryColors=pp$boundary_colors)
            allHaloAnnotations[[s]] <- sampHaloAnn
        }
        for(s in names(allHaloAnnotations)){
            for(fov in names(sampHaloAnn)){
               for(bt in names(sampHaloAnn[[fov]])){
                    ann <- sampHaloAnn[[s]][[fov]][[bt]]
                    ann$Sample <- s
                    ann$FOV <- fov
                    allHaloAnnTXT <- allHaloAnnTXT %>% bind_rows(ann)
               }
            }
        }
        write.table(allHaloAnnTXT, file="haloAnnotations.txt", col.names=T, row.names=T, quote=F, sep="\t")
        updated=TRUE
    } else {
        print("No boundary annotations found. Can not run exclusions OR infiltration analysis.")
    }

    return(list(dat=allHaloAnnotations, updated=updated, pp=pp))
}


cellDive.writeMarkerComboTables <- function(pp, allDat, whichTableFile, markerConfig, updatedExclusions){
    allCountTbls <- NULL
    updated <- FALSE
    mCfg <- markerConfig
    if(!is.null(pp[[whichTableFile]]) && file.exists(pp[[whichTableFile]]) && !updatedExclusions){
        print(paste0("Marker combo table: ",pp[[whichTableFile]]," already exists. No update needed"))
    } else {

        indivMarkers <- unique(mCfg$all_cell_type_markers)
        allCountTbls <- markerComboCounts(dat=allDat, indivMarkers, oneSheet=T)
        flatMeta <- flattenMetaData(metaDir=pp$meta_dir, metaFiles=pp$meta_files)
        cellTypes <- flatMeta$MarkerCombinationAnnotation

        #for(x in 1:length(allCountTbls)){
            #s <- names(allCountTbls)[x]
            #tbl <- allCountTbls[[s]]
            tbl <- allCountTbls[[1]]
            iTbl <- interpretMarkerCombos(tbl, cellTypes)
            allCountTbls[["AllCounts"]] <- iTbl %>%
                                 select(-c(comboString,cellTypeNum,Cell_Type,Cell_Subtype)) %>%
                                 select(Markers, Cell_Type=cellType, Cell_Subtype=subtype, Total, everything()) 


            #allCountTbls[[s]] <- iTbl %>%
            #        unite("Interpretation", cellType, cellTypeNum, comboString, sep = " ") %>%
            #        dplyr::select(Count, PCT, CUM.PCT, Interpretation, everything())

            #totalsByType <- as.tibble(iTbl) %>%
            #                filter(CUM.PCT <= 0.95) %>%
            #                group_by(cellType) %>%
            #                summarize(Count=sum(Count),PCT=sum(PCT))

            #rep_na <- as.list(rep("",ncol(iTbl)))
            #names(rep_na) <- names(iTbl)
            #allCountTbls[[paste0(s,"_by_cell_type")]] <- iTbl %>%
            #                bind_rows(totalsByType) %>%
            #                replace_na(rep_na) %>%
            #                arrange(cellType) %>%
            #                unite("Interpretation", cellType, cellTypeNum, comboString, sep = " ") %>%
            #                dplyr::select(Count, PCT, CUM.PCT, Interpretation, everything())

        #}
        #uniqCT <- sort(unique(cellTypes$CellTypes$Cell_type))
        #mcw <- markerComboWorkbook(allCountTbls, uniqCT)
        if(is.null(pp[[whichTableFile]])){
            xlsxFile <- file.path(pp$study_dir, paste0(gsub("_file","",whichTableFile),"_",format(Sys.time(), "%Y%m%d"),".xlsx"))
        } else {
            xlsxFile <- pp[[whichTableFile]]
        }
        #saveWorkbook(as.data.frame(iTbl), xlsxFile)
        writeXLSXsheet(as.data.frame(allCountTbls[[1]]), xlsxFile, sheetName="AllCountTbls", append = FALSE)
        pp[[whichTableFile]] <- xlsxFile
        updated <- TRUE
    }

    return(list(dat=allCountTbls, updated=updated, pp=pp))
}


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

cellDive.calculateTotalFOVdensity <- function(allDat, markerConfig, fovAreas, pp, updatedFOVareas, updatedComboTable, fovDensityPrefix=""){
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
            fileName <- "All_samples_density.csv"
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

cellDive.calculateInfiltrationArea <- function(allDat, pp, haloAnnotations, updatedExclusions){
    updated <- FALSE

    if(is.null(pp$interface_area_file) || length(pp$interface_area_file) == 0 || !file.exists(pp$interface_area_file) ||
       is.null(pp$interface_band_assignments_file) || length(pp$interface_band_assignments_file) == 0 || 
       !file.exists(pp$interface_band_assignments_file) || updatedExclusions){

        interfaceBins <- (-(pp$max_distance_from_interface/pp$band_width):(pp$max_distance_from_interface/pp$band_width))*pp$band_width

        samp_ia <- calculateInterfaceArea(allDat, haloAnnotations=allHaloAnnotations, writeCSVfiles=pp$write_csv_files,
                                          maxG=pp$max_g, outDir=pp$infiltration_area_dir, interfaceBins=interfaceBins)

        ia <- samp_ia$area
        ba <- samp_ia$bandAssignments

        pp$interface_area_file <- file.path(pp$interface_area_dir, "all_samples_interface_area.csv")
        pp$interface_band_assignments_file <- file.path(pp$interface_area_dir, "all_samples_band_assignments.csv")
        updated <- TRUE
    } else {
        ba <- as.tibble(read.delim(pp$interface_band_assignments_file,header=T,sep=",",check.names=FALSE))
        ia <- as.tibble(read.delim(pp$interface_area_file,header=T,sep=",",check.names=FALSE))
    }

    return(list(dat=list(area=ia, bandAssignments=ba), pp=pp, updated=updated))
}

cellDive.fillInMarkerComboInterpretations <- function(xlfile=NULL, markerComboTable=NULL, outDir=getwd(), allCellTypes=NULL){
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


cellDive.calculateInfiltrationDensity <- function(markerConfig, infiltrationAreas, bandAssignments, 
                                                  pp, updatedInfiltrationAreas, 
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

    if(is.null(pp[[file_key]]) || length(pp[[file_key]])==0 || !file.exists(pp[[file_key]]) || updatedInfiltrationAreas || updatedComboTable){
        infiltrationDensity <- NULL
        for(s in unique(bandAssignments$Sample)){
            #dat <- allDat %>% filter(Sample == s)
            ia <- infiltrationAreas %>% filter(Sample == s)
            ba <- bandAssignments %>% filter(Sample == s)
            den <- calculateInfiltrationDensity(ia, ba, markerCombos, writeCSVfiles=pp$write_csv_files, 
                                                outDir=pp$infiltration_density_dir, statsByBand=TRUE)
            infiltrationDensity <- bind_rows(infiltrationDensity, den)
        }
        if(pp$write_csv_files){
            fileName <- "All_samples_density.csv"
            if(!infiltrationDensityPrefix==""){
                fileName <- paste0(infiltrationDensityPrefix,"_",fileName)
            }
            pp[[file_key]] <- file.path(pp$infiltration_density_dir,fileName)
        }
        updated <- TRUE
    } else {
        infiltrationDensity <- as.tibble(read.delim(pp[[file_key]], header=T, sep=","))
    }

    return(list(dat=infiltrationDensity, pp=pp, updated=updated))
}

