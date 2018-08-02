#' Build a list of all marker combinations described in "complex" combination
#' format in cell types spreadsheet (described in docs)
#'
#' Based on combination type ("ALL", "ANY", "<2", ">3", etc.), generate all
#' possible combinations of markers given
#' 
#' @param comboType  tells what kind of combinations to build ("ALL" indicates
#'                   all markers given must stay together, "ANY" means any combination
#'                   of given markers, including all different lengths, and 
#'                   types like ">2" or "<=3" tells how many markers each combination must 
#'                   contain
#' @param markers    vector of markers to combine
#' @return a list of all possible combinations matching criteria in comboType
getComplexMarkerCombos <- function(comboType, markers){
    #if(markers == "TOTAL"){
        ## handle differently
    #}

    markers <- unlist(strsplit(markers,","))
    if(comboType == "ALL"){
        nums <- length(markers)
    } else if(comboType == "ANY") {
        nums <- c(1:length(markers))
    } else if( length(grep(">=",comboType)) > 0 ){
        num <- as.integer(gsub(">=","",comboType))
        nums <- c(num:length(markers))
    } else if( length(grep("<=",comboType)) > 0 ){
        num <- as.integer(gsub("<=","",comboType))
        nums <- c(1:num)
    } else if( length(grep(">",comboType)) > 0 ){
        num <- as.integer(gsub(">","",comboType))
        nums <- c((num+1):length(markers))
    } else if( length(grep("<",comboType)) > 0 ){
        num <- as.integer(gsub("<","",comboType))
        nums <- c(1:(num-1))
    }
    mSets <- list()
    for(i in nums){
        cmb <- t(combn(sort(markers),m=i))
        for(c in 1:nrow(cmb)){
            mSets[[length(mSets)+1]] <- cmb[c,]
        }
    }
    return(mSets)
}

#' Get all possible marker combinations for a certain cell type/subtype
#'
#' Create list of all marker combinations that define a given cell type and/or
#' subtype
#' 
#' @param cellTypeTable  parsed table of all marker combinations, cell types and subtypes
#' @param cellType       the cell type to filter for
#' @param cellSubtype    the cell subtype to filter for
#' @return a vector of cell type marker combinations
#' @export
getAllCellTypeCombinations <- function(cellTypeTable, cellType, cellSubtype=NULL){
    if(is.null(cellSubType)){
        combos <- cellTypeTable$MarkerCombination %>% filter(Cell_type == cellType)
    } else {
        combos <- cellTypeTable$MarkerCombination %>% filter(Cell_type == cellType & Subtype == cellSubtype)
    }
    return(combos)
}

#' Make list of all known cell types (marker combinations) given in meta data files
#' 
#' Make list of all known cell types (marker combinations) given in meta data files
#' 
#' @param simpleTypes   straight-forward markercombination-cellType pairs
#' @param complexTypes  table describing combinations that are not exactly straight-forward
#'                      (described in meta data format documentation)
#' @param unreasonableCombinations **NOTE: WILL MOST LIKELY BE ELIMINULLTED**
#' @return list of two elements: CellTypes and Conditional, where CellTypes are ALL known 
#'           markercombination-cellType pairs and Conditional is a vector describing cell types
#'           that are conditional based on the data itself
#' @export
getCellTypes <- function(simpleTypes, complexTypes){ #, unreasonableCombinations){
    #unCombos <- unreasonableCombinations
    conditionalTypes <- NULL

    cellTypes <- simpleTypes

    ## sort combos from simple sheet file
    for(x in 1:nrow(cellTypes)){
        combo <- unlist(strsplit(cellTypes$Marker_combination[x],","))
        negs <- unlist(strsplit(cellTypes$Negatives[x],","))
        cellTypes$Marker_combination[x] <- paste(c(sort(combo),negs),collapse=",")
    }

    ## parse complex types into "simple" types by generating all possible
    ## combinations for the type given
    for(x in 1:nrow(complexTypes)){
        ct <- complexTypes[x,]
        ## pull out conditional types to be parsed later
        if("TOTAL" %in% c(ct$OF,ct$OF.1)){
            conditionalTypes <- rbind(conditionalTypes,ct)
            next
        }
        set1 <- getComplexMarkerCombos(ct$IF,ct$OF)
        set2 <- getComplexMarkerCombos(ct$MATCHES,ct$OF.1)
        allCombos <- expand.grid(set1,set2)
        finalCombos <- c()
        for(i in 1:length(allCombos[[1]])){
            combo <- unique(sort(unlist(c(allCombos[[1]][i],allCombos[[2]][i]))))
            finalCombos <- c(finalCombos,paste(combo,collapse=","))
        }
        finalCombos <- unique(paste0(finalCombos,",",ct$NEGATIVES))
        overlaps <- finalCombos[which(finalCombos %in% cellTypes$Marker_combination)]
        if(length(overlaps) > 0){
            stop(paste0("\n\nThe following 'Complex' marker combinations also occur in the list of 'Simple' combinations: ",paste(overlaps,collapse="; "),". Please correct and rerun.\n\n"))
        }
        cts <- data.frame(Marker_combination = finalCombos,
                          Negatives = rep(ct$NEGATIVES,length(finalCombos)),
                          With_OR_Without = "", 
                          Cell_type = rep(ct$CELL_TYPE,length(finalCombos)),
                          Subtype = rep(ct$CELL_TYPE,length(finalCombos)))
        cellTypes <- rbind(cellTypes,cts)
    }

    ## parse unreasonable combos matrix to get the "Unknown" combinations
    #rownames(unCombos) <- unCombos[,1]
    #unCombos <- unCombos[,-1]
    #allUnknowns <- c()
    #for(x in rownames(unCombos)){
    #    unknowns <- colnames(unCombos)[which(unCombos[x,] == "Unknown")]
    #    if(!is.null(unknowns) & length(unknowns) > 0){
    #        combos <- unlist(lapply(unknowns,function(y){ paste(sort(c(y,x)),collapse=",") }))
    #        overlaps <- combos[which(combos %in% cellTypes$Marker_combination)]
    #        if(length(overlaps) > 0){
    #            stop(paste0("\n\nThe following 'Unknown' marker combinations also occur in the CellTypes xlsx file: ",paste(overlaps,collapse="; "),". Please correct and rerun.\n\n"))
    #        }
    #        allUnknowns <- c(allUnknowns,combos)
    #    }
    #}
    #uc <- data.frame(Marker_combination = sort(allUnknowns),
    #                 Cell_type = rep("Unknown", length(allUnknowns)),
    #                 Subtype = rep("Unknown", length(allUnknowns)))
    #cellTypes <- rbind(cellTypes,uc)

    return(list(CellTypes=cellTypes,Conditional=conditionalTypes))
}

#' Based on MarkerCombinationAnnotation table from flattenMetaData(), assign a cell type
#' to a combination in marker combination table
#' 
#' Based on MarkerCombinationAnnotation table from flattenMetaData(), assign a cell type
#' to a combination in marker combination table
#' 
#' @param combinationRow   table row from marker combination table
#' @param cellTypes        table describing known marker combinations, containing at least
#'                         Marker_combination and Cell_type columns
#' @param conditionalTypes table of conditional cell types (*NOTE* this is not implemented yet)
#' @return character string indicating cell type of given combination
#' @export
assignCellType <- function(combinationRow,cellTypes,conditionalTypes,totalCellCount){
    
    if(names(combinationRow)[1] == "Markers"){
        allMarkers <- unique(gsub("-","",unlist(strsplit(cellTypes$Marker_combination,","))))
        combo <- combinationRow[1,"Markers"] %>% pull
        markers <- sort(unlist(strsplit(combo,",")))
        if(length(which(allMarkers %in% markers)) > 0){
            negs <- allMarkers[-which(allMarkers %in% markers)] 
        } else {
            negs <- NULL
        }
    } else {
        markerVals <- combinationRow[,4:ncol(combinationRow)]
        markers <- names(markerVals)[which(markerVals == "X")]
        negs <- names(markerVals)[-which(markerVals == "X")]
        combo <- paste0(sort(c(markers, paste0(negs,"-"))),collapse=",")
    }

    allCTcombos <- lapply(cellTypes$Marker_combination, function(x){ sort(unlist(strsplit(x,","))) })

    cellType <- NULL
    if(is.null(markers) || length(markers) == 0 || markers == "AllNeg"){
        return(list(cellType="NEGATIVE",subtype=""))
    }
    for(x in 1:length(allCTcombos)){
        actc <- sort(allCTcombos[[x]])
        if(length(grep("-",actc)) > 0){
            pos <- actc[-grep("-",actc)]
            requiredNeg <- gsub("-","",actc[grep("-",actc)])
        } else {
            pos <- actc
            requiredNeg <- NULL
        }
        if(!is.null(negs) && !is.null(requiredNeg) && length(negs) > 0 && length(requiredNeg) > 0){
            if(all(requiredNeg %in% negs) && all(pos %in% markers)){
                if(!is.null(cellType)){
                    print(paste0("Warning: ",combo," is ambiguous according to cell types file"))
                }               
                cellType=cellTypes$Cell_type[x]
                subtype=cellTypes$Subtype[x]    
            }
        }
    }
    #### TO DO #####
    #if(is.null(cellType)){
        # check if combo is in conditional combinations
    #    for(cct in conditionalTypes){
            
    #    } 
    #}
    ###############
    if(is.null(cellType)){
        return(list(cellType="Unreasonable",subtype=""))
    } else {
        return(list(cellType=cellType,subtype=subtype))
    }

}

#' Flatten meta data for a single study
#'
#' Given all available meta data files for a single study, flatten all data into three
#' tables: 
#'   (1) SlideAnnotation - all information regarding samples and FOV
#'   (2) MarkerAnnotation - all information describing individual markers
#'   (3) MarkerCombinationAnnotation - all information describing important marker
#'       combinations
#'
#' @param metaDir            directory containing all available meta data files in XLSX format; default=NULL
#' @param metaFiles          vector of meta data files; an alternative to providing metaDir; default=NULL
#' @param sampAnnFile        sample annotation file in XLSX format; if given, this value will override a sample
#'                           annotation file found in metaDir or metaFiles; default=NULL 
#' @param fovAnnFile         FOV annotation file in XLSX format; if given, this value will override a FOV
#'                           annotation file found in metaDir or metaFiles; default=NULL
#' @param cellTypesFile      XLSX file describing all known cell types (marker combinations); if given, this 
#'                           value will override a cell types file found in metaDir or metaFiles; default=NULL
#' @param haloPanelsFile     XLSX file describing all halo panels; if given, this value will override a halo
#'                           panels file found in metaDir or metaFiles; default=NULL
#' @param markerDescFile     XLSX file describing individual markers; if given, this value will override a 
#'                           marker description file found in metaDir or metaFiles; default=NULL
#' @param unreasonalbesFile  #### MOST LIKELY WILL BE REMOVED
#' @return a list containing three data frames: SlideAnnotation, MarkerAnnotation and MarkerCombinationAnnotation
#' @export
flattenMetaData <- function(metaDir=NULL,metaFiles=NULL, sampAnnFile=NULL, fovAnnFile=NULL,
                            cellTypesFile=NULL, haloPanelsFile=NULL, markerDescFile=NULL){ #,
                            #unreasonablesFile=NULL){
    
    if(is.null(metaFiles)){
        if(is.null(metaDir)){
            stop("Need either metaFiles or metaDir to NOT be NULL")
        }
        metaFiles <- file.path(metaDir, dir(metaDir)[grep(".xlsx",dir(metaDir))])
    }
    sampAnnFile        <- if(is.null(sampAnnFile)) { metaFiles[grep("_SampleAnnotations.xlsx",metaFiles)] }
    fovAnnFile         <- if(is.null(fovAnnFile)) { metaFiles[grep("_FOVannotations.xlsx",metaFiles)] }
    cellTypesFile      <- if(is.null(cellTypesFile)) { metaFiles[grep("_CellTypes.xlsx",metaFiles)] }
    haloPanelsFile     <- if(is.null(haloPanelsFile)) { metaFiles[grep("_HaloPanels.xlsx",metaFiles)] }
    markerDescFile     <- if(is.null(markerDescFile)) { metaFiles[grep("_MarkerDescriptions.xlsx",metaFiles)] }
    #unreasonablesFile  <- if(is.null(unreasonablesFile)) { metaFiles[grep("_UnreasonableCombinations_Matrix.xlsx",metaFiles)] }

    sampAnn          <- read.xlsx(sampAnnFile, 1, stringsAsFactors=FALSE)
    fovAnn           <- read.xlsx(fovAnnFile, 1, stringsAsFactors=FALSE)
    cellTypes        <- read.xlsx(cellTypesFile, sheetName="Simple",stringsAsFactors=FALSE)
    cellTypesComplex <- read.xlsx(cellTypesFile, sheetName="Complex",stringsAsFactors=FALSE)
    haloPanels       <- read.xlsx(haloPanelsFile, 1, stringsAsFactors=FALSE)
    markerDesc       <- read.xlsx(markerDescFile, 1, stringsAsFactors=FALSE)
    #unreasonables    <- read.xlsx(unreasonablesFile, 1, stringsAsFactors=FALSE)

    slideInfo       <- full_join(sampAnn, fovAnn, by=c("CELL_DIVE_ID"))
    indivMarkerInfo <- full_join(haloPanels, markerDesc, by=c("Marker_name"))
    markerCombos    <- getCellTypes(cellTypes, cellTypesComplex)#, unreasonables)  
 
    return(list(SlideAnnotation=slideInfo, MarkerAnnotation=indivMarkerInfo, MarkerCombinationAnnotation=markerCombos))
   
}

#' Generate a template marker configuration based on meta data
#' 
#' Given either a table of cell types or a cell types file to be parsed, generate in
#' list format a configuration of marker combinations to be used for downstream 
#' analysis. The purpose of this template is to provide the beginnings of any and all more complex
#' configuration needed to calculate statistics and create plots of halo data
#'
#' @param allCellTypes    flattened table containing all cell type file found in meta data files; if NULL,
#'                        user must provide cellTypesFile to be parsed; default=NULL
#' @param cellTypesFile   cell types file in XLSX format structured according to meta data documentation rules;
#'                        if NULL, user must provide a pre-generated table consisting of the data found
#'                        in meta data file
#' @param writeYAML       logical indicating whether to write configuration to file; default=TRUE
#' @param outDir          directory where YAML file should be written; default=getwd()    
#' @return all marker information and default marker configuration values
#' @export
getTemplateMarkerConfig <- function(allCellTypes=NULL, cellTypesFile=NULL, writeYAML=TRUE, outDir=getwd()){
    if(is.null(allCellTypes)){
print("allCellTypes is null")
        cellTypes <- read.xlsx(cellTypesFile, sheetName="Simple",stringsAsFactors=FALSE)
        cellTypesComplex <- read.xlsx(cellTypesFile, sheetName="Complex",stringsAsFactors=FALSE)   
        #unreasonables    <- read.xlsx(unreasonablesFile, 1, stringsAsFactors=FALSE)
        allCellTypes <- getCellTypes(cellTypes, cellTypesComplex)#, unreasonables)
    } 

    mCfg <- list()

    mCfg$all_cell_type_markers <- gsub("-","",unique(unlist(strsplit(allCellTypes$Marker_combination,","))))
    mCfg$cell_type_marker_combinations <- unique(allCellTypes$Marker_combination)
    mCfg$marker_sets <- list()

    for(ct in unique(allCellTypes$Cell_type)){
        ctName <- gsub("\\s+","_",ct)
        mCfg$marker_sets[[ctName]] <- list()
        mCfg$marker_sets[[ctName]][["populations"]] <- unique(allCellTypes$Marker_combination[which(allCellTypes$Cell_type == ct)])
        mCfg$marker_sets[[ctName]][["functional"]] <- NULL
        mCfg$marker_sets[[ctName]][["cell_type_labels"]] <- list()
        mCfg$marker_sets[[ctName]][["cell_type_labels"]][["remove_negative_markers"]] <- TRUE
        mCfg$marker_sets[[ctName]][["cell_type_labels"]][["remove_populations_markers"]] <- FALSE
    }

    if(writeYAML){
        write(as.yaml(mCfg, indent=4, indent.mapping.sequence=TRUE), file=file.path(outDir, "template_marker_config.yaml"))
    } 

    return(mCfg)
}

#' Generate a template plot configuration based on marker configuration
#' 
#' Given either a marker configuration in list format, or a configuration file in YAML
#' format, create a template plot configuration to be modified manually by user for plotting
#'
#' @param markerConfig     marker configuration in nested list format (most likely previously-read
#'                         from a YAML file); if NULL, user must provide a marker config file in YAML
#'                         format; default=NULL
#' @param markerConfigFile marker configuration YAML file; if NULL, user must provide configuration
#'                         in list format; default=NULL
#' @param writeYAML        logical indicating whether to write configuration to file; default=TRUE
#' @param outDir           directory where YAML file should be written; default=getwd()    
#' @return all plot information and default plot configuration values
#' @export
getTemplatePlotConfig <- function(markerConfig, writeYAML=TRUE, outDir=getwd()){

    mCfg <- markerConfig

    pCfg <- list()

    pCfg[["plot_options"]] <- list()
    pCfg$plot_options[["marker_colors"]] <- list()
    pCfg$plot_options$marker_colors[["Negative"]] <- "#caced1"

    allCombos <- c()

    for(ms in names(mCfg$marker_sets)){
        msCombos <- c()
        for(p in mCfg$marker_sets[[ms]][["populations"]]){

            msCombos <- c(msCombos,mCfg$marker_sets[[ms]][["populations"]])

            if("functional" %in% names(mCfg$marker_sets[[ms]])){
                for(f in mCfg$marker_sets[[ms]][["functional"]]){
                    if("functional_negatives" %in% names(mCfg$marker_sets[[ms]])){
                        all_func <- unlist(strsplit(mCfg$marker_sets[[ms]][["functional"]], ","))
                        f <- addNegatives(f, all_func)
                    }
                    if("remove_population_markers" %in% names(mCfg$marker_sets[[ms]][["cell_type_labels"]])){
                        msCombos <- c(msCombos, f)
                    } else {
                        msCombos <- c(msCombos, paste0(p,",",f))
                    }
                }

            }
        }  

        if("cell_type_labels" %in% names(mCfg$marker_sets[[ms]])){
            if("remove_negative_markers" %in% names(mCfg$marker_sets[[ms]][["cell_type_labels"]])){
                for(x in seq(msCombos)){
                    msCombos[x] <- removeNegativeMarkers(msCombos[x])
                }
            }
        }   

        allCombos <- c(allCombos, msCombos)       
    }

    for(c in unique(allCombos)){
        if(!c %in% names(pCfg$plot_options$marker_colors)){
            pCfg$plot_options$marker_colors[[c]] <- "gray"
        }
    }

    if(writeYAML){
        write(as.yaml(pCfg, indent=4, indent.mapping.sequence=TRUE), file=file.path(outDir, "template_plot_config.yaml"))
    } 

    return(pCfg)
}

#' Create a study config template
#' 
#' Create a study config template, returning config in list format and optionally writing
#' config to yaml file. User may choose to set all parameters to NULL or to set any available
#' default values.
#'
#' @param studyName
#' @param setDefaultDirectoryStruct 
getTemplateStudyConfig <- function(studyDirectory=getwd(), setDefaultDirectoryStruct=TRUE, studyName=NULL, writeYAML=TRUE){
     ## set defaults
    sCfg <- list(raw_data_dir                = NULL,
                 raw_data_files              = NULL,
                 data_dir                    = NULL,
                 data_files                  = NULL,
                 meta_dir                    = NULL,
                 meta_files                  = NULL,
                 drift_dir                   = NULL,
                 drift_files                 = NULL,
                 study_dir                   = NULL, 
                 log                         = NULL,
                 marker_config_file          = NULL,
                 plot_config_file            = NULL,
                 raw_marker_combo_table_file = NULL,
                 celltype_config_file        = NULL,
                 marker_analysis_config_file = NULL,
                 pad                         = 20,
                 drift_threshold             = 0.1,
                 max_g                       = 5,
                 band_width                  = 10,
                 by_band                     = TRUE,
                 max_distance_from_interface = 360,
                 debug                       = TRUE,
                 debug_dir                   = NULL,
                 infiltration_dir            = NULL,
                 infiltration_density_dir    = NULL,
                 infiltration_density_file   = NULL,
                 infiltration_area_dir       = NULL,
                 infiltration_area_file         = NULL,
                 infiltration_band_assignments_file = NULL,
                 fov_stats_dir               = NULL,
                 fov_density_dir             = NULL,
                 fov_density_file            = NULL,
                 fov_area_dir                = NULL,
                 fov_area_file               = NULL)

    if(setDefaultDirectoryStruct){
        if(is.null(studyName)){
            stop("Need studyName in order to set default directory structure.")
        }
        sCfg$study_dir                <- file.path(studyDirectory, studyName)
        sCfg$config_dir               <- file.path(sCfg$study_dir, "config")
        sCfg$debug_dir                <- file.path(sCfg$study_dir, "debug")
        sCfg$data_dir                 <- file.path(sCfg$study_dir, "objectAnalysisData")
        sCfg$infiltration_dir         <- file.path(sCfg$study_dir, "infiltration")
        sCfg$infiltration_density_dir <- file.path(sCfg$infiltration_dir, "density")
        sCfg$infiltration_area_dir    <- file.path(sCfg$infiltration_dir, "area")
        sCfg$fov_stats_dir            <- file.path(sCfg$study_dir, "fov_data")
        sCfg$fov_density_dir          <- file.path(sCfg$fov_stats_dir, "density")
        sCfg$fov_area_dir             <- file.path(sCfg$fov_stats_dir, "area")
        sCfg$log                      <- file.path(sCfg$study_dir, paste0(studyName,".log"))
    }

    if(writeYAML){
        write(as.yaml(sCfg, indent=4, indent.mapping.sequence=TRUE), file=file.path(sCfg$config_dir, "study_config.yaml"))
    }
    return(sCfg)
}

#' Configure study
#'
#' Set study parameters and write YAML file to be used as study config
#'
#' @param studyName
#' @param studyDir
#' @param setDefaultDirectoryStructure
#' @param configDir
#' @param studyConfigFile
#' @param raw_data_dir
#' @param raw_data_files
#' @param data_dir
#' @param data_files
#' @param meta_dir
#' @param meta_files
#' @param annotations_dirs
#' @param annotations_files
#' @param drift_dir
#' @param drift_files
#' @param study_dir
#' @param cohort_dir
#' @param infiltration_dir
#' @param fov_stats_dir
#' @param infiltration_density_dir
#' @param infiltration_area_dir
#' @param infiltration_density_file
#' @param infiltration_area_file
#' @param fov_density_dir
#' @param fov_density_file
#' @param fov_area_dir
#' @param fov_area_file
#' @param marker_config_file
#' @param plot_config_file
#' @param marker_analysis_config_file
#' @param celltype_config_file
#' @param log
#' @param debug
#' @param raw_marker_combo_table_file
#' @param pad
#' @param drift_threshold
#' @param updateExistingConfigFiles
#' @param write_csv_files
#' @param max_g
#' @param band_width
#' @param maximum_distance_from_interface
#' @return nothing 
#' @export
configureStudy <- function(studyName=NULL, studyDir=NULL, setDefaultDirectoryStruct=TRUE, configDir=NULL, 
                           studyConfigFile=NULL, raw_data_dir=NULL, raw_data_files=NULL, 
                           data_dir=NULL, data_files=NULL, meta_dir=NULL, meta_files=NULL, 
                           annotations_dirs=NULL, annotations_files=NULL,
                           drift_dir=NULL, drift_files=NULL, study_dir=NULL, cohort_dir=NULL, 
                           infiltration_dir=NULL, fov_stats_dir=NULL, infiltration_density_dir=NULL, 
                           infiltration_density_file=NULL, infiltration_area_dir=NULL, infiltration_band_assignments_file=NULL,
                           infiltration_area_file=NULL, fov_density_dir=NULL, fov_density_file=NULL,
                           fov_area_dir=NULL, fov_area_file=NULL, marker_config_file=NULL,
                           plot_config_file=NULL, marker_analysis_config_file=NULL, celltype_config_file=NULL, log=NULL, debug=TRUE, raw_marker_combo_table_file=NULL,
                           pad=20, drift_threshold=0.1, updateExistingConfigFiles=TRUE, write_csv_files=TRUE,
                           max_g=5, band_width=10, maximumDistanceFromInterface=360){

    allConfig <- c("studyName", "studyDir", "configDir", "raw_data_dir", "raw_data_files",
                   "data_dir", "data_files", "meta_dir", "meta_files", "drift_dir", "drift_files",
                   "annotations_dirs", "annotations_files", "write_csv_files",
                   "study_dir", "cohort_dir", "infiltration_dir", "fov_stats_dir", "infiltration_density_dir", 
                   "infiltration_area_dir", "infiltration_density_file", "infiltration_area_file", "infiltration_band_assignments_file",
                   "fov_density_dir", "fov_density_file", "fov_area_dir", "fov_area_file", 
                   "log", "debug", "pad", "drift_threshold", "plot_config_file", "marker_config_file",
                   "celltype_config_file", "marker_analysis_config_file", "raw_marker_combo_table_file",
                   "max_g", "band_width", "maximum_distance_from_interface")

    if(is.null(studyDir)){ studyDir <- getwd() }

    ## create template
    print("Creating template study config")
    tmpCfg <- getTemplateStudyConfig(studyName=studyName, studyDir=studyDir,
                                     setDefaultDirectoryStruct=setDefaultDirectoryStruct,
                                     writeYAML=FALSE) 
    if(!is.null(studyConfigFile)){
        print("Reading study config")
        sCfg <- read_yaml(studyConfigFile)

        ## if a default setting exists that is not in the given config at all, 
        ## add it. if it does exist in the given config and it is a directory setting, 
        ## if setDefDirStruct==TRUE and the *_files setting that goes with the dir setting
        ## is NULL, set the default 
        for(c in names(tmpCfg)){
            if(!c %in% names(sCfg)){
                sCfg[[c]] <- tmpCfg[[c]]
            } else {
                if(grepl("_dir",c) && is.null(sCfg[[gsub("_dir","_files",c)]]) && setDefaultDirectoryStruct){
                    sCfg[[c]] <- tmpCfg[[c]]
                }
            }
        }
    } else {
        sCfg <- tmpCfg
    }
    ## set any other parameters passed to this function
    for(param in allConfig){
        if(exists(param) && !is.null(get(param))){
print(paste0("setting param ",param," to ",get(param)))
            sCfg[[param]] <- get(param)
        } 
    }

    ## create directories
    for(d in names(sCfg)[grep("_dir",names(sCfg))]){
print(d)
        if(!is.null(sCfg[[d]])){
print(paste0("Creating dir: ",sCfg[[d]]))
            dir.create(sCfg[[d]], recursive=T, showWarnings=T)
        }
    }

    if(("meta_dir" %in% names(sCfg) && !is.null(sCfg$meta_dir) && sCfg$meta_dir != "") | 
       ("meta_files" %in% names(sCfg) && !is.null(sCfg$meta_files) && length(sCfg$meta_files) > 0)){
       print("flattening meta data")
       flatMeta <- flattenMetaData(metaDir=sCfg$meta_dir, metaFiles=sCfg$meta_files)
       print("getting marker config")
       tmpMarkerCfg <- getTemplateMarkerConfig(allCellTypes=flatMeta$MarkerCombinationAnnotation$CellTypes,
                                               writeYAML=F, outDir=sCfg$config_dir)
       if(!is.null(sCfg$marker_config_file) && file.exists(sCfg$marker_config_file)){
           mCfg <- read_yaml(sCfg$marker_config_file)
           if(updateExistingConfigFiles && !all(names(tmpMarkerCfg) %in% names(mCfg))){
               mCfg <- c(mCfg, tmpMarkerCfg[-which(names(tmpMarkerCfg) %in% names(mCfg))])
               for(p in names(mCfg)){
                   if(is.null(mCfg[[p]]) && p %in% names(tmpMarkerCfg)){
                       mCfg[[p]] <- tmpMarkerCfg[[p]]
                   }
               }
               print(paste0("WARNING: Overwriting existing marker config file ",sCfg$marker_config_file))
               write(as.yaml(mCfg, indent=4, indent.mapping.sequence=TRUE), file=sCfg$marker_config_file)
           }
       } else {
           ## write marker config template if one does not exist
           print("Writing marker config template")
           mCfg <- tmpMarkerCfg
           sCfg$marker_config_file <- file.path(sCfg$config_dir, "template_marker_config.yaml") 
           write(as.yaml(mCfg, indent=4, indent.mapping.sequence=TRUE), file=sCfg$marker_config_file)
       }

       print("setting up plot config")
       tmpPlotCfg <- getTemplatePlotConfig(mCfg, writeYAML=F, outDir=sCfg$config_dir)
       if(!is.null(sCfg$plot_config_file) & file.exists(sCfg$plot_config_file)){
           pCfg <- read_yaml(sCfg$plot_config_file)
           if(updateExistingConfigFiles && !all(names(tmpPlotCfg) %in% names(pCfg))){
               pCfg <- c(pCfg, tmpPlotCfg[-which(names(tmpPlotCfg) %in% names(pCfg))])
               for(p in names(pCfg)){
                   if(is.null(pCfg[[p]]) && p %in% names(tmpPlotCfg)){
                       pCfg[[p]] <- tmpPlotCfg[[p]]
                   }
               }
               print(paste0("WARNING: Overwriting existing plot config file ",sCfg$plot_config_file))
               write(as.yaml(pCfg, indent=4, indent.mapping.sequence=TRUE), file=sCfg$plot_config_file)
           }
       } else {
           ## write plot config template if one does not exist
           print("Writing plot config template")
           pCfg <- tmpPlotCfg 
           sCfg$plot_config_file <- file.path(sCfg$config_dir, "template_plot_config.yaml") 
           write(as.yaml(pCfg, indent=4, indent.mapping.sequence=TRUE), file=sCfg$plot_config_file)
       }
    } else {
        warning("No meta dir or files given. Could not create marker or plot configs")
    }

    ## validate config

    ## write study config
    write(as.yaml(sCfg, indent=4, indent.mapping.sequence=TRUE), file=file.path(sCfg$config_dir, "study_config.yaml"))

}


