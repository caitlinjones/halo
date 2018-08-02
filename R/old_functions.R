#' Load data from *.rda file
#'
#' Load data for one subsample from one cancer type, and optionally
#' one FOV. Default is to load ALL FOV for a subsample 
#' 
#' @param    cancerType    cancer type (default: melanoma)
#' @param    sample        sample name, as it is called in *.rda file
#' @param    FOV (SPOT)    integer; FOV/SPOT number (default: NULL)
#' @return   tibble saved in *.rda file
#' @export
loadHalo <- function(cancerType, sample, FOV=NULL){
    if(!is.null(FOV)){
        fileName <- paste0(cancerType, "/", sample, "_by_FOV/FOV_", FOV, ".rda")
    } else {
        fileName <- paste0(cancerType, "/", sample, "_ObjectAnalysisData_MegaTableV2.rda")
    }
    return(readRDS(system.file("extdata", fileName, package="halo")))
}

#' Remove from data markers that are not in marker file
#' 
#' Remove from data markers that are not in marker file
#' 
#' @param dat       tibble containing marker data
#' @param m2        list of individual markers in marker file
#' @return  a tibble with extra markers filtered out
#' @export
removeExtraMarkers <- function(dat, m2){
    datMarkers <- unique(dat$MarkerName)
    if(length(setdiff(datMarkers,m2)) > 0){
        flog.warn(paste0("Removing arkers in data but NOT in markerFile: ",
                          paste(setdiff(datMarkers,m2),collapse=", ")))
        warning(paste0("Removing arkers in data but NOT in markerFile: ",
                          paste(setdiff(datMarkers,m2),collapse=", ")))
        ## remove markers in data but not in marker file
        finalMarkers <- datMarkers[which(datMarkers %in% m2)]
        dat <- filter(dat, MarkerName %in% finalMarkers)
    }
    return(dat)
}

#' Add missing markers
#' 
#' Add to data markers that are in marker file but are missing
#' from given data set
#'
#' @param dat       tibble containing marker data
#' @param m2        list of individual markers in marker file
#' @return  a tibble with missing markers added, with "Positive" value
#'          set to NA
addMissingMarkers <- function(dat, m2){
    datMarkers <- unique(dat$MarkerName)
    if(length(setdiff(m2, datMarkers)) > 0){
        flog.warn(paste0("Markers in markerFile but NOT in data: ",
                   paste(setdiff(m2,datMarkers),collapse=", ")))
        warning(paste0("Markers in markerFile but NOT in data: ",
                   paste(setdiff(m2,datMarkers),collapse=", ")))
        ## add markers in marker file but not data
        ## if there is no pos or neg marker, set Value = NA; 
        ## otherwise, set to zero
        for(m in setdiff(m2, datMarkers)){
            tmp <- distinct(dat, UUID, Sample, FOV, SLICE)
            tmp$MarkerName <- m
            if(length(grep("-",m))==0){
                if(paste0(m,"-") %in% datMarkers){
                    flog.debug("%s is negative for ALL cells in ALL FOV. Setting counts to zero.",m)
                    tmp$Value <- 0
                } else{
                    flog.debug("%s is MISSING from data. Setting counts to NA.",m)
                    tmp$Value <- NA
                }
            } else {
                if(gsub("-","",m) %in% datMarkers){
                    flog.debug("%s is negative for ALL cells in ALL FOV. Setting counts to zero.",m)
                    tmp$Value <- 0
                } else{
                    flog.debug("%s is MISSING from data. Setting counts to NA.",m)
                    tmp$Value <- NA
                }
            }
            dat <- bind_rows(dat, tmp)
        }
    }
    return(dat)
}

#' Remove cells near edge of FOV from dataset
#' 
#' Given a tibble and an amount to trim around edges, remove cells
#' that fall outside the inner border
#'
#' @param dat    dataset (tibble)
#' @param pad    amount to trim around edges (pixels, double; will be converted
#'               to um)
#' @return tibble with cells that fall outside padding filtered out
#' @export
trimImage <- function(dat, pad=0){
    ## remove cells around edge of FOV (???)
    if(pad > 0){
        bb <- list(X0=min(dat$XMax),Y0=min(dat$YMax),
                   X1=max(dat$XMin),Y1=max(dat$YMin))
        padPx <- pad/pixel2um
        dat <- filter(dat,
            XMax > (bb$X0+padPx)
            & YMax > (bb$Y0+padPx)
            & XMin < (bb$X1-padPx)
            & YMin < (bb$Y1-padPx))
    }
    return(dat)
}

#' Generate name of output file for marker stats based on 
#' input file names
#' 
#' Use name of marker file and either a single data file name
#' or the directory of data files to generate output file name
#'
#' @param markerFile    File containing list of marker names
#' @param dataFiles     vector of data file(s)
#' @param pad           amount that will be trimmed from FOV
#' @return  file name to be used for saving counts
#' @export
markerStatsFile <- function(markerFile,dataFiles,pad){
    if(length(dataFiles)==1) {
        outFile <- paste("markerTable",
            gsub("_MegaTableV2.rda","",basename(dataFiles[1])),
            file_path_sans_ext(basename(markerFile)),
            ".xlsx",sep="_")
    } else {
        outFile <- paste("markerTable",
            basename(dirname(dataFiles[1])),
            substr(digest(sort(dataFiles)),1,8),"__",
            file_path_sans_ext(basename(markerFile)),
            ".xlsx",sep="_")
    }
    if(pad > 0){
        outFile <- gsub("\\.xlsx",paste0("__PAD_",pad,".xlsx"),outFile)
    }
    return(outFile)
}

#' Clean data table
#' 
#' Alter column names, ensure correct data types for certain columns
#' 
#' @param  dat    tibble containing marker data
#' @return   modified tibble
cleanData <- function(dat){
    flog.debug("changing SPOT to FOV and SubSample to Sample")
    dat <- rename(dat, FOV = SPOT)
    dat$FOV <- as.integer(dat$FOV)
    dat$SLICE <- as.integer(dat$SLICE)
    dat$Sample <- gsub("_ObjectAnalysisData","",dat$SubSample)
    return(dat)
}

#' Get order of workbook sheets for final excel file
#' 
#' Order sheets according to whether fractions and/or medians are to be run
#' 
#' @param  runCounts     include raw counts in final output; DEFAULT=TRUE
#' @param  runMedians    include a sheet of medians for each counts and fraction sheet; 
#'                       DEFAULT=TRUE
#' @param  runFracTotal  include a sheet of fractions of total cells; DEFAULT=FALSE
#' @param  altBases      include a sheet of fractions of cells positive for each of these 
#'                       markers; DEFAULT=NULL
#' @param  debug         print debug messages; DEFAULT=FALSE
getSheetOrder <- function(runCounts=TRUE, runMedians=TRUE, runFracTotal=FALSE, altBases=NULL,
                   debug=FALSE){
    sheetOrder <- c("Counts")
    if(runMedians){ sheetOrder <- c(sheetOrder, "Median.Counts") }
    if(runFracTotal == TRUE){ sheetOrder <- c(sheetOrder, "Frac.TotalCounts") }
    if(!is.null(altBases) & length(altBases) > 0){
        for(s in altBases){
            ## only count cells if DAPI is positive also 
            ## create mock marker "ab,DAPI"
            s = paste0("Frac.",s)
            if(s == "Frac.DAPI"){
                sheetOrder <- c(sheetOrder, s)
            } else {
                sheetOrder <- c(sheetOrder, paste(s,"DAPI",sep=","))
            }
            if(runMedians){
                sheetOrder <- c(sheetOrder,paste0("Median.",s))
            }
        }
    }
    return(sheetOrder)
}

#' Generate a XLSX file of counts for all markers in a given file
#' 
#' Given a marker file and a directory of *.rda files each containing Halo
#' data for a single sample, generate a XLSX file with multiple sheets of
#' counts, median counts for each marker in each sample, and fractions of counts
#' using alternate "bases" (e.g., DAPI)
#' 
#' By default, the file will contain at minimum a sheet of counts and a sheet
#' containing fractions of total counts. In all sheets, a column represents a marker 
#' and a row represents a single FOV from a single sample. 
#'
#' If a vector of markers is provided (e.g., c("DAPI","CD3")), another sheet will 
#' be created for each of them containing the fraction of cells using those markers as 
#' "bases". IMPORTANT NOTE: Any altBases in addition to DAPI must be counted ONLY IF DAPI
#' IS POSITIVE AS WELL (e.g., CD3+DAPI+). The names in the output .xlsx file, however
#' will exclude the DAPI+ to avoid redundancy. 
#'
#' @param markerFile       File containing list of marker names
#' @param dataDir          directory of *.rda files containing Halo data 
#' @param pad              amount that will be trimmed from FOV
#' @param countsXLSXFile   name of XLSX file to write to; if NULL, name will
#'                         be automatically generated according to input file names and padding
#' @param countsRDAFile    name of RDA file to write to; if NULL, RDA file will NOT be written
#' @param writeXLSXfile    write counts to XLSX file; DEFAULT=TRUE
#' @param saveRDSfile      save counts table to RDS file; DEFAULT=TRUE
#' @param runCounts        include counts sheet in output; DEFAULT=TRUE
#' @param runFracTotal     include fractions of total counts sheet in output; DEFAULT=FALSE
#' @param runMedians       include medians of counts for each marker in each sample; DEFAULT=TRUE
#' @param altBases         vector of additional markers for which fractions 
#'                         of counts should be calculated; one sheet will be generated
#'                         for each element 
#' @export
countMarkers <- function(markerFile, dataDir, pad=0, countsXLSXFile=NULL,
               countsRDAFile=NULL, writeXLSXfile=TRUE, saveRDSfile=TRUE, runCounts=TRUE,
               runFracTotal=FALSE, runMedians=TRUE, altBases=NULL){

    flog.debug("Reading marker file")
    markerNames <- cleanMarkers(markerFile,altBases)
    dataFiles <- file.path(dataDir,dir(dataDir)[grep("\\.rda$",dir(dataDir))])

    flog.info("Generating counts for %s and %s files in %s", basename(markerFile),
                  length(dataFiles), dataDir)

    if(writeXLSXfile & is.null(countsXLSXFile)){
        countsXLSXFile <- projectFileName(markerFile,dataFiles,pad,"xlsx")
    }
    if(saveRDSfile & is.null(countsRDAFile)){
        countsRDAFile <- projectFileName(markerFile,dataFiles,pad,"rda")
    }

    ###              
    ## initialize tables that will be written to separate xlsx sheets
    ###              
    flog.debug("Initializing tables")
    allTbls <- list()
    sheetOrder <- getSheetOrder(runCounts=runCounts, runFracTotal=runFracTotal,
                                runMedians=runMedians, altBases=altBases, debug=debug)
    for(tableName in sheetOrder){ allTbls[[tableName]] <- tibble() }

    ###              
    ## generate counts for each file and add them to the final tables
    ###              
    for(f in dataFiles){

        ## initialize tables for this one sample
        sampTbls <- list()
        for(tableName in sheetOrder){ sampTbls[[tableName]] <- tibble() }

        flog.info(file_path_sans_ext(basename(f)))
        flog.debug("Reading rda file")
        dat <- readRDS(f)

        flog.info("Trimming %s pixels from FOV",pad)
        dat <- trimImage(dat,as.numeric(pad))

        dat <- cleanData(dat)
        samp <- unique(dat$Sample)
        flog.debug("there is/are %s sample(s) in this file: ",samp)
        ## add new column containing marker names that alone will
        ## indicate whether a cell is positive or negative for that
        ## marker (e.g., SOX10 with Value == 0 becomes SOX10-)
        flog.debug("adding negative marker columns and adjusting pos/neg values")
        dat <- dat %>%
            filter(ValueType == "Positive") %>%
            dplyr::select(UUID,Sample,FOV,SLICE,Marker,Value) %>%
            mutate(Sign = ifelse(Value == 1, "", "-")) %>%
            unite(MarkerName, Marker, Sign, sep="")

        ## change all values to 1, as any nonsensical relationships
        ## will be changed to NA and then 0 by spread(); 1 will indicate positive, so if
        ## a marker with "-" has val 1, the pos marker is negative
        flog.debug("changing all values to 1 and spreading table, leaving NAs for nonsensical relationships")
        dat$Value <- 1

        ## sync markers in marker file and data
        dat <- removeExtraMarkers(dat, unique(unlist(strsplit(markerNames, ","))))
        dat <- addMissingMarkers(dat, unique(unlist(strsplit(markerNames, ","))))

        ###              
        ## process counts for each marker
        ###
        flog.info("    Processing counts for each marker")
        for(i in 1:length(markerNames)){
            x <- markerNames[i]
            flog.debug(paste0("    [",i,"] ",x))
            indivMarkers <- unlist(strsplit(x,","))

            tmp <- filter(dat, MarkerName %in% indivMarkers) %>%
                   spread(MarkerName, Value, drop = FALSE, fill = 0)
            tmp[[x]] <- ifelse(rowSums(tmp[,indivMarkers]) == length(indivMarkers),1,0)
            tmp <- dplyr::select(tmp, Sample, FOV, SLICE, x) %>%
                     gather(x, key=MarkerNames, value="tmp") %>%
                     group_by(Sample, FOV, SLICE, MarkerNames) %>%
                     summarize(Counts = sum(tmp)) %>%
                     spread(MarkerNames, Counts)
            if(length(sampTbls[["Counts"]]) == 0){
                sampTbls[["Counts"]] <- tmp
            } else {
                sampTbls[["Counts"]] <- left_join(sampTbls[["Counts"]], tmp, by=c("Sample","FOV","SLICE"))
            }
        }
        ###              
        ## calculate medians
        ## add to main table as well as a separate table of only medians
        ###              
        medRow <- medianRow(sampTbls[["Counts"]],samp)
        sampTbls[["Counts"]] <- bind_rows(sampTbls[["Counts"]],medRow)

        if(length(sampTbls[["Median.Counts"]]) == 0){
            sampTbls[["Median.Counts"]] <- medRow
        } else {
            sampTbls[["Median.Counts"]] <- bind_rows(sampTbls[["Median.Counts"]],medRow)
        }

        ###              
        ## calculate fractions
        ###       
        if(!is.null(altBases) & length(altBases) > 0){
            for(ab in altBases){
                m <- ifelse(ab == "DAPI", ab, paste0(ab,",DAPI"))
                tName <- paste0("Frac.",m)
                flog.debug(paste0("getting fraction of cells based on total ",ab," cells"))
                vals <- sampTbls[["Counts"]][[m]]
                tmp <- sampTbls[["Counts"]]
                tmp[,4:ncol(tmp)] <- tmp[,4:ncol(tmp)]/as.numeric(vals)
                sampTbls[[tName]] <- tmp

                medRow <- medianRow(sampTbls[[tName]],samp)
                medTbl <- paste0("Median.Frac.",ab)
                if(length(sampTbls[[medTbl]]) == 0){
                    sampTbls[[medTbl]] <- medRow
                } else {
                    sampTbls[[medTbl]] <- bind_rows(sampTbls[[medTbl]],medRow)
                }
            }
        }
        ###
        ## add to final tables
        ###
        for(s in sheetOrder){
            if(length(allTbls[[s]]) > 0){
                allTbls[[s]] <- bind_rows(allTbls[[s]], sampTbls[[s]])
            } else {
                allTbls[[s]] <- sampTbls[[s]]
            }
        }

        if(saveRDSfile){
            if(is.null(countsRDAFile)){
                countsRDAFile <- projectFileName(markerFile,dataFiles,pad,"rda")
            }
            saveRDS(allTbls,countsRDAFile)
        }
    }

    if(writeXLSXfile){
        if(is.null(countsXLSXFile)){
            countsXLSXFile <- projectFileName(markerFile,dataFiles,pad,"xlsx")
        }
        for(s in 1:length(sheetOrder)){
            append <- ifelse(s == 1,FALSE,TRUE)
            writeXLSXsheet(as.data.frame(allTbls[[sheetOrder[s]]]), countsXLSXFile,
                           sheetName=sheetOrder[s], append=TRUE)
        }
    }
    return(allTbls)

}

#' Remove exclusion boundaries that are contained in another one
#' 
#' Remove data from boundaries table that represent exclusion boundaries
#' that are completely surrounded by another exclusion boundary
#' 
#' @param boundaries            table generated by readHaloAnnotations
#'                             
cleanBoundaries_OLD <- function(boundaries){

    regionTable <- boundaries %>% bind_rows %>% count(RegionNum,RegionCode)
    excB <- boundaries[regionTable$RegionCode!="Tum"]
    tumB <- boundaries[regionTable$RegionCode=="Tum"]

    if(len(excB) > 1){
        ## make list of spatial polygons, one element for each exclusion boundary
        spExcB <- seq(excB) %>% map(function(b){boundaryToSpatialPolygon(excB[[b]],b)})

        ## are any contained within another?
        containedBoundary <- rep(FALSE,len(excB))
        for(i in seq(len(excB))){
            containedBoundary[i] <- spExcB[-i] %>% map(gContains,spExcB[[i]]) %>% unlist %>% any
        }
        ## remove those
        excB <- excB[!containedBoundary]
    }
    return(list(excB=excB,tumB=tumB))
}

#' Get a table of marker combination densities
#' 
#' Generate a table of cell type densities
#'
#' @param dat                       tibble containing halo data, specifically Sample, UUID, SPOT, Marker,
#'                                  Xmin/max, Ymin/max, ValueType, Value
#' @param dataFiles                 vector of *.rda files containing Halo data for a single sample
#' @param annotationsDirs           vector of directories containing *.annotations files, XML files of boundary
#'                                  annotations from Halo; order must correspond to order of dataFiles
#' @param pad                       amount to trim from each FOV
#' @param funcMarker                character string of functional marker to plot on top of every other marker
#' @param sortByMarker              sort FOVs by density of this marker (cell type)in ascending order;
#'                                  default=NULL
#' @param writeCSVfiles             logical; write CSV files, one of density values and one of area values; 
#'                                  default=TRUE
#' @param maxG                      maximum factor of 10 to use for generating random points on grid [ REWORD??? ]
#' @param outDir                    if writeCSVfiles=T, directory to which these files will be written
#' @param byBand                    logical indicating whether to break down density by distance from tumor interface; 
#'                                  default=FALSE
#' @param bandWidth                 width of each band around tumor interface, to be used when byBand is TRUE; 
#'                                  default=10
#' @param maxDistanceFromInterface  include only those cells within this distance from tumor interface
#' @return tibble containing density for all markers in all samples
#' @export
getMarkerDensityTable <- function(markers, dat=NULL, outDir=getwd(), dataFiles=NULL, annotationsDirs=NULL,
                                  haloAnnotations=NULL, areaFiles=NULL, densityFiles=NULL, sampleOrder=NULL,
                                  writeCSVfiles=FALSE, pad=20, byBand=TRUE, maxDistanceFromInterface=360,
                                  bandWidth=10, markerSetName=NULL, sortByMarker=NULL, funcMarker=NULL,
                                  maxG=5, calcFOVarea=FALSE){

    interfaceBins <- (-(maxDistanceFromInterface/bandWidth):(maxDistanceFromInterface/bandWidth))*bandWidth

    ia  <- NULL  ## interface area
    ba  <- NULL  ## band assignments
    rho <- NULL  ## density

    if(!is.null(densityFiles)){
        print("reading density files")
        ## read density files
        for(df in densityFiles){
            if(grepl("\\.rda$",df)){
                rho <- rho %>% bind_rows(readRDS(df))
            } else if(grepl("\\.csv$",df)){
                rho <- rho %>% bind_rows(read.delim(df,header=T,sep=","))
            }
        }
    } else {

        if(is.null(dat)){
            print("reading data files")
            for(df in dataFiles){
                ddat <- readRDS(df)
                ddat <- ddat %>% select(Sample, UUID, SPOT, Marker, matches("XM|YM"), ValueType, Value)
                dat <- dat %>% bind_rows(ddat)
            }
        }
        #dat$Sample <- gsub("_ObjectAnalysisData","",dat$Sample)

        if(!is.null(areaFiles)){
            print("reading area files")
            ## read area files
            for(af in areaFiles){
                if(grepl("\\.rda$",af)){
                    a <- readRDS(af)
                    ia <- ia %>% bind_rows(a$area)
                    ba <- ba %>% bind_rows(a$bandAssignments)
                } else if(grepl("\\.csv$",af)){
                    ia <- ia %>% read.delim(af,header=T,sep=",")
                    #ba <- TO DO: ADD THIS
                }
            }
            ia$Sample <- gsub("_ObjectAnalysisData","",ia$Sample)
            ba$Sample <- gsub("_ObjectAnalysisData","",ba$Sample)
        } else {
            print("reading data and annotations files")
            ## calculate area
            if(is.null(dat) | (is.null(annotationsDirs) && is.null(haloAnnotations))){
                stop(paste0("Need either density files, area files + dat, or dat +",
                            " annotations directories"))
            }

            ### calculate area then add it to the area tibble for entire analysis
            print(paste0("calculating area for all samples"))
            if(is.null(haloAnnotations)){
                aFiles <- file.path(annotationsDir,dir(annotationsDir)[grep("\\.annotations$",
                                                               dir(annotationsDir))])
            } else {
                aFiles <- NULL
            }
            flog.debug("calculating interface area")
            samp_ia <- calculateInterfaceArea(dat, haloAnnotations=haloAnnotations, aFiles=aFiles, writeCSVfiles=writeCSVfiles,
                                              maxG=maxG, outDir=outDir, statsByBand=byBand,
                                              interfaceBins=interfaceBins)

            if(!is.null(samp_ia$area) & !is.null(samp_ia$bandAssignments)){
                ia <- samp_ia$area
                ba <- samp_ia$bandAssignments
                saveRDS(samp_ia, file=file.path(outDir,"all_samples_area.rda"))
            }
        }
        ## calculate density
        if(!is.null(ia) & !is.null(ba)){
            for(s in unique(ia$Sample)){
                print(paste0("calculating density for ",s))
                sdat <- dat %>% filter(Sample == s)
                sia <- ia %>% filter(Sample == s)
                sba <- ba %>% filter(Sample == s)
                den <- calculateDensity(sdat, sia, sba,
                                    markers, pad, funcMarker=funcMarker,
                                    sortByMarker=sortByMarker,writeCSVfiles=writeCSVfiles,
                                    outDir=outDir, statsByBand=byBand)
                if(!is.null(den)){
                    if(byBand){
                        den <- den %>% dplyr::select(-one_of(names(den)[grep("Count",names(den))]))
                    } else {
                        den <- den %>% dplyr::select(Sample,SPOT,CellType,Total)
                    }
                    rho <- rho %>% bind_rows(den)
                    saveRDS(rho, file=file.path(outDir, "All_density.rda"))
                }
            }
        } else {
            print("WARNING: Interface area could not be calculated.")
            flog.warn("Interface area could not be calculated.")
        }
    }

    if(!"Density" %in% names(rho)){
        rho <- rho %>% dplyr::select(Density=Total, everything())
    }

    return(rho)
}

#' Summarize infiltration by getting total density of all FOV
#' in each band
#' 
#' Summarize infiltration by getting total density of all FOV
#' in each band
#'
#' @param markers    vector of markers for which density is to be calculated
#' @param ia         tibble of interface area; required if areaFiles not given; default=NULL
#' @param ba         tibble containing band assignments of markers; required if areaFiles 
#'                   not given; default=NULL
#' @param areaFiles  vector of files containing area and band assignments for each sample; 
#'                   required if ia or ba is not given; default=NULL
#' @return tibble containing total density for each band
#' @export
getInfiltrationDensitySummary <- function(markers, ia=NULL, ba=NULL, areaFiles=NULL){

    if(is.null(ia) | is.null(ba)){
        ## read areaFiles
        for(af in areaFiles){
            a <- readRDS(af)
            ia <- ia %>% bind_rows(a$area)
            ba <- ba %>% bind_rows(a$bandAssignments)
        }
    }
    ia$Sample <- gsub("_ObjectAnalysisData","",ia$Sample)
    ba$Sample <- gsub("_ObjectAnalysisData","",ba$Sample)
#    bandAreas <- ia %>% gather(3:ncol(ia), key="Band", value="Area")
    bandAreas <- ia

    markerCounts <-computeMultiMarkerTable(ba,markers) %>%
                   select(Sample,SPOT,Band,markers) %>%
                   filter(!is.na(Band))

    markerCounts <- markerCounts %>%
                    gather(4:ncol(markerCounts), key="CellType", value="Positive") %>%
                    group_by(Sample,SPOT,Band,CellType) %>%
                    summarise(Count=sum(Positive))

    summaryDat <- left_join(markerCounts, bandAreas) %>%
               group_by(Sample,Band,CellType) %>%
               summarise(TotalCountAllFOV=sum(Count), TotalAreaAllFOV=sum(Area)) %>%
               mutate(TotalDensity=TotalCountAllFOV/TotalAreaAllFOV)

    return(summaryDat)
}

getTotalSampleDensitySummary <- function(dat, markers, ia=NULL, areaFiles=NULL){
    if(is.null(ia)){
        if(is.null(areaFiles)){
            warning("No area data or area files provided. Can not calculate total FOV density summary")
            return(NULL)
        }
        for(af in areaFiles){
            a <- read.delim(af, header=T, sep=",")
            ia <- ia %>% bind_rows(a)
        }
    }

    markerCounts <- computeMultiMarkerTable(dat,markers)  %>%
                   select(Sample,SPOT,markers)

    markerCounts <- markerCounts %>%
                    gather(4:ncol(markerCounts), key="CellType", value="Positive") %>%
                    group_by(Sample,SPOT,CellType) %>%
                    summarise(Count=sum(Positive))

    summaryDat <- left_join(markerCounts, dat) %>%
               group_by(Sample,Band,CellType) %>%
               summarise(TotalCountAllFOV=sum(Count), TotalAreaAllFOV=sum(Area)) %>%
               mutate(TotalDensity=TotalCountAllFOV/TotalAreaAllFOV)

    return(summaryDat)
}


