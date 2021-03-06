options(stringsAsFactors=FALSE)

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


#' Make sure list of unique FOVs contains only one element
#'
#' If the list is longer than one element, print error message and EXIT 
#'
#' @param uniqFOVs    vector containing list of unique FOVs in dataset
#' @export
validateSingleFOVFile <- function(uniqFOVs){
    ## ensure that file contains only ONE FOV
    if(length(uniqFOVs) > 1){
        flog.fatal("Multiple FOVs found in this file.")
        stop("MULTIPLE FOVs found in file. Please split file by FOV and rerun")
    }
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


#' Shortcut to write a XLSX file 
#' 
#' Write a XLSX file with the following defaults:
#'    row.names = F
#'    append = TRUE
#'    sheetName = NULL
#'
#' @param dat        data to be written, either tibble, matrix or dataframe
#' @param outFile    name of output file
#' @param sheetName  name of Excel Worksheet
#' @param append     logical indicating whether to append to the file or
#'                   overwrite; Default=TRUE
writeXLSXsheet <- function(dat, outFile, sheetName=NULL, append=TRUE){
    write.xlsx(as.data.frame(dat), outFile, row.names=F, sheetName = sheetName, append = append)
}

#' Clean marker file
#' 
#' Remove spaces, new lines, tabs from file
#' 
#' @param markerFile   File condatining list of marker names
#' @param altBases     a vector of alternate markers that will be considered the baseline
#'                     for all other counts
#' @return  a vector of markers
cleanMarkers <- function(markerFile,altBases=NULL){
    markers <- scan(markerFile, "", sep="\n")
    markers <- gsub("[[:space:]]", "", markers)
    if(!is.null(altBases) & length(altBases) > 0){
        for(ab in altBases){
            if(ab != "DAPI"){
                m <- paste0(ab,",DAPI")
                if(! m %in% markers){ 
                    markers <- c(markers,m)
                }
            }
        }
    }
    return(markers)
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

#' Create median row
#' 
#' Given a tibble with the first three columns Sample, FOV, SLICE, 
#' caclulate medians of all remaining columns and return row to be added to tibble
#' 
#' @param  dat    table of counts where first three columns are Sample, 
#'                  FOV, SLICE, and the remaining columns are marker counts or fractions
#' @param  samp   sample name
#' @return   a row where Sample is the unique sample in input table, FOV and SLICE are NA,
#'           and the remaining values are medians of each column
medianRow <- function(dat,samp){
    meds <- as.list(apply(dat[,4:ncol(dat)],2,median))
    medRow <- list(Sample=samp, FOV=NA, SLICE=NA)   
    medRow <- c(medRow, meds)
    return(as.tibble(medRow))
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
#' @param runMedians	   include medians of counts for each marker in each sample; DEFAULT=TRUE
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

computeMultiMarkerTable<-function(mt,markerCols) {
    for(ii in which(grepl(",",markerCols))) {
        markers=strsplit(markerCols[ii],",")[[1]]
        mt[[markerCols[ii]]] = mt %>%
            dplyr::select(one_of(markers)) %>%
            mutate(XX=ifelse(rowSums(.)==len(markers),1,0)) %$%
            as.vector(XX)
    }

    ## Remove any DAPI negative cells
    ##mt %<>% filter(DAPI==1)
    mt
}
