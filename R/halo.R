#' Log message to file and screen if specified
#' 
#' Print date() and log message to file, and if
#' verbose is set to TRUE, also print to screen
#' 
#' @param msg        message to be printed
#' @param v          verbose - when set to TRUE (default), messages
#'                   will be printed to screen as well as to file
#' @param logFile    file to write log messages to; DEFAULT=NULL
#' @export
logMsg <- function(msg, v=TRUE, logFile=NULL){
    ## write log message to file, and if verbose, print
    ## to stdout
    if(v){
        write(paste("[INFO]",date(),msg,sep="\t"),stdout())
    }
    if(!is.null(logFile)){
        write(msg, logFile, append=TRUE)
    }
}

#' Build list of project parameters and values by reading
#' project manifest
#' 
#' Given a tab-delimited file of key-value pairs, read
#' them into list. File may contain comments both on top and in columns
#' beyond the second, but will be ignored here
#'
#' TO DO: Add to description list of required and optional keys in file
#' 
#' @param    file    manifest file to be read
#' @return   a list of all project parameters
#' @export
projectParams <- function(file){
    pp <- list(markers=NULL,
               dataDir=NULL,
               cancer=NULL, 
               sample=NULL, 
               pad=NULL,
               FOV=NULL,
               log=NULL,
               verbose=TRUE,
               altBases=c())

    man <- read.delim(file, comment.char="#", header=FALSE)[,c(1,2)]
    for(x in nrow(man)){
        if("," %in% man[x,2]){
            pp[[man[x,1]]] <- trimws(unlist(strsplit(man[x,2])))
        } else {
            pp[[man[x,1]]] <- man[x,2]
        }
    }
    return(pp)
}

#' Load data from *.rda file
#'
#' Load data for one subsample from one cancer type, and optionally
#' one FOV. Default is to load ALL FOV for a subsample 
#' 
#' @param    cancerType    cancer type (default: melanoma)
#' @param    subSample     subsample name, as it is called in *.rda file
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

#' Split data in *.rda file by FOV (SPOT)
#' 
#' Take in *.rda file containing data for multiple
#' FOV and output one *.rda file for each one. Input
#' tibble must have FOV column named "SPOT". Output
#' file will be named "FOV_[$SPOT].rda"
#' 
#' @param dat      input tibble/data.frame
#' @param outDir   output directory
#' @export
split_by_fov <- function(dat,outDir){
    for(i in levels(as.factor(dat$SPOT))){
        fname <- file.path(outDir,paste0("FOV_",i))
        fov <- filter(dat, SPOT==i)
        saveRDS(fov, paste0(fname,".rda"))
    }
}

#' Compare list of unique markers in a data set with the list
#' of markers in a marker file
#' 
#' Print a warning with a list of any markers that are in m1 but
#' NOT in m2
#'
#' @param dat       tibble containing marker data
#' @param m2        list of individual markers in marker file
#' @param v         verbose - when set to TRUE, messages
#'                  will be printed to screen as well as to file;
#'                  DEFAULT = TRUE
#' @param logFile   file to write log messages to; DEFAULT=NULL
#' @return  a tibble with extra markers filtered out
#' @export
validateMarkers <- function(dat,m2,v=FALSE,logFile=NULL){
    datMarkers <- unique(dat$MarkerName)
    finalMarkers <- intersect(datMarkers, m2)
    if(length(setdiff(datMarkers,m2)) > 0){
        msg <- paste0("WARNING: Markers in data but NOT in markerFile: ",paste(setdiff(datMarkers,m2),collapse=", "))
        if(!is.null(logFile)){
            logMsg(msg,v=FALSE,logFile)
        } else {
            warning("Markers in data but NOT in markerFile: ",paste(setdiff(datMarkers,m2),collapse=", "))
        }
    }
    if(length(setdiff(m2, datMarkers)) > 0){
        msg <- paste0("WARNING: Markers in markerFile but NOT in data: ",paste(setdiff(m2,datMarkers),collapse=", "))
        if(!is.null(logFile)){
            logMsg(msg,v=FALSE,logFile)
        } else {
            warning("Markers in markerFile but NOT in data: ",paste(setdiff(m2,datMarkers),collapse=", "))
        }
    }
    return(filter(dat, MarkerName %in% finalMarkers))
}

#' Make sure list of unique FOVs contains only one element
#'
#' If the list is longer than one element, print error message and EXIT 
#'
#' @param uniqFOVs    vector containing list of unique FOVs in dataset
#' @param v           verbose - when set to TRUE, messages
#'                    will be printed to screen as well as to file;
#'                    DEFAULT = TRUE
#' @param logFile     file to write log messages to; DEFAULT=NULL
#' @export
validateSingleFOVFile <- function(uniqFOVs, v=FALSE, logFile=NULL){
    ## ensure that file contains only ONE FOV
    if(length(uniqFOVs) > 1){
        logMsg("ERROR: Multiple FOVs found in file.",v,logFile) 
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
    pixel2um <- 0.293
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


#' Generate a XLSX file of counts for all markers in a given file
#' 
#' Given a marker file and a directory of *.rda files each containing Halo
#' data for a single sample, generate a XLSX file with multiple sheets of
#' counts, median counts for each marker in each sample, and fractions of counts
#' using alternate "bases" (e.g., DAPI)
#' 
#' By default, the file will contain at minimum a sheet of counts and a sheet
#' containing fractions of total counts. In all sheets, a column represents a marker 
#' and a row represents a single FOV from a single sample. If a vector of markers
#' is provided (e.g., c("DAPI","CD3")), another sheet will be created for each of
#' them containing the fraction of cells using those markers as "bases" 
#'
#' @param markerFile    File containing list of marker names
#' @param dataFiles     vector of data file(s)
#' @param pad           amount that will be trimmed from FOV
#' @param v             verbose - when set to TRUE, messages
#'                      will be printed to screen as well as to file;
#'                      DEFAULT = TRUE
#' @param logFile       file to write log messages to; DEFAULT=NULL
#' @param outFile       name of XLSX file to write to; if NULL, name will
#'                      be automatically generated according to input file names and padding
#' @param runCounts     include counts sheet in output; DEFAULT=TRUE
#' @param runFracTotal  include fractions of total counts sheet in output; DEFAULT=TRUE
#' @param altBases      vector of additional markers for which fractions 
#'                      of counts should be calculated; one sheet will be generated
#'                      for each element 
countMarkers <- function(markerFile, dataDir, lf=NULL, v=TRUE, pad=0, outFile=NULL,
               runCounts=TRUE, runFracTotal=TRUE, altBases=c()){
    ############
    ### TO DO: change this function to take in cancer type, subsample and FOV instead of 
    ### data directory so that we can just load data from within package????
    ###########

    markerNames <- scan(markerFile,"")
    dataFiles   <- file.path(dataDir,dir(dataDir)[grep("\\.rda$",dir(dataDir))])
    outFile     <- markerStatsFile(markerFile,dataFiles,pad)

    ## initialize tables that will be written to separate xlsx sheets
    allTbls     <- list()
    sheetOrder <- c("Counts","Median.Counts")
    for(s in c("TotalCells",altBases)){ sheetOrder <- c(sheetOrder, s, paste0("Median.Frac.",s)) }
    for(tableName in sheetOrder){ allTbls[[tableName]] <- tibble() }

    for(f in dataFiles){
        fov <- file_path_sans_ext(basename(f))
        logMsg(fov,v,lf)
        dat <- readRDS(f)

        logMsg(paste0("    Validating file ",f),v,lf)

        logMsg(paste0("Trimming ", pad, " pixels from image"),v,lf)
        dat <- trimImage(dat,pad)

        ## add new column containing marker names that alone will
        ## indicate whether a cell is positive or negative for that
        ## marker (e.g., SOXP10 with Value == 0 becomes SOXP10-)
        dat <- dat %>%
            filter(ValueType == "Positive") %>%
            select(UUID,SubSample,SPOT,SLICE,Marker,Value) %>%
            mutate(Sign = ifelse(Value == 1, "", "-")) %>%
            unite(MarkerName, Marker, Sign, sep="")

        logMsg(paste0("    Validating markers ",markerFile),v,lf)
        dat <- validateMarkers(dat, unique(unlist(strsplit(markerNames,","))),v,lf)
        dat$SPOT <- as.numeric(dat$SPOT)

        ## change all values to 1, as any nonsensical relationships
        ## will be changed to NA and then 0 by spread(), with fill=0 and drop=F.
        ## ultimately 1 will indicate positive for every marker (including
        ## negative markers)
        ## e.g., for every zero in the original data, a marker name with "-"
        ## and a value of 1 will be in the new data
        ## 
        ## create a single column for each Marker (negative and positive) with
        ## 1 indicating it should be counted and 0 indicating it should NOT be
        ## counted
        dat$Value <- 1 
        dat <- spread(dat, MarkerName, Value,drop=FALSE,fill = 0)

        ## Some markers are eliminated completely if they have all zeros. Add them back
        for(x in setdiff(unlist(strsplit(markerNames,",")),names(dat))){
            warning(paste0(x, " is NOT in data set!!!!"),immediate.=FALSE)
            dat[[x]] = 0
        }

        logMsg("    Calculating stats",v,lf)
        for(x in markerNames){
            tmp <- select(dat, unlist(strsplit(x,",")))
            dat[[x]] <- ifelse(rowSums(tmp) == ncol(tmp),1,0)
        }
        tblRow <- select(dat,SubSample, SPOT, SLICE, markerNames) %>%
                 gather(markerNames, key=MarkerNames, value="tmp") %>% 
                 group_by(SubSample, SPOT, SLICE, MarkerNames) %>% 
                 summarise(Counts = sum(tmp)) %>%
                 spread(MarkerNames, Counts)
        allTbls[["Counts"]] <- bind_rows(allTbls[["Counts"]],tblRow)

        fracTotal <- tblRow
        fracTotal[1,4:ncol(fracTotal)] <- fracTotal[1,4:ncol(fracTotal)]/nrow(dat)        

        allTbls[["TotalCells"]] <- bind_rows(allTbls[["TotalCells"]],fracTotal)

        for(ab in altBases){
            tmp <- tblRow
            tmp[1,4:ncol(tmp)] <- tmp[,4:ncol(tmp)]/as.numeric(tmp[1,ab])
            allTbls[[ab]] <- bind_rows(allTbls[[ab]],tmp)
        } 
    }

    for(tblName in names(allTbls)){
        if(length(grep("Median",tblName)) > 0){
            row <- allTbls[[gsub("Median\\.|Frac\\.","",tblName)]] %>%
                   gather(markerNames, key=MarkerNames, value="tmp") %>%
                   group_by(SubSample, MarkerNames) %>%
                   summarise(Median = median(tmp)) %>%
                   spread(MarkerNames, Median)
            row$SPOT <- "NA"
            row$SLICE <- "NA"
            allTbls[[tblName]] <- bind_rows(allTbls[[tblName]],row) 
        }
    }

    logMsg("Writing output",v,lf)
    ## write output
    logMsg(paste0("Writing file ",outFile), v, lf)
    for(t in sheetOrder){
        append <- TRUE
        sheetName <- paste0("Frac.",t)
        if(t == "Counts"){
            append <- FALSE
            sheetName <- t
        }
        logMsg(paste0("Writing sheet ",sheetName), v, lf)
        tbl <- select(allTbls[[t]],SubSample,SPOT,SLICE,markerNames) %>%
                 arrange(SubSample,SPOT)
        writeXLSXsheet(tbl, outFile, sheetName = sheetName, append = append) 
    }
}
