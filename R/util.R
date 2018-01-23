#' Log message to file and screen if specified
#' 
#' Print date() and log message to file, and if
#' verbose is set to TRUE, also print to screen
#' 
#' @param msg        message to be printed
#' @param v          verbose - when set to TRUE (default), messages
#'                   will be printed to screen as well as to file
#' @param logFile    file to write log messages to; DEFAULT=NULL
#' @param type      message type (INFO, DEBUG, WARN); DEFAULT=INFO
#' @export
logMsg <- function(msg, v=TRUE, logFile=NULL, type="INFO"){
    ## write log message to file, and if verbose, print
    ## to stdout
    if(v){
        write(paste0("[",type,"]","\t",date(),"\t",msg),stdout())
    }
    if(!is.null(logFile)){
        write(paste0("[",type,"]","\t",date(),"\t",msg), logFile, append=TRUE)
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
               data_dir=NULL,
               cancer=NULL,
               sample=NULL,
               pad=NULL,
               fov=NULL,
               log=NULL,
               verbose=TRUE,
               alt_bases=c(),
               run_counts=TRUE,
               run_frac_total=TRUE,
               run_medians=TRUE,
               pie_charts=FALSE,
               counts_rda_file=NULL,
               counts_xlsx_file=NULL,
               cell_type_markers=NULL,
               exclude_sample_fov=NULL)

    man <- read.delim(file, sep="\t", comment.char="#", header=FALSE, stringsAsFactors=FALSE)[,c(1,2)]
    for(x in 1:nrow(man)){
        if(length(grep(",",man[x,2])) > 0){
            pp[[man[x,1]]] <- trimws(unlist(strsplit(man[x,2],",")))
        } else if(man[x,2] %in% c("TRUE","FALSE")){
            pp[[man[x,1]]] <- ifelse(man[x,2] == "TRUE",TRUE,FALSE)
        } else if(man[x,2] == "NULL"){
            pp[[man[x,1]]] <- NULL
        } else {
            pp[[man[x,1]]] <- man[x,2]
        }
    }
    if(is.null(pp$log)){
        df <- file.path(pp$data_dir,dir(pp$data_dir)[grep("\\.rda$",dir(pp$data_dir))]) 
        pp$log <- projectFileName(pp$markers,df,pp$pad,"log")
    }
    return(pp)
}

#' Remove from counts tibble any FOV to be excluded
#'
#' Filter data to exclude specific FOV for specific samples
#' 
#' @param dat         counts tibble from countMarkers()
#' @param exclusions  a string of exlusions in the form: Sample1:3+5+9,Sample2:1+16+22
#' @param v           verbose - when set to TRUE, messages
#'                              will be printed to screen as well as to file;
#'                              DEFAULT = TRUE
#' @param debug       print debug messages; default=FALSE
#' @export
remove_exclusions <- function(dat,exclusions,v=TRUE,debug=FALSE){
    exclude_sample_fov <- trimws(unlist(strsplit(exclusions,",")))
    for(ex in exclude_sample_fov){
        samp <- unlist(strsplit(ex,":"))[1]
        fovEx <- unlist(strsplit(ex,":"))[2]
        fovEx <- unlist(strsplit(fovEx,"\\+"))
        if(debug){ logMsg(paste0("Removing ",samp," FOV ",fovEx)) }
        dat <- filter(dat, !(Sample == samp & FOV %in% fovEx))
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
#' @param type          file type ("rda","txt","xlsx"); essentially
#'                      the file extension
#' @return  file name to be used for saving counts
#' @export
projectFileName <- function(markerFile,dataFiles,pad,type){
    if(length(dataFiles)==1) {
        outFile <- paste("markerTable",
            gsub("_MegaTableV2.rda","",basename(dataFiles[1])),
            file_path_sans_ext(basename(markerFile)),sep="_")
        if(pad > 0){
            outFile <- paste0(outFile, "__PAD_",pad)
        }
    } else {
        outFile <- paste("markerTable",
            basename(dirname(dataFiles[1])),
            substr(digest(sort(dataFiles)),1,8),"__",
            file_path_sans_ext(basename(markerFile)),sep="_")
    }
    if(pad > 0){
        outFile <- paste0(outFile,"__PAD_",pad)
    }
    outFile <- paste0(outFile,".",type)
    return(outFile)
}

