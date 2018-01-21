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
    return(pp)
}

