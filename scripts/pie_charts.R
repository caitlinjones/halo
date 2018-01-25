#!/home/byrne/R/R-3.4.3/bin/R

options( java.parameters = c("-Xss2560k", "-Xmx8g") )
options('python_cmd'='/opt/common/CentOS_6-dev/python/python-3.5.1/bin/python')

write("Loading libraries...",stdout())
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("halo"))

###
# Parse user input
###
parser <- ArgumentParser()
## args are preferrably all in manifest
parser$add_argument("-m", "--manifest", type="character", default=NULL, help="file containing all project parameters; run ?projectParams for details")
## args required if manifest not given
parser$add_argument("--cell_type_markers", type="character", default=NULL, help="comma-separated string of markers, for each of which a pie chart will be generated; e.g., 'CD3,CD4,CD8,CD20'")
parser$add_argument("-c", "--counts_rda_file", type="character", help="*.rda file containing list of counts-related tables (output of counts.R)")
parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="print extra output")
parser$add_argument("-l", "--log", type="character", default=gsub(" ","_",date()), help="log file")
parser$add_argument("--debug", action="store_true", default=FALSE, help="print extra output for debugging")

args <- parser$parse_args()
####################################################


logParams <- function(x,action="Generating pie charts"){
    #for(n in names(x)){ if(is.null(x[[n]])){ x[[n]] <- NA } }
    write(date(),file=x$log)
    write(paste0("\n",action," with params:\n"),file=x$log,append=TRUE)
    for(n in names(x)){
        write(paste(n," = ",paste(x[[n]],collapse=",")), file=x$log, append=TRUE)
    }
    write("\n#######################################################################\n",file=x$log,append=TRUE)
}

pp <- NULL

## get all params from manifest if it exists, otherwise get them from command line
print(args$manifest)
if(!is.null(args$manifest)){
    pp <- projectParams(args$manifest)
} else {
    pp <- args
}

if(!is.null(pp)){
    logParams(pp)
}

if(is.null(pp$cell_type_markers)){
    stop("Need list of cell type markers in order to make pie charts.")
}

if(is.null(pp$counts_rda_file)){
    stop("Please provide counts with 'counts_rda_file' param")
}

if(is.null(pp$pdf_pie_charts_by_sample)){
    pdfFile <- gsub(".rda","_pie_charts.pdf",pp$counts_rda_file)
} else {
    pdfFile <- pp$pdf_pie_charts_by_sample
}

allTbls <- readRDS(pp$counts_rda_file)
countsTbl <- allTbls[["Counts"]]
markerNames <- trimws(unlist(strsplit(pp$cell_type_markers,",")))
plot_marker_percentages(countsTbl,markerNames,type="pie",other_threshold=as.numeric(pp$other_threshold),
                          exclude_sample_fov=pp$exclude_sample_fov, pdfFile=pp$pdf_pie_charts_by_sample,
                          v=pp$verbose,logFile=pp$log,debug=args$debug,custom_colors=pp$custom_colors
                       )
