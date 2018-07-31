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
parser$add_argument("-m", "--manifest", type="character", default=NULL, 
                    help="file containing all project parameters; run ?initializeProject for details")
parser$add_argument("--debug", action="store_true", default=FALSE, 
                    help="print extra output for debugging")

args <- parser$parse_args()
####################################################

usage <- function(){
    stop("Usage: Rscript pie_charts.R -m manifest.txt")
}

pp <- NULL

## get all params from manifest if it exists, otherwise get them from command line
print(args$manifest)
if(!is.null(args$manifest)){
    pp <- initializeProject(args$manifest,type="pie_charts")
} else {
    usage()
}


## validate input
if(is.null(pp$cell_type_markers)){
    stop("Need list of cell type markers in order to make pie charts.")
}
if(is.null(pp$counts_rda_file)){
    stop("Please provide counts with 'counts_rda_file' param")
}

## log parameters
if(!is.null(pp)){
    logParams(pp,"Making pie charts")
}

if(is.null(pp$pdf_pie_charts_by_sample)){
    pdfFile <- gsub(".rda","_pie_charts.pdf",pp$counts_rda_file)
} else {
    pdfFile <- pp$pdf_pie_charts_by_sample
}

allTbls <- readRDS(pp$counts_rda_file)
countsTbl <- allTbls[["Counts"]]
markerNames <- trimws(unlist(strsplit(pp$cell_type_markers,",")))
plotMarkerPercentages(countsTbl,markerNames,type="pie",other_threshold=as.numeric(pp$other_threshold),
                          exclude_sample_fov=pp$exclude_sample_fov, pdfFile=pp$pdf_pie_charts_by_sample,
                          custom_colors=pp$custom_colors
                       )
