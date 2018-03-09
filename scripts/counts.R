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
    stop("Usage: Rscript counts.R -m manifest.txt")
}

pp <- NULL

## get all params from manifest if it exists, otherwise get them from command line
print(args$manifest)
if(!is.null(args$manifest)){
    pp <- initializeProject(args$manifest,type="counts")
} else {
    usage()
}

## validate input
if(is.null(pp$markers) || is.null(pp$data_dir)){
    stop("The following args need to be set in order to run counts: 'markers' and 'data_dir'")
}

if(!is.null(pp)){
    logParams(pp,"Generating counts")
}

allTbls <- countMarkers(pp$markers, pp$data_dir, pad=pp$pad, altBases=pp$alt_bases,
                     countsXLSXFile=pp$counts_xlsx_file, countsRDAFile=pp$counts_rda_file, 
                     runCounts=pp$run_counts, runFracTotal=pp$run_frac_total,
                     runMedians=pp$run_medians, saveRDSfile=pp$save_rds_file, 
                     writeXLSXfile=pp$write_xlsx_file)

