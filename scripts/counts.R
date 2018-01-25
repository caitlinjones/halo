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
parser$add_argument("-mf", "--markers", type="character", help="comma separated file containing marker information")
parser$add_argument("-d", "--data_dir", type="character", default=".", help="directory containing *.rda files of image data")

## optional arguments
parser$add_argument("-p", "--pad", type="integer", default=0, help="add this much padding around all sides of FOV")
parser$add_argument("--run_frac_total", action="store_true", help="add sheet to counts file containing fractions of total cells")
parser$add_argument("--run_medians", action="store_true", help="add sheet to counts file containing median counts")
parser$add_argument("--alt_bases", default=NULL, help="comma separated string of markers to use as alternate 'baselines'")
parser$add_argument("--counts_xlsx_file", type="character", help="*.xlsx file of count")

parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="print extra output")
parser$add_argument("-l", "--log", type="character", default=gsub(" ","_",date()), help="log file")
parser$add_argument("--debug", action="store_true", default=FALSE, help="print extra output for debugging")

args <- parser$parse_args()
####################################################

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


### run counts first, if needed
if(is.null(pp$markers) || is.null(pp$data_dir)){
    stop("The following args need to be set in order to run counts: 'markers' and 'data_dir'")
}
allTbls <- countMarkers(pp$markers, pp$data_dir, lf=pp$log, v=pp$verbose, pad=pp$pad,
                     countsXLSXFile=pp$counts_xlsx_file,
                     countsRDAFile=pp$counts_rda_file, runCounts=pp$run_counts,
                     runFracTotal=pp$run_frac_total, runMedians=pp$run_medians,
                     altBases=pp$alt_bases, debug=args$debug,
                     saveRDSfile=pp$save_rds_file, writeXLSXfile=pp$write_xlsx_file)

