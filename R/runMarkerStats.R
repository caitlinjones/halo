#!/opt/common/CentOS_6-dev/R/R-3.3.1/bin/R

.libPaths(c("/home/byrne/R_libs","/home/socci/lib/R"))
options( java.parameters = c("-Xss2560k", "-Xmx8g") )

write("Loading libraries...",stdout())
library("halo")
suppressPackageStartupMessages(library("argparse"))

###
# Parse user input
###

parser <- ArgumentParser()

parser$add_argument("-m", "--manifest", type="character", default=NULL,
    help="file containing all project parameters; run ?projectParams for details"
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE, 
    help="print extra output [default]")
parser$add_argument("-m", "--markerFile", type="character", 
    help="comma separated file containing marker information")
parser$add_argument("-d", "--dataDir", type="character", default=".",
    help="directory containing *.rda files of image data")
parser$add_argument("-p", "--pad", type="integer", default=0,
    help="add this much padding around all sides of FOV")
parser$add_argument("-l", "--logFile", type="character", default=gsub(" ","_",date()),
    help="log file")

args <- parser$parse_args()


#### if manifest is given, default to using that regardless of any other arguments
#### if not, use arguments given on command line and set defaults for others
if(!is.null(args$manifest)){
    pp <- projectParams(args$manifest)
    if(!is.null(pp)){
        countMarkers(pp$markers, pp$data_dir, lf=pp$log, v=pp$verbose, pad=pp$pad, outFile=pp$out_file,
                 runCounts=pp$run_counts, runFracTotal=pp$run_frac_total, runMedians=pp$run_medians,
                 altBases=pp$alt_bases)
    }
} else {
    countMarkers(args$markerFile, args$dataDir, lf=args$logFile, 
         v=args$verbose, pad=args$pad, outFile=NULL, runCounts=TRUE, 
         runFracTotal=TRUE, runMedians=TRUE, altBases=c("CD3","DAPI"))
}
