options( java.parameters = c("-Xss2560k", "-Xmx8g") )
options('python_cmd'='/opt/common/CentOS_6-dev/python/python-3.5.1/bin/python')

write("Loading libraries...",stdout())
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("halodev"))

source("/home/byrne/halo/dev/halodev/R/constants.R")
source("/home/byrne/halo/dev/halodev/R/process_meta.R")
source("/home/byrne/halo/dev/halodev/R/util.R")
source("/home/byrne/halo/dev/halodev/R/spatial_util.R")
source("/home/byrne/halo/dev/halodev/R/marker_combo_table.R")
source("/home/byrne/halo/dev/halodev/R/marker_counts.R")
source("/home/byrne/halo/dev/halodev/R/plots.R")
source("/home/byrne/halo/dev/halodev/R/interface.R")

###
## Parse user input
###
parser <- ArgumentParser()

## args are preferrably all in manifest
parser$add_argument("-i", "--input", type="character", default=NULL,
                    help="XLSX file containing table of marker combinations and counts")
parser$add_argument("-o", "--outDir", type="character", default=getwd(),
                    help="output directory")
parser$add_argument("-m", "--metaDir", type="character", default=NULL,
                    help="directory of study meta data files (XLSX files created by lab member")
args <- parser$parse_args()

print(args)

###
## check for meta data files
###
if(!is.null(args$metaDir)){
    metaFiles <- file.path(args$metaDir,dir(args$metaDir)[grep("\\.xlsx",dir(args$metaDir))])
    print("Meta files:")
    print(paste0("  ",metaFiles))
} else {
    stop("No meta data given.")
}
###
## remove files that someone currently has open, indicated by "~" in basename
###
if(length(grep("~",metaFiles)) > 0){
    metaFiles <- metaFiles[-grep("~",metaFiles)]
}

flatMeta <- flattenMetaData(metaFiles=metaFiles)
allCellTypes <- flatMeta$MarkerCombinationAnnotation

tbl <- cellDive.fillInMarkerComboInterpretations(xlfile=args$input, outDir=args$outDir, allCellTypes=allCellTypes)
