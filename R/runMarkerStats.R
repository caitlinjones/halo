#!/opt/common/CentOS_6-dev/R/R-3.3.1/bin/R

###
# Parse user input
###
suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()

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

###############
### THIS WILL ALL BE CHANGED TO USE PACKAGE
###############
write("Loading libraries...",stdout())
lib.loc="/home/socci/lib/R"

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse,lib.loc=lib.loc))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(rJava,lib.loc=lib.loc))
options( java.parameters = c("-Xss2560k", "-Xmx8g") )
suppressPackageStartupMessages(library(xlsxjars,lib.loc=lib.loc))
suppressPackageStartupMessages(library(xlsx,lib.loc=lib.loc))
suppressPackageStartupMessages(library(digest))
suppressPackageStartupMessages(library(tools))
suppressPackageStartupMessages(library(kimisc,lib.loc=lib.loc))

source("halo.R")

pp <- projectParams("MelanomaV2_MANIFEST_20171227.txt")
if(!is.null(pp)){
    countMarkers(pp$markers, pp$dataDir, lf=pp$log, v=pp$verbose, 
}


countMarkers(args$markerFile, args$dataDir, lf=args$logFile, 
     v=args$verbose, pad=args$pad, outFile=NULL, runCounts=TRUE, 
     runFracTotal=TRUE, altBases=c("CD3","DAPI"))
