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
parser$add_argument("-m", "--manifest", type="character", default=NULL, help="file containing all project parameters; run ?initializeProject for details")
## args required if manifest not given
parser$add_argument("--cell_type_markers", type="character", default=NULL, help="comma-separated string of markers, for each of which a pie chart will be generated; e.g., 'CD3,CD4,CD8,CD20'")
parser$add_argument("-c", "--counts_rda_file", type="character", help="*.rda file containing list of counts-related tables (output of runMarkerCounts.R)")
parser$add_argument("--pie_charts", action="store_true", default=FALSE, help="generate pie charts for cell type markers")
parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="print extra output")
parser$add_argument("-l", "--log", type="character", default=gsub(" ","_",date()), help="log file")
parser$add_argument("--debug", action="store_true", default=FALSE, help="print extra output for debugging")
## args required if manifest not given and counts not yet generated
parser$add_argument("-mf", "--markers", type="character", help="comma separated file containing marker information")
parser$add_argument("-d", "--data_dir", type="character", default=".", help="directory containing *.rda files of image data")
parser$add_argument("-p", "--pad", type="integer", default=0, help="add this much padding around all sides of FOV")
parser$add_argument("--run_frac_total", action="store_true", help="add sheet to counts file containing fractions of total cells")
parser$add_argument("--run_medians", action="store_true", help="add sheet to counts file containing median counts")
parser$add_argument("--alt_bases", default=NULL, help="comma separated string of markers to use as alternate 'baselines'")
parser$add_argument("--counts_xlsx_file", type="character", help="*.xlsx file of count")

args <- parser$parse_args()
####################################################

pp <- NULL

## get all params from manifest if it exists, otherwise get them from command line
print(args$manifest)
if(!is.null(args$manifest)){
    pp <- initializeProject(args$manifest,type="pipeline")
} else {
    pp <- args
}

if(!is.null(pp)){
    logParams(pp)
}

countsTbl <- NULL

### run counts first, if needed
if(pp$run_counts){
    if(is.null(pp$markers) || is.null(pp$data_dir)){
        stop("The following args need to be set in order to run counts: 'markers' and 'data_dir'")
    }
    allTbls <- countMarkers(pp$markers, pp$data_dir, lf=pp$log, v=pp$verbose, pad=pp$pad,
                     countsXLSXFile=pp$counts_xlsx_file,
               countsRDAFile=pp$counts_rda_file, runCounts=pp$run_counts,
                     runFracTotal=pp$run_frac_total, runMedians=pp$run_medians,
                     altBases=pp$alt_bases,debug=args$debug)
    countsTbl <- allTbls[["Counts"]]
}

if(pp$pie_charts){
    pdfFile <- "pie_charts.pdf"
    if(is.null(pp$cell_type_markers)){
        stop("Need list of cell type markers in order to make pie charts.")
    }
    if(is.null(countsTbl)){
        if(is.null(pp$counts_rda_file)){
            stop("Something is not right. Counts table does not exist and *.rda file was not given")
        }
        allTbls <- readRDS(pp$counts_rda_file)
        countsTbl <- allTbls[["Counts"]]
        pdfFile <- gsub(".rda","_pie_charts.pdf",pp$counts_rda_file)
    } else {
        pdfFile <- gsub(".txt","_pie_charts.pdf",pp$markers)
    }
    markerNames <- trimws(unlist(strsplit(pp$cell_type_markers,",")))
    plotMarkerPercentages(countsTbl,markerNames,type="pie",other_threshold=as.numeric(pp$other_threshold),
                              exclude_sample_fov=pp$exclude_sample_fov, pdfFile=pp$pdf_pie_charts_by_sample,
                              v=pp$verbose,logFile=pp$log,debug=args$debug,custom_colors=pp$custom_colors)
}
