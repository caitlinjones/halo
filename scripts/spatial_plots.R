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
parser$add_argument("--cell_types_file", type="character", default=NULL, help="file containing list of cell type markers to be included on plots; each line is one marker combination and may contain a single marker or a comma-separated list of markers (e.g., 'CD3' or 'CD3,CD8-,SOX10-'")
parser$add_argument("-d", "--data_file", type="character", help="*.rda files, containing halo data for one sample")
parser$add_argument("-a", "--annotations_dir", type="character", help="directory of Halo *.annotations files in XML format")
## optional args
parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="print extra output")
parser$add_argument("-l", "--log", type="character", default=gsub(" ","_",date()), help="log file")
parser$add_argument("--debug", action="store_true", default=FALSE, help="print extra output for debugging")

args <- parser$parse_args()
####################################################

pp <- NULL

## get all params from manifest if it exists, otherwise get them from command line
print(args$manifest)
if(!is.null(args$manifest)){
    pp <- initializeProject(args$manifest,type="spatial_plots")
    flog.info("Pulled params from manifest %s",args$manifest)
} else {
    pp <- args
}

if(!is.null(pp)){
    logParams(pp,"Making spatial plots")
}

## validate input
if(is.null(pp$cell_types_file)){
    stop("Please specify 'cell_types_file'.") 
}
if(is.null(pp$data_file)){
    stop("Please specify 'data_file'.")
}
if(is.null(pp$annotations_dir)){
    stop("Please specify 'annotations_dir'.")
}

flog.info("Plotting cell type locations")
pdfFile <- gsub("\\.rda","_cellTypeLocations.pdf",basename(pp$data_file))
plotCellTypeLocations(pp$data_file, pp$annotations_dir, pp$cell_types_file, pp$fov_bb, pp$plot_bb,
                      pp$boundary_colors, pp$cell_type_colors, pad=pp$pad, pdfFile=pdfFile)

flog.info("Plotting total density")
pdfFile <- gsub("\\.rda","_totalDensity.pdf",basename(pp$data_file))
plotTotalDensity(pp$data_file, pp$annotations_dir, pp$cell_types_file, pp$cell_type_name, 
                  pp$fov_bb, pp$pad, pp$plot_bb, pdfFile, logPlot=pp$log_plot,
                  funcMarker=pp$func_marker, sampleColor=pp$sample_color, sampleColorDark=pp$sample_color_dark,
                  exclude_sample_fov=pp$exclude_sample_fov, sortByMarker=pp$sort_by_marker,
                  exclude_sample_marker=pp$exclude_sample_marker,ymax=pp$ymax)

