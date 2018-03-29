#!/home/byrne/R/R-3.4.3/bin/R

options( java.parameters = c("-Xss2560k", "-Xmx8g") )
options('python_cmd'='/opt/common/CentOS_6-dev/python/python-3.5.1/bin/python')

write("Loading libraries...",stdout())
suppressPackageStartupMessages(library("argparse"))
#suppressPackageStartupMessages(library("halo"))
suppressPackageStartupMessages(library("halodev"))

###
# Parse user input
###
parser <- ArgumentParser()
## args are preferrably all in manifest
parser$add_argument("-m", "--manifest", type="character", default=NULL, 
                    help="file containing all project parameters; run ?initializeProject for details")
parser$add_argument("--cellTypeLocations", action="store_true", default=FALSE, 
                    help="plot cell type locations for each sample FOV")
parser$add_argument("--plotDensity", action="store_true", default=FALSE, 
                    help="plot total density of cells within a set distance from a tumor boundary")
parser$add_argument("--debug", action="store_true", default=FALSE, 
                    help="print extra output for debugging")

args <- parser$parse_args()
####################################################

usage <- function(){
    stop("Usage: Rscript spatial_plots.R -m manifest.txt")
}

pp <- NULL

## eventually set up defaults for everything, then allow everything to be
## passed on command line; for now require manifest
print(args$manifest)
if(!is.null(args$manifest)){
    pp <- initializeProject(args$manifest,type="spatial_plots")
    flog.info("Pulled params from manifest %s",args$manifest)
} else {
    usage()
}

if(!is.null(pp)){
    logParams(pp,"Making spatial plots")
}


#flog.info("Plotting cell type locations")
#pdfFile <- gsub("\\.rda","_cellTypeLocations.pdf",basename(pp$data_file))
#plotCellTypeLocations(pp$data_file, pp$annotations_dir, pp$cell_types_file, pp$fov_bb, pp$plot_bb,
#                      pp$boundary_colors, pp$cell_type_colors, pad=pp$pad, pdfFile=pdfFile)

flog.info("Plotting total density")

plotDensity(pp$data_files, 
            pp$annotations_dirs, 
            pp$marker_set_files, 
            pp$cell_type_name, 
            pp$fov_bb, 
            pp$pad, 
            pp$pdf_file, 
            densityFiles=pp$density_files,
            logPlot=pp$log_plot,
            funcMarker=pp$func_marker, 
            funcPosColor=pp$functional_pos_color, 
            funcNegColor=pp$functional_neg_color,
            sampleColor=pp$sample_color, 
            sampleColorDark=pp$sample_color_dark, 
            exclude_sample_fov=pp$exclude_sample_fov, 
            sortByMarker=pp$sort_by_marker, 
            sampleOrder=pp$sample_order,
            outDir=pp$out_dir,
            writeCSVfiles=pp$write_csv_files, 
            maxG=pp$max_g,
            byBand=pp$by_band,
            bandWidth=pp$band_width, 
            maxDistanceFromInterface=pp$max_distance_from_interface,
            plotPer=pp$plot_per,
            filePer=pp$file_per,
            includeSortByMarker=pp$include_sort_by_marker,
            consistentYScaling=pp$consistent_y_scaling,
            boxPlots=pp$box_plots)

