library(contoureR)
library(digest)
library(futile.logger)
library(grid)
library(gridExtra)
library(kimisc)
library(magrittr)
library(plotrix)
library(randtoolbox)
library(raster)
library(RColorBrewer)
library(rgeos)
library(rJava)
library(scales)
library(SearchTrees)
library(sp)
library(tidyverse)
library(tools)
library(xlsx)
library(xlsxjars)
library(XML)
library(yaml)

source("/home/byrne/halo/dev/halo/R/melanoma_spatial.R")
source("/home/byrne/halo/dev/halo/R/util.R")
source("/home/byrne/halo/dev/halo/R/spatial_util.R")
source("/home/byrne/halo/dev/halo/R/marker_counts.R")
source("/home/byrne/halo/dev/halo/R/plots.R")

manifest <- "/home/byrne/halo/dev/halo/example/spatial_plots/yaml/PR_manifest.txt"
pp <- initializeProject(manifest)

dataFile = pp$data_file
annotationsDir = pp$annotations_dir
cellTypeFile = pp$cell_types_file
cellTypeName = pp$cell_type_name
fovBB = pp$fov_bb
pad=30
plotBB = pp$plot_bb
pdfFile = 'testing.pdf'
ymax = pp$ymax
log_plot=pp$log_plot
funcMarker = pp$func_marker
sortByMarker = pp$sort_by_marker
sampleColor = pp$sample_color
sampleColorDark = pp$sample_color_dark
exclude_sample_fov = pp$exclude_sample_fov
writeCSVfiles =pp$write_csv_files
haloInfiltrationDir = pp$halo_infiltration_dir
maxG = pp$max_g
outDir = pp$out_dir
byBand = pp$by_band
bandWidth = pp$band_width
maxDistanceFromInterface=pp$max_distance_from_interface

