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
parser$add_argument("--plotDensityByBand", action="store_true", default=FALSE,
                    help="plot infiltration density by set distance bands from tumor boundaries")
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

## read data files 
dat <- NULL
if(!is.null(pp$data_files)){
    flog.info("Reading data files")
    for(df in pp$data_files){
        dd <- readRDS(df)
        dat <- dat %>% bind_rows(dd)
    }
    #dat$Sample <- gsub("_ObjectAnalysisData","",dat$Sample)
    ### TEMPORARY UNTIL EXCLUSIONS ARE DONE UPSTREAM
    #dat <- filter(dat, !(Sample == "Untreated" & Marker=="CD68" & SPOT %in% c(5, 11,12,13,14, 18, 19, 20, 21, 22, 23, 24)))
    dat <- dat %>% filter(EXCLUDE == "")
}

## read area files if they exist
ia <- ba <- NULL
if(!is.null(pp$area_files)){
    flog.info("Reading area files")
    ia <- ba <- NULL
    for(af in pp$area_files){
        a <- readRDS(af)
        ia <- ia %>% bind_rows(a$area)
        ba <- ba %>% bind_rows(a$bandAssignments)
    }
    ia$Sample <- gsub("_ObjectAnalysisData","",ia$Sample)
    ba$Sample <- gsub("_ObjectAnalysisData","",ba$Sample)
    ba <- ba[!is.na(ba$Band),]
    if(!pp$by_band){
        ia <- ia %>% mutate(`(-360,360]`=rowSums(ia[,3:ncol(ia)])) %>%
               select(Sample,SPOT,`(-360,360]`)
        ba$Band <- "(-360,360]"
    }
}

## plot cell type locations
if(args$cellTypeLocations){
    flog.info("Plotting cell type locations")
    pdfFile <- gsub("\\.rda","_cellTypeLocations.pdf",basename(pp$data_file))
    plotCellTypeLocations(pp$data_file, pp$annotations_dir, pp$cell_types_file, pp$fov_bb, pp$plot_bb,
                          pp$boundary_colors, pp$cell_type_colors, pad=pp$pad, pdfFile=pdfFile)
}

## plot infiltration density
if(args$plotDensityByBand){
    flog.info("Plotting density by distance interval")
    if(is.null(pp$marker_config_file)){
        stop("Can not plot infiltration density. Plot config required.")
    }
    ## get marker config
    markerCfg <- halodev::getMarkerConfig(pp$marker_config_file,pp$plot_config_file)
    markers <- unique(markerCfg$CellType)

    ## get density
    flog.info(paste0("Getting density for ",length(markers)," marker combinations...\n"))
    den <- halodev::getMarkerDensityTable(markers, dat=dat, dataFiles=pp$data_files, annotationsDirs=pp$annotations_dirs, 
                                 outDir=pp$out_dir, areaFiles=pp$area_files, densityFiles=pp$density_files, 
                                 sampleOrder=pp$sample_order, writeCSVfiles=pp$write_csv_files, pad=pp$pad, 
                                 byBand=pp$by_band, maxDistanceFromInterface=pp$max_distance_from_interface,
                                 bandWidth=pp$band_width, markerSetName=pp$marker_set_name, 
                                 sortByMarker=pp$sort_by_marker, funcMarker=pp$func_arker,
                                 maxG=pp$max_g, calcFOVarea=pp$calc_fov_area)

    #### TEMPORARY UNTIL FILTERING IS DONE UPSTREAM
    #den <- den %>% filter(!( Sample == "Untreated" & SPOT %in% c(5, 11,12,13,14, 18, 19, 20, 21, 22, 23, 24) & (grepl("^CD68$",CellType) | grepl("CD68,",CellType))))

    ## plot plots
    flog.info("Printing plots")
    halodev::printDensityPlots(den, pp$band_width, markerCfg, yScaleConsistency=pp$infiltration_y_scale_consistency, 
                      absoluteDensity=pp$infiltration_absolute_density,
                      densityPercentage=pp$infiltration_density_percentages, 
                      byFOV=pp$infiltration_by_fov, summarize=pp$infiltration_summary,
                      sampleOrder=pp$sample_order, stacked=pp$infiltration_stacked,
                      separateLegend=pp$infiltration_legend_separate)
    
}

## plot total infiltration density
if(args$plotDensity){
    flog.info("Plotting total density")

    plotDensity(pp$data_files, 
                pp$annotations_dirs, 
                pp$marker_set_files, 
                pp$cell_type_name, 
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
}
