options( java.parameters = c("-Xss2560k", "-Xmx8g") )

write("Loading libraries...",stdout())
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("halodev"))

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
    #pp <- initializeProject(args$manifest,type="counts")
    pp <- read_yaml(args$manifest)
    print(pp)
} else {
    usage()
}

## check for data
if(!is.null(pp$raw_data_files)){
    dataFiles <- pp$raw_data_files
} else if(!is.null(pp$raw_data_dir)){
    dataFiles <- file.path(pp$raw_data_dir,dir(pp$raw_data_dir)[grep("\\.rda",dir(pp$raw_data_dir))])
} else {
    stop("No data given.")
}
print("Data files:")
print(paste0("  ",dataFiles))

## check for drift/loss summaries
if(!is.null(pp$drift_files)){
    driftFiles <- pp$drift_files
    print("Drift files:")
    print(paste0("  ",driftFiles))
} else if(!is.null(pp$drift_dir)){
    driftFiles <- file.path(pp$drift_dir,dir(pp$drift_dir)[grep(".*summary.*txt",dir(pp$drift_dir))])
    print("Drift files:")
    print(paste0("  ",driftFiles))
} else {
    driftFiles <- NULL
    warning("No drift summaries given.")
}

## check for meta data files
if(!is.null(pp$meta_files)){
    metaFiles <- pp$meta_files
    print("Meta files:")
    print(paste0("  ",metaFiles))
} else if(!is.null(pp$meta_dir)){
    metaFiles <- file.path(pp$meta_dir,dir(pp$meta_dir)[grep("\\.xlsx",dir(pp$meta_dir))])
    print("Meta files:")
    print(paste0("  ",metaFiles))
} else {
    stop("No meta data given.")
}

## read study annotations
sampleAnnFile <- metaFiles[grep("SampleAnnotations.xlsx",metaFiles)]
if(length(sampleAnnFile) == 0){
    stop("No sample annotation found.")
}
sampAnn <- as.tibble(read.xlsx(sampleAnnFile,1))

fovAnnFile <- metaFiles[grep("FOVannotations.xlsx",metaFiles)]
if(length(fovAnnFile) == 0){
    stop("No FOV annotation found.")
}
fovAnnotations <- as.tibble(read.xlsx(fovAnnFile,1))

##### some minimal validation
borderPad <- pp$pad
if(is.null(borderPad)){
    error("No pad given. Please set in manifest file") 
}
driftThreshold <- pp$drift_threshold
if(is.null(driftThreshold)){
    error("No drift threshold given. Please set in manifest file") 
}

## figure out where to write rda files with exclusion column
outDir <- getwd()
if("data_dir" %in% names(pp) && !is.null(pp$data_dir)){
    outDir <- pp$data_dir
} else if("out_dir" %in% names(pp) && !is.null(pp$out_dir)){
    outDir <- pp$out_dir
}
if(!file.exists(outDir)){
    dir.create(outDir, recursive=T, showWarnings=FALSE)
}

for(x in seq(dataFiles)){
    df <- dataFiles[x]

    print(paste0("Reading file ",df))
    dat <- readRDS(df)
    #dat$Sample <- gsub("_.*","",dat$Sample) #### FIGURE OUT A BETTER WAY

    ## to map FOV annotations to samples
    samp <- unique(dat$Sample)
    cellDiveId <- sampAnn %>% filter(Sample_name == samp) %>% pull(CELL_DIVE_ID) %>% as.character()
    #cellDiveId <- samp   #### RIGHT NOW THIS IS INCONSISTENT BETWEEN COHORTS

    drift <- NULL
    ## figure out which drift file to use
    if(!is.null(driftFiles)){
        dlFile <- driftFiles[grep(paste0(cellDiveId,"_"), driftFiles)]
        if(length(dlFile) == 0){
            warning("No drift file found. Skipping drift exclusions")
        } else {
            print(paste0("Reading drift file",dlFile))
            drift <- as.tibble(read.delim(dlFile, header=T, sep="\t"))
        }
    }
    print(paste0("Marking exclusions for ",samp))
    dat <- markExclusions(samp, dat, drift, fovAnnnotations, cellDiveId, printPlots=TRUE)

    saveRDS(dat, file=file.path(outDir,gsub("\\.rda","_Excl.rda",basename(df))))
}

