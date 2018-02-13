# plotFOV, plotFOV0, plotFOV1 ... all versions of the same thing?

#' Plot location of given cell types within one FOV
#'
#' Given a halo data *.rda file with XML annotations and a list of cell types, 
#' plot the X and Y coordinates of those cells showing locations relative to
#' each other and to tissue boundary
#'
#' @param dataFile         *.rda file containing ObjectAnalysisData
#' @param annotationsDir   directory of *.annotations XML files from Halo for one sample
#' @param cellTypesFile    text file containing a list of cell type markers, one on each line; 
#'                         each cell type marker may consist of an indivdual marker or a comma-
#                          separated list of markers (e.g., 'CD3,CD8,SOX10-')
#' @param pad              amount to trim from FOV; default=30
#' @param pdfFile          name of output PDF file
#' @param verbose          print progress messages; default=FALSE
#' @param logFile          log file; default=NULL
#' @param debug            print debug messages; default=FALSE
#' @return  nothing  
#' @export
plotCellTypeLocations <- function(dataFile, annotationsDir, cellTypesFile, pad=30, 
                                 pdfFile=NULL, verbose=FALSE, logFile=NULL, debug=FALSE){

    ## CONSTANTS
    pixel2um <- 0.293
    p2tomm <- pixel2um^2/(1000^2)
    p2toum <- pixel2um^2

    bbFOV0 <-list(X0=1,Y0=-3343,X1=5363,Y1=1) ## solid gray rectangle
    #bbFOV0 <- list(X0=1,Y0=-3377,X1=5363,Y1=1)
    bbPlot <- bbFOV0 ## outermost black rectangle
    bbPlot$X0 <- bbPlot$X0-750 ## add 500 padding to left
    bbPlot$Y1 <- bbPlot$Y1+500 ## add 500 padding on top

    bbData <- bbFOV0
    fovTag <- "ORIG"
    clipSize <- 0

    bndColors=list(Tum="darkred",Exc="grey",Epi="darkgreen")
    cls <- c("#2166ac","#b2182b","#c2a5cf")
    cls <- c("#377eb8","#e41a1c","#cab2d6") ## colors as input?

    aFiles <- dir(annotationsDir)[grep("\\.annotations",dir(annotationsDir))]
    epiFiles <- grep("epi",aFiles)
    if(length(epiFiles) > 0){
        aFiles <- aFiles[-epiFiles]
    }
    aFiles <- file.path(annotationsDir,aFiles)
    aFileSpots <- as.numeric(gsub(".*_Spot|\\.annotations","",aFiles))

    ## process rda file
    logMsg(paste0("Reading data file ",dataFile),verbose,logFile)
    dd <- readRDS(dataFile)
    sampleName <- gsub("_.*","",basename(dataFile))

    if(debug){ logMsg("Getting FOVs",verbose,logFile) }
    spots <- dd %>% select(SPOT) %>% distinct(SPOT) %>% arrange(SPOT) %>% pull(SPOT) 

    if(pad > 0){
        fovTag <- paste0("clip",pad)
        pdfFile <- gsub("\\.pdf",paste0("_",fovTag,".pdf"),pdfFile)
        clipSize <- pad
        bbData <- padBoundingBox(bbFOV0,-clipSize/pixel2um)
    }

    cellTypes <- scan(cellTypesFile,"")
    names(cls) <- cellTypes

    pdf(pdfFile,width=11,height=8.5)
    for(i in seq(spots)){
        spot = spots[i]
        logMsg(paste0("Getting data for spot ",spot,"..."),verbose,logFile) 
        ds=dd %>%
            filter(SPOT==spot & ValueType=="Positive") %>%
            mutate(X=(XMax+XMin)/2,Y=-(YMax+YMin)/2) %>%
            select(Sample,SPOT,UUID,X,Y,Marker,Value) %>%
            spread(Marker,Value)
        logMsg(paste0("  plotting..."),verbose,logFile) 
        ## bbPlot = black solid (entire plot)
        ## bbFOV0 = gray solid  (fov max/min from Halo)
        ## bbData = gray dotted (fov minus padding)
        plotFOV1(bbData,sampleName,spot,bbPlot,bbFOV0)
        ## draw boundaries
        sTag <- paste(sampleName, paste0("Spot",spot,".annotations"),sep="_")
        aai <- grep(sTag,aFiles)
        if(len(aai)>0){
            if(debug){ logMsg(paste0("  Found boundary annotation file ",aFiles[aai],". Adding boundaries.."),verbose,logFile) }
            aFile=aFiles[aai]
            boundaries=readHaloAnnotations(aFile,verbose,logFile)
            boundaries %>%
                walk(function(x){
                    color=ifelse(x$RegionCode=="Tum",1,8);
                    with(x,lines(X,Y,col=color,lwd=2))
                    }
                )
        }

        ## draw locations of all cells
        with(ds,points(X,Y,col="#EEEEEE",pch=4,cex=.35))

        ## draw locations of cells matching cell types in file
        for(ci in seq(cellTypes)) {
            ds %>%
                filter_(markerStringToPredicate(cellTypes[ci])) %$%
                points(X,Y,pch=16,col=cls[cellTypes[ci]],cex=.8)
        }
        legend(bbPlot$X0,bbPlot$Y1,c("DAPI",cellTypes),col=c(8,cls[cellTypes]),pch=c(4,16,16,16),bg="#FFFFFF")

    }
    dev.off()
}
