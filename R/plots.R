#' Extract a legend from a gplot
#' 
#' Extract a legend from a plot generated with ggplot
#' 
#' @param plt    a plot generated with ggplot 
getLegend <- function(plt){
    tmp <- tryCatch({ 
               ggplot_gtable(ggplot_build(plt))
             }, error=function(e){
               plt
             }) 
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
    legend <- tmp$grobs[[leg]] 
    return(legend)
}

#' Remove a legend from a gplot
#' 
#' Remove a legend from a gplot
#'
#' @param plt    a plot generated with ggplot 
removeLegend <- function(plt){
    tmp <- ggplot_gtable(ggplot_build(plt))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    tmp$grobs[[leg]] <- NULL
    return(tmp)
}

gridArrangeSharedLegend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {

  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)

  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))

  grid.newpage()
  grid.draw(combined)

  # return gtable invisibly
  invisible(combined)

}


######### "SPATIAL" PLOTS ########

### FIGURE OUT WHAT EACH OF THESE DOES AND RENAME
plotFOV<-function(bbData,bnd,sampleName,spot,interfaceSize=boundaryInterfaceSize) {
    bb=getBoundingBoxL(bnd %>% bind_rows)
    bb=joinBoundingBoxes(bb,bbFOV0)
    bb=padBoundingBox(bb,1.1*interfaceSize)

    spMax=SpatialPolygons(list(
        Polygons(list(Polygon(boundingBoxToRect(bb))),"maxBB")))

    spFOV=SpatialPolygons(list(
        Polygons(list(Polygon(boundingBoxToRect(bbFOV0))),"FOV")))

    ## type="n": blank canvas
    plot(1,type="n",xlim=c(bb$X0,bb$X1),ylim=c(bb$Y0,bb$Y1),
        xlab="",ylab="",main=paste(sampleName,spot))

    ## col 8 = gray
    lines(spFOV,col=8,lwd=2)

    spData=SpatialPolygons(list(
        Polygons(list(Polygon(boundingBoxToRect(bbData))),"Data")))

    lines(spData,col=1,lty=2)

    lines(spMax,col=8,lty=3)

}

plotFOV0<-function(bbData,sampleName,spot,bbPlot,bbFOV0) {
    spFOV=SpatialPolygons(list(
        Polygons(list(Polygon(boundingBoxToRect(bbFOV0))),"FOV")))
    bb=bbFOV0
    plot(1,type="n",xlim=c(bb$X0,bb$X1),ylim=c(bb$Y0,bb$Y1),
        xlab="",ylab="",main=paste(sampleName,spot))

    ## solid gray rectangle
    lines(spFOV,col=8,lwd=2)

    spData=SpatialPolygons(list(
        Polygons(list(Polygon(boundingBoxToRect(bbData))),"Data")))

    ## black dotted line
    lines(spData,col=1,lty=2)

}

plotFOV1<-function(bbData,sampleName,spot,bbPlot,bbFOV0) {
    ## solid gray rectangle
    spFOV=SpatialPolygons(list(
        Polygons(list(Polygon(boundingBoxToRect(bbFOV0))),"FOV")))
    bb=bbFOV0
    plot(1,type="n",xlim=c(bbPlot$X0,bbPlot$X1),ylim=c(bbPlot$Y0,bbPlot$Y1),
        xlab="",ylab="",main=paste(sampleName,spot))

    lines(spFOV,col=8,lwd=2)

    spData=SpatialPolygons(list(
        Polygons(list(Polygon(boundingBoxToRect(bbData))),"Data")))

    lines(spData,col=1,lty=2)

}


#' Plot location of given cell types within one FOV
#'
#' Given a halo data *.rda file with XML annotations and a list of cell types, 
#' plot the X and Y coordinates of those cells showing locations relative to
#' each other and to tissue boundary
#'
#' @param dataFiles        *.rda file containing ObjectAnalysisData
#' @param annotationsDirs  directory of *.annotations XML files from Halo for one sample
#' @param cellTypesFile    text file containing a list of cell type markers, one on each line; 
#'                         each cell type marker may consist of an indivdual marker or a comma-
#                          separated list of markers (e.g., 'CD3,CD8,SOX10-')
#' @param pad              amount to trim from each FOV
#' @param boundaryColors   a list of hexadecimal color codes for Tum, Exc, and Epi boundaries 
#' @param cellTypeColors   a list of hexadecimal color codes for each cell type in cellTypesFile
#' @param pdfFile          name of output PDF file
#' @return  nothing  
#' @export
plotCellTypeLocations <- function(dataFiles, annotationsDirs, cellTypesFile,
                                  boundaryColors, cellTypeColors, pad=30, 
                                  pdfFile=NULL){

    aFiles <- dir(annotationsDir)[grep("\\.annotations",dir(annotationsDir))]
    epiFiles <- grep("epi",aFiles)
    if(length(epiFiles) > 0){
        aFiles <- aFiles[-epiFiles]
    }
    aFiles <- file.path(annotationsDir,aFiles)
    aFileSpots <- as.numeric(gsub(".*_Spot|\\.annotations","",aFiles))

    ## process rda file
    flog.info("Reading data file %s",dataFile)
    dd <- readRDS(dataFile)
    sampleName <- getSampleFromFileName(dataFile)

    flog.debug("Getting FOVs")
    spots <- dd %>% dplyr::select(SPOT) %>% distinct(SPOT) %>% arrange(SPOT) %>% pull(SPOT) 


    cellTypes <- scan(cellTypesFile,"")

    pdf(pdfFile,width=11,height=8.5)
    for(i in seq(spots)){
        spot = spots[i]
        flog.debug("Getting data for spot %s",spot)
        ds=dd %>%
            filter(SPOT==spot & ValueType=="Positive") %>%
            mutate(X=(XMax+XMin)/2,Y=-(YMax+YMin)/2) %>%
            dplyr::select(Sample,SPOT,UUID,X,Y,Marker,Value) %>%
            spread(Marker,Value)

         bbFOV <- getBoundingBoxL(ds)
         bbPlot <- list(X0=bbFOV$X0-500, X1=bbFOV$X1+500, Y0=bbFOV$Y0-500, Y1=bbFOV$Y1+500) 
         bbData <- bbFOV
         if(pad > 0){
           bbData <- padBoundingBox(bbFOV,-pad/pixel2um)

         }       
 
        flog.debug("  plotting...")
        ## bbPlot = black solid (entire plot)
        ## bbFOV0 = gray solid  (fov max/min from Halo)
        ## bbData = gray dotted (fov minus padding)
        plotFOV1(bbData,sampleName,spot,bbPlot,bbFOV)
        ## draw boundaries
        sTag <- paste(sampleName, paste0("Spot",spot,".annotations"),sep="_")
        aai <- grep(sTag,aFiles)
        if(len(aai)>0){
            flog.debug("Found boundary annotation file %s. Adding boundaries..",aFiles[aai])
            aFile=aFiles[aai]
            boundaries=readHaloAnnotations(aFile)
            boundaries %>%
                walk(function(x){
                    color=ifelse(x$RegionCode=="Tum",boundaryColors[["Tum"]],boundaryColors[["Exc"]]);
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
                points(X,Y,pch=16,col=cellTypeColors[[cellTypes[ci]]],cex=.8)
        }
        legend(bbPlot$X0,bbPlot$Y1,c("DAPI",cellTypes),col=c(8,unlist(cellTypeColors)[cellTypes]),pch=c(4,rep(16,len(cellTypes))),bg="#FFFFFF")

    }
    dev.off()
}

#' Get legend from a gplot
#' 
#' Given a plot generated with ggplot(), separate and return only
#' the legend component
#' 
#' @param a.gplot  a plot generated by ggplot()
#' @return  the legend of a.gplot 
g_legend <- getLegend 

#' Get plot theme for a given type of halo plot
#' 
#' Depending on what type of plot you're making, get a theme()
#' object with the elements standard to that plot type
#'
#' @param plotType  character string indicating which type of plot youre making;
#'                  options: ["infiltration density"]
#' @return theme() object with standard theme elements
#' @export
getPlotTheme <- function(plotType){

    if(plotType == "infiltration density"){
        return(theme(legend.title = element_blank(),
                  legend.position = "right",
                  legend.text = element_text(size=12),
                  axis.text = element_text(size=12),
                  axis.text.x = element_text(size=8),
                  axis.title.y = element_text(size=14),
                  axis.title.x = element_text(size=14),
                  axis.ticks = element_blank(),
                  axis.line = element_line(color = "black", size=0.75),
                  strip.text.x = element_text(size=12, margin=margin(.1, .1, .1, .1, "cm")),
                  strip.placement = "outside",
                  plot.background = element_blank(),
                  plot.title = element_text(size=20, margin=margin(b=0.5, t=0.2, r=0, l=0.1, "cm")),
                  panel.background = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()
                  )
               )
    }

    if(plotType == "infiltration density single page"){
        return(theme(legend.position="none",
                     axis.text.x = element_text(size = 8),
                     axis.text.y = element_text(size = 8),
                     axis.title.x = element_text(size = 8),
                     axis.title.y = element_text(size = 8),
                     plot.title = element_text(size = 12)
                    )
              )
    }

    if(plotType == "total FOV density"){
        return(theme(legend.title = element_text(size=12, face="bold"),
                  legend.position = "right",
                  legend.text = element_text(size=12),
                  axis.text = element_text(size=12),
                  axis.text.x = element_text(size=8,hjust=1,angle=45),
                  axis.title.y = element_text(size=14),
                  axis.title.x = element_text(size=14),
                  axis.ticks = element_blank(),
                  axis.line = element_line(color = "black", size=0.75),
                  strip.text.x = element_text(size=12, margin=margin(.1, .1, .1, .1, "cm")),
                  strip.placement = "outside",
                  plot.background = element_blank(),
                  plot.title = element_text(size=20, margin=margin(b=0.5, t=0.2, r=0, l=0.1, "cm")),
                  panel.background = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()
                  )
               )
    }

}

#' Plot and print absolute density or density percentages of cell types near tumor
#' interface
#' 
#' Plot and print absolute density or density percentages of cell types near tumor
#' interface
#'
#' @param mDen             tibble containing 
#' @param densityMarkers   a vector of markers for which density is to be plotted
#' @param bandWidth        if plotting by distance intervals from infiltration boundary, this is the size of 
#'                         each interval
#' @param clrs             a list of colors, named by either by CellType from mDen tibble, 
#'                         or cellTypeLabels if specified
#' @param plotTitle        title of plot
#' @param cellTypeLabels   a list of labels for each marker being plotted, if different
#'                         than CellType in mDen tibble; default: NULL
#' @param sampleOrder      order in which samples should appear on the plot; default=NULL
#' @param yMax             if specified, y axis will span from zero to this number; default=NULL
#' @param separateLegend   logical; when TRUE, legend will be separated from plot. if printLegend is TRUE, 
#'                         will be printed first, before plot. if printLegend is FALSE, will not be printed at all; default=FALSE
#' @param printLegend      logical; when separateLegend is TRUE, if printLegend is FALSE, no legend will be printed at all; default=TRUE
#' @param legendOnly       logical; when separateLegend is TRUE, if legendOnly is TRUE, only legend will be printed, no plot
#' @param yCol             character string specifying which of the columns in mDen is to be used for y values; default="Density"
#' @param pct              logical; values being plotted are percentages; default=FALSE
#' @param facetByFOV       logical; plot each FOV separately; default=TRUE
#' @param facetByCellType  logical; plot each cell type separately; default=FALSE
#' @param yAxisTitle       character string to be used as title of the y axis
#' @param xAxisTitle       character string to be used as title of the x axis
#' @return nothing
#' @export
plotInfiltrationDensity <- function(mDen, densityMarkers, bandWidth, clrs, plotTitle="", cellTypeLabels=NULL, sampleOrder=NULL,
                               yMax=NULL, separateLegend=FALSE, printLegend=TRUE, legendOnly=FALSE, yCol="Density",
                               pct=FALSE, facetByFOV=TRUE, facetByCellType=FALSE, yAxisTitle="Density (counts/mm^2)",
                               xAxisTitle="Distance to Tumor Interface"){

    selVars <- c("Sample","Band","BandLabels","CellTypeLabels",yCol) 
    if(facetByFOV){
        selVars <- c(selVars, "SPOT")
    }
    flog.debug("setting up labels")
    mDen$CellTypeLabels <- unlist(mDen$CellType)
    if(!is.null(cellTypeLabels)){
        for(x in names(cellTypeLabels)){
            mDen$CellTypeLabels[which(mDen$CellType == x)] <- cellTypeLabels[[x]]
        }
    }

    mDen$BandLabels <- factor(as.vector(sapply(as.vector(mDen$Band), function(x){
                         tmp <- unlist(strsplit(x,",",fixed=TRUE))
                         lbl <- as.numeric(gsub("\\(|\\]","",tmp[2]))
                         if(!is.null(bandWidth)){ lbl <- lbl - bandWidth/2 } ## shift half the width of a bar to make zero 
                                                                             ## the exact middle of plot
                     })))

    if(is.null(clrs)){
        clrs <- rep("grey",length(unique(mDen$CellTypeLabels)))
        names(clrs) <- unique(mDen$CellTypeLabels)
    } 

    ## spread density table to get to get a column for each cell type
    ## TO DO: MAKE SURE UPSTREAM THAT LABELS WILL BE UNIQUE; UNTIL THEN, JUST
    ## USE CELLTYPE AND ABANDON LABELS
    mt <- tryCatch({ mDen %>%
                     dplyr::select_(.dots = selVars) %>%
                     spread_("CellTypeLabels",yCol) },
              error = function(){ mDen$CellTypeLabels <- mDen$CellType
                                  mDen %>%
                                  dplyr::select_(.dots = selVars) %>%
                                  spread_("CellTypeLabels",yCol)
         })

    ## if faceting by FOV, mt includes SPOT but if not, it doesn't
    firstCol <- ifelse(facetByFOV, 5, 4)
    mtd <- gather(mt, firstCol:ncol(mt), key="CellType", value="Density")

    ## sort samples and FOV
    if(!is.null(sampleOrder)){
        flog.debug("Sorting samples in this order: %s", paste(sampleOrder,   collapse=", "))
        mtd$Sample2 <- factor(mtd$Sample, levels=sampleOrder)
    }
    if(facetByFOV){
        mtd$SPOT2 <- factor(mtd$SPOT, levels=sort(unique(as.vector(mtd$SPOT))))
    }

    if(!is.null(cellTypeLabels)){
        mtd$CellTypeLabels <- factor(mtd$CellType, levels=names(clrs))
    } else {
        mtd$CellTypeLabels <- factor(mtd$CellType)
    }

    ## finally, PLOT! 
    flog.debug("getting plot theme")
    plotTheme <- getPlotTheme("infiltration density")
    flog.debug("plotting")
    numBars <- length(unique(mtd$BandLabels))
    middle <- numBars/2
    p1 <- ggplot(mtd, aes(x=BandLabels, y=Density, fill=CellTypeLabels)) +
         geom_vline(xintercept=seq(floor(numBars/10),numBars,10)-0.5, color="white") +
         geom_vline(xintercept=middle+0.5, linetype="dotted", size=1) +
         geom_bar(stat="identity", position="stack", size=0.05, width=0.75) +
         plotTheme +
         scale_fill_manual(values=clrs, labels=names(clrs)) +
         scale_x_discrete(breaks=c(-305,-205,-105,-5,95,195,295),
                          labels=c("-300","-200","-100","0","100","200","300"),
                          expand = c(0.04, 0)) +
         ylab(yAxisTitle) +
         xlab(xAxisTitle) +
         labs(title=plotTitle) + 
         theme(axis.text.x = element_text(hjust=0), 
               panel.grid.major.x=element_blank()) 

    if(facetByFOV & facetByCellType){
        p1 <- p1 + facet_wrap(SPOT2 ~ CellType, ncol=3)
    } else if(facetByFOV){
        p1 <- p1 + facet_wrap(~SPOT2, ncol=3)
    } else if(facetByCellType){
        p1 <- p1 + facet_wrap(~CellTypeLabel, ncol=3)
    } else {
        p1 <- p1 + theme(axis.text.x = element_text(size=14))
    } 

    if(pct){ 
        p1 <- p1 + scale_y_continuous(labels = scales::percent, expand = c(0,0))
    } else {
        if(is.null(yMax)){
            yMax <- max(mtd$Density)*1.1
        }
        p1 <- p1 + scale_y_continuous(limits = c(0,yMax), expand=c(0,0))
    }

    return(p1)
}

#' Plot and print absolute density or density percentages of cell types in entire FOVs
#' 
#' Plot and print absolute density or density percentages of cell types in entire FOVs
#'
#' @param mDen             tibble containing density data 
#' @param densityMarkers   a vector of markers for which density is to be plotted
#' @param clrs             a list of colors, named by either by CellType from mDen tibble, 
#'                         or cellTypeLabels if specified
#' @param plotTitle        title of plot
#' @param cellTypeLabels   a list of labels for each marker being plotted, if different
#'                         than CellType in mDen tibble; default: NULL
#' @param sampleOrder      order in which samples should appear on the plot; default=NULL
#' @param yMax             if specified, y axis will span from zero to this number; default=NULL
#' @param separateLegend   logical; when TRUE, legend will be separated from plot. if printLegend is TRUE, 
#'                         will be printed first, before plot. if printLegend is FALSE, will not be printed at all; default=FALSE
#' @param printLegend      logical; when separateLegend is TRUE, if printLegend is FALSE, no legend will be printed at all; default=TRUE
#' @param legendOnly       logical; when separateLegend is TRUE, if legendOnly is TRUE, only legend will be printed, no plot
#' @param yCol             character string specifying which of the columns in mDen is to be used for y values; default="Density"
#' @param pct              logical; values being plotted are percentages; default=FALSE
#' @param yAxisTitle       character string to be used as title of the y axis
#' @param xAxisTitle       character string to be used as title of the x axis
#' @param annotationCols   vector of column names to use for annotation
#' @return plot 
#' @export
plotTotalFOVMarkerDensity <- function(mDen, densityMarkers, clrs, plotTitle="", cellTypeLabels=NULL, sampleOrder=NULL,
                               yMax=NULL, separateLegend=FALSE, printLegend=TRUE, legendOnly=FALSE, xCol="Sample", 
                               yCol="Total", pct=FALSE, groupBy=NULL, yAxisTitle="Density (counts/mm^2)", 
                               xAxisTitle="Sample", annotationCols=NULL){

    ####
    ## PREP DATA
    ####
    flog.debug("setting up labels")
    mDen$CellTypeLabels <- as.vector(unlist(mDen$CellType))
    if(!is.null(cellTypeLabels)){
        for(x in names(cellTypeLabels)){
            mDen$CellTypeLabels[which(mDen$CellType == x)] <- cellTypeLabels[[x]]
        }
    }
    if(!is.null(annotationCols) && any(grepl(" ",annotationCols))){
        for(x in grep(" ",annotationCols)){
            chars <- unlist(strsplit(annotationCols[x],""))
            if(chars[1] != "`" && chars[length(chars)] != "`"){
                annotationCols[x] <- paste0("`",annotationCols[x],"`") 
            }
        }
    }
    selVars <- paste0("`",names(mDen),"`")
    if(is.null(clrs)){
        clrs <- rep("grey",length(unique(mDen$CellTypeLabels)))
        names(clrs) <- unique(mDen$CellTypeLabels)
    }
 
    mt <- tryCatch({ mDen %>%
                     dplyr::select_(.dots = selVars)
                   }, error = function(err){ mDen$CellTypeLabels <- mDen$CellType
                                  mDen %>%
                                  dplyr::select_(.dots = selVars)
                  })
    if(!yCol == "Total"){
        mt$Total <- mt[[yCol]]
    }

    ## summarize by Sample, CellTypeLabels, and any group specified
    if(!is.null(groupBy)){
        annCols <- gsub("`","",annotationCols)
        groupByL <- tolower(groupBy)
        groupVars <- c("Sample","CellTypeLabels")
        addlGrpVars <- switch(groupByL,
                              "fov_type" = unique(c("SPOT","FOV_type", annCols)),
                              "patient"  = unique(c("Patient", annCols)),
                              "lesion response" = unique(c("Lesion Response", annCols)), 
                              flog.warn(paste0("Code does not currently support grouping/annotating by ",groupBy))
                       )
        groupVars <- unique(c(groupVars, addlGrpVars))

        mt <- mt %>% group_by_at(groupVars) %>% 
                     summarize(SampleTotal=sum(Total))

        mtd <- mt %>%
               select(Density=SampleTotal, everything())
    } else {
        mtd <- mt %>%
               select(Density=Total, everything())
    }
    
    ## factor all vars except density, making any of them ready for facetting
    for(nm in names(mtd)[-which(names(mtd) == "Density")]){
        mtd[[nm]] <- factor(mtd[[nm]])
    }

    ## set colors for cell type labels
    if(!is.null(cellTypeLabels)){
        mtd$CellTypeLabels <- factor(mtd$CellTypeLabels, levels=names(clrs))
    }

    ## sort samples and FOV
    if(!is.null(sampleOrder)){
        mtd$Sample <- factor(mtd$Sample, levels=sampleOrder)
    }

    ######
    ## finally, PLOT! 
    ######
    flog.debug("getting plot theme")
    plotTheme <- getPlotTheme("total FOV density")
    flog.debug("plotting")
    
    p1 <- ggplot(mtd, aes_string(x=xCol, y="Density", fill="CellTypeLabels")) +
          geom_bar(stat="identity", position="stack", size=0.5) +
          plotTheme +
          scale_fill_manual("Cell type", values=clrs, labels=names(clrs)) +
          scale_linetype_manual(values=c("solid","twodash","solid")) +
          scale_color_manual(values=c("white","black","black")) +
          ylab(yAxisTitle) +
          xlab(xAxisTitle) +
          labs(title=plotTitle, fill="Cell type") +
          theme(plot.margin=margin(10,0,10+30*length(annotationCols),60),
                axis.text.x = element_text(size=14), 
                strip.text.x = element_text(size=9, face="bold"),
                strip.background.x = element_rect(fill="#F8F8F8"))

    if(!is.null(groupBy)){
        p1 <- p1 + facet_grid(as.formula(paste0(". ~ `",groupBy, "`")), scales="free_x", space="free", switch="both")
    }

    if(pct){
        p1 <- p1 + scale_y_continuous(labels = scales::percent, expand = c(0,0))
    } else {
        if(is.null(yMax)){
            yMax <- mtd %>% group_by_at(groupBy) %>% summarise(TotDen=sum(Density)) %>% pull(TotDen) %>% max()*1.1
        }
        p1 <- p1 + scale_y_continuous(limits = c(0,yMax), expand=c(0,0))
    } 

    gt <- NULL
    if(!is.null(annotationCols)){

        p1 <- p1 + theme(axis.text.x = element_blank())

        mtdL <- mtd[1,]
        mtdL[[groupBy]] <- ""
        mtdL[["lbl"]][1] <- groupBy
        
        size <- 3.5
        y <- yMax *.025
        p1 <- p1 + geom_text(data=mtdL, aes(x=-1.2, y=y, label=lbl), vjust=3.5, hjust=1, size=size, fontface="bold")

        size <- 3.5
        for(x in seq(annotationCols)){
            mtdL[["lbl"]] <- ""
            mtdL[["lbl"]][1] <- gsub("`","",annotationCols[x])
            vjust <- 3.5 + 2*x
            annCol <- ifelse((grepl(" ",annotationCols[x]) && !grepl("`",annotationCols[x])), paste0("`",annotationCols[x],"`"), annotationCols[x])
            p1 <- p1 + geom_text(data=mtd, aes_string(x="Sample", y=y, label=annCol), vjust=vjust, hjust=0.5, size=size)
            p1 <- p1 + geom_text(data=mtdL, aes(x=-1.2, y=y, label=lbl), vjust=vjust, hjust=1, size=size, fontface="bold")
        }
    
        gt <- ggplot_gtable(ggplot_build(p1))

        ## change the color of the NA strip to white 
        strps <- which(grepl("strip-",gt$layout$name))
        fills <- rep("#F8F8F8", length(strps))
        fills[1] <- "white"

        k <- 1
        for (i in strps) {
            j <- which(grepl('rect', gt$grobs[[i]]$grobs[[1]]$childrenOrder))
            gt$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
            k <- k+1
        }
        
        # Code to override clipping
        gt$layout$clip[grep("panel",gt$layout$name)] <- "off"
        #grid.draw(gt)

    } else {
        gt <- ggplot_gtable(ggplot_build(p1))
    }
    return(gt)
    #return(p1)
}

#' Box plots of densities showing differences between samples for different cell types
#'
#' Box plots of densities showing differences between samples for different cell types
#
#' @param markers     vector of markers, for each of which a plot will be generated
#' @param density     tibble of density data, including columns for at least Sample, CellType
#'                    and [Total|Density]
#' @param sampleOrder if specified, sample boxes will appear in this order; default=NULL
#' @param starSize    text size, specifically for significance starts; default=8
#' @param pdfFile     if specified, plots will be printed to this file; default=NULL 
#' @return nothing
densityBoxPlots <- function(markers, density, sampleOrder=NULL, starSize=8, pdfFile=NULL){

    plotTheme = theme(legend.title=element_blank(),
                     axis.text=element_text(size=8),
                     axis.ticks = element_blank(),
                     strip.text.x = element_text(face="bold", size=12, margin = margin(.4, .2, .4, .2, "cm")),
                     strip.text.y = element_text(margin = margin(.4, .2, .4, .2, "cm")),
                     strip.placement = "outside",
                     plot.background = element_blank(),
                     panel.background = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.line = element_line(color = "black", size=0.5),
                     )

    ## density may have "Total" or "Density" as column name; make sure it's Density
    if(!"Density" %in% names(density)){
        if(!"Total" %in% names(density)){
            warning("Data passed to densityBoxPlot() does not contain either 'Density' or 'Total' columns. Skipping")
            return()
        }
        density <- rename(density, Density=Total)
    }


    if(!is.null(pdfFile)){
        pdf(pdfFile,height=8.5,width=11)
    }
    for(marker in markers){
        if(!(marker %in% density$CellType)){
            warning(paste0("Marker ",marker," not in data"))
        }
        den <- filter(density, CellType == marker)
        if(!is.null(sampleOrder)){
             den$Sample <- factor(den$Sample, levels=sampleOrder)
        }
        ## we want to compare means of all pairwise combinations
        comparisons <- as.tibble(combn(as.vector(unique(den$Sample)),m=2))
        g <- ggplot(den, aes(Sample, Density)) +
           geom_boxplot(width=0.5) +
           stat_compare_means(comparisons=comparisons,label="p.signif",hide.ns=TRUE) +
           plotTheme +
           ggtitle(marker) +
           xlab("") +
           ylab("")
        ## change size of significance stars
        g$layers[[2]]$aes_params$textsize <- starSize
        print(g)
    }
    if(!is.null(pdfFile)){
        dev.off()
    }

}


#' Generate all infiltration density plots for a single sample
#' 
#' Generate all infiltration density plots for a single sample
#'
#' @param den                tibble with columns for Sample,SPOT,Band,Counts,Density
#' @param config             parsed marker config (return value from getMarkerConfig())
#' @param markers            vector of marker combinations to plot
#' @param bandWidth          width of each band around tumor interface
#' @param sampleOrder        vector of sample names in order they should appear on plots
#' @param plotTitle          title to appear on plot
#' @param yMax               maximum y value for all plots
#' @param byFOV              logical indicating whether to generate plots for each individual FOV (default=TRUE)
#' @param absoluteDensity    logical indicating whether to plot density values on absolute scale (default=TRUE)
#' @param densityPercentage  logical indicating whether to plot density percentages (default=TRUE)
#' @param summarize          logical indicating whether to generate summary plots of total band densities (ALL FOVs; default=TRUE)
#' @param stacked            logical indicating whether to generate plots of all cell types for a single marker set
#'                           in a stacked bar chart (default=TRUE)
#' @param summaryTitle       title to appear on summary plots 
#' @return list of sample plots
#' @export
makeSampleInfiltrationDensityPlots <- function(den, config, markers, bandWidth, 
                                               sampleOrder=NULL, plotTitle=NULL, yMax=NULL,sumYmax=NULL, 
                                               byFOV=TRUE, absoluteDensity=TRUE, 
                                               densityPercentage=TRUE, summarize=TRUE, 
                                               summaryTitle=NULL, forceStack=FALSE){
    plotList <- list()

    ## set up labels and colors based on config
    ctLabels <- config$Label
    names(ctLabels) <- config$CellType
    clrLbls <- unique(config[,c("Color","Label")])
    ctClrs <- pull(clrLbls[,"Color"])
    names(ctClrs) <- pull(clrLbls[,"Label"])

    ## first, plot individual FOV
    tryCatch({

        if(absoluteDensity){
            flog.info("      plotting absolute infiltration density")
            if(is.null(yMax)){ yMax <- max(den$Density) * 1.25 }
            plotList[["sep_fov"]][["absolute"]] <- plotInfiltrationDensity(den, markers, bandWidth, ctClrs, 
                                                              sampleOrder=sampleOrder,
                                                              plotTitle=plotTitle, yMax=yMax, separateLegend=FALSE,
                                                              facetByFOV=TRUE, facetByCellType=FALSE, 
                                                              cellTypeLabels=ctLabels)
        }
        if(densityPercentage && length(unique(den$CellType)) > 1){
            if(length(which(den$Density > 0)) > 1){
                ctToStack <- stackable(unique(den$CellType))
                if(!is.null(ctToStack) && length(ctToStack) > 1){
                    den <- den %>% filter(CellType %in% ctToStack)
                    flog.info("      plotting infiltration density percentages")
                    tmp <- den %>% select(Sample,SPOT,Band,CellType,Counts) %>%
                           group_by(Sample,SPOT,Band) %>% summarise(TotalCountsPerBandPerFOV=sum(Counts)) %>% ungroup()
                    sDenPct <- den %>% left_join(tmp, by=c("Sample","SPOT","Band")) %>%
                               mutate(Percent=Counts/TotalCountsPerBandPerFOV)      
                    sDenPct$Percent[is.na(sDenPct$Percent)] <- 0
                    plotList[["sep_fov"]][["percentage"]] <- plotInfiltrationDensity(sDenPct, markers, 
                                                           bandWidth, ctClrs, sampleOrder=sampleOrder,
                                                           plotTitle=plotTitle, yMax=yMax, separateLegend=FALSE,
                                                           printLegend=TRUE, yCol="Percent", facetByFOV=TRUE,
                                                           facetByCellType=FALSE, yAxisTitle="Percent Total Density",
                                                           cellTypeLabels=ctLabels,pct=TRUE)
                }
            }
        }
      }, error = function(){
            warning(paste0("Could not plot sample ",s,", for cell population ",pop))
      }
    )

    ## then, summarize each band over all FOV
    if(summarize){
        tryCatch({
          if(!"BandArea" %in% names(den) && "Area" %in% names(den)){
              den <- den %>% select(BandArea=Area, everything())
          }

          ssum <- den %>% 
                  mutate(Density=SampleBandCellTypeDensity) %>% 
                  select(Sample,Band,CellType,SampleBandCellTypeCounts,Density) %>% 
                  unique()

          ssum[is.na(ssum)] <- 0
          if(is.null(sumYmax)){
              sumYmax <- max(ssum$Density) * 1.25
          }
          if(length(which(!is.na(ssum$Density))) == 0){
              flog.info(paste0("Could not plot summary for sample ",s,", cell population ",pop))
          } else {
              flog.info("      plotting infiltration summary")
              if(absoluteDensity){
                  plotList[["summary"]][["absolute"]] <- plotInfiltrationDensity(ssum, densityMarkers,
                                                                             bandWidth, ctClrs, plotTitle=plotTitle, yMax=sumYmax,
                                                                             separateLegend=FALSE, printLegend=TRUE, facetByFOV=FALSE,
                                                                             facetByCellType=FALSE, yCol="Density",cellTypeLabels=ctLabels)
              }
              if(densityPercentage && length(unique(ssum$CellType)) > 1){
                  ctToStack <- stackable(unique(ssum$CellType))
                  if(!is.null(ctToStack) && length(ctToStack) > 1){
                      flog.info("      plotting infiltration summary percentages")
                      tmp <- ssum %>% group_by(Sample,Band) %>% 
                                      summarise(TotalSampleBandCounts=sum(SampleBandCellTypeCounts)) %>% 
                                      ungroup()
                      ssumPct <- ssum %>% left_join(tmp, by=c("Sample","Band")) %>%
                                 mutate(Percent=SampleBandCellTypeCounts/TotalSampleBandCounts) %>% unique()
                      plotList[["summary"]][["percentage"]] <- plotInfiltrationDensity(ssumPct, densityMarkers,
                                                                                       bandWidth, ctClrs, sampleOrder=sampleOrder,
                                                                                       plotTitle=plotTitle, yMax=sumYmax, separateLegend=FALSE,
                                                                                       printLegend=TRUE, yCol="Percent", facetByFOV=FALSE,
                                                                                       facetByCellType=FALSE, yAxisTitle="Percent Total Density",
                                                                                       cellTypeLabels=ctLabels,pct=TRUE)
                  }
              }
          }
         }, error = function(){
               warning(paste0("Could not plot summary for sample ",s,", cell population ",pop))
        })
    } # end if summarize
    return(plotList)     
}

#' Get maximum density value from all Sample+Band combos in a data set
#' 
#' Get maximum density value from all Sample+Band combos in a data set
#' 
#' @param den       density tibble
#' @param markers   markers to consider
#' @param ia        infiltration area tibble (Sample,Band,Area)
#' @return  double; maximum density value from all Sample+band groups
#' @export
getSampleBandMaxY <- function(den, markers, ia){
    sDen <- den %>% filter(CellType %in% markers) %>%
            group_by(Sample,Band) %>%
            summarise(SampleBandCounts=sum(Counts))
    sDen <- sDen %>% full_join(ia, by=c("Sample","Band")) %>%
            mutate(SampleBandDensity=SampleBandCounts/TotalSampleBandArea)
    return(max(sDen$SampleBandDensity))         
}


#' Determine whether a set of markers is 'stackable' in a bar chart
#'
#' Extract from a set of markers only those combinations that are mutually exlusive and therefore 
#' 'stackable' in a bar chart
#' 
#' @param markers  vector of marker combinations to be tested for stackability
#' @return subset of markers including only the ones that ARE stackable 
#' @export
stackable <- function(markers){

    pos <- lapply(markers, function(x){ unlist(strsplit(x,","))[-grep("-",unlist(strsplit(x,",")))] })
    posMarkers <- unique(unlist(pos))

    neg <- lapply(markers, function(x){ unlist(strsplit(x,","))[grep("-",unlist(strsplit(x,",")))] })
    negMarkers <- unique(unlist(neg))

    unstackable <- c()

    for(m in posMarkers){

        ps <- pos[grep(m,pos)]
        ns <- neg[grep(m,pos)]
 
        for(x in 1:(length(ps))){
            combo <- ps[[x]]
            comboNegs <- gsub("-","",ns[[x]])
            others <- ps[-x]
            if(length(others) == 0){ next }
            othersNegs <- ns[-x]
            for(o in 1:length(others)){
                addlMarkers <- others[[o]][-which(others[[o]] %in% combo)]
                if(!all(addlMarkers %in% comboNegs)){
                    warning(paste0("Cell types [",paste(c(ps[[x]],ns[[x]]),collapse=","), 
                                   "] and [",   
                                   paste(c(others[[o]],othersNegs[[o]]),collapse=","), 
                                   "] are not stackable"))
                    unstackable <- c(unstackable, paste(c(ps[[x]],ns[[x]]),collapse=","))
                }
            }
        }
    }

    ## for any markers that are never positive, make sure they are negative
    ## in ALL populations
    for(m in negMarkers){
        if(gsub("-","",m) %in% posMarkers){ next }
        for(x in 1:length(ns)){
            if(!m %in% ns[[x]]){
                warning(paste0("Marker ",m," is negative in only some combos and is never positive. Must be negative in ALL. Excluding all combos with this marker from stacking."))
                unstackable <- c(unstackable, markers[grep(m,markers)])
            }
        }
    }

    if(!is.null(unstackable) && length(unstackable) > 0){
        return(markers[-which(markers %in% unstackable)])
    } else {
        return(markers)
    }
}

#' Create and store all infiltration plots in a list to be printed later
#'
#' Create and store all infiltration plots in a list to be printed later
#'
#' @param den                tibble with columns for Sample,SPOT,Band,Counts,Density
#' @param bandWidth          width of each band around tumor interface
#' @param config             parsed marker config (return value from getMarkerConfig())
#' @param yScaleConsistency  possible values: "population", "all", or "sample" indicating which
#'                           subset of values to use to set maximum for all plots
#' @param indivPopulations   plot populations individually; default=TRUE; if set to FALSE, ONLY stacked plots
#'                           will be created
#' @param absoluteDensity    logical indicating whether to plot density values on absolute scale (default=TRUE)
#' @param densityPercentage  logical indicating whether to plot density percentages (default=TRUE)
#' @param byFOV              logical indicating whether to generate plots for each individual FOV (default=TRUE)
#' @param summarize          logical indicating whether to generate summary plots of total band densities (ALL FOVs; default=TRUE)
#' @param stacked            logical indicating whether to generate plots of all cell types for a single marker set
#'                           in a stacked bar chart (default=TRUE)
#' @param sampleOrder        vector of sample names in order they should appear on plots
#' @param infiltrationAreas  area data for each sample/spot/band
#' @return list of all plots
#' @export
getAllInfiltrationPlots <- function(den, bandWidth=NULL, config, yScaleConsistency="population", absoluteDensity=TRUE,
                              densityPercentage=TRUE, byFOV=TRUE, summarize=TRUE, stacked=TRUE,
                              sampleOrder=NULL, infiltrationAreas=NULL, indivPopulations=TRUE, forceStack=FALSE){

    allPlots <- list()

    ## get sample/band level density info
    den <- getAllInfiltrationDensityValues(den, infiltrationAreas)

    if(!"Density" %in% names(den) && "Total" %in% names(den)){
        den <- den %>% select(Density=Total, everything())
    }

    yMax <- NULL
    if(yScaleConsistency=="all"){
        yMax <- max(den$Density)
    }
    for(ms in unique(config$MarkerSet)){
        flog.info(paste0("MARKER SET: ",ms))

        msCfg <- config %>% filter(MarkerSet == ms)
        densityMarkers <- unique(msCfg$CellType)
        msDen <- den %>% filter(CellType %in% densityMarkers)
        msYmax <- msDen %>% 
                  group_by(Sample,Band) %>% 
                  mutate(TotalSampBandDensity=sum(Counts)/SampleBandArea) %>% 
                  ungroup() %>% 
                  summarise(max(TotalSampBandDensity)) %>% 
                  pull()

        if(indivPopulations){
            for(pop in unique(msCfg$Population)){
                flog.info(paste0("  POPULATION: ",pop))

                popCfg <- msCfg %>% filter(Population == pop)
                ## extract density for this population
                densityMarkers <- unique(popCfg$CellType)
                popDen <- msDen %>% filter(CellType %in% densityMarkers)
                popYmax <- popDen %>% 
                           group_by(Sample,Band) %>% 
                           mutate(TotalSampBandDensity=sum(Counts)/SampleBandArea) %>% 
                           ungroup() %>% 
                           summarise(max(TotalSampBandDensity)) %>% 
                           pull()
                for(s in unique(popDen$Sample)){
                    flog.info(paste0("    SAMPLE: ",s))
                    sDen <- popDen %>% filter(Sample == s) %>% ungroup()
                    plotTitle <- paste0(s, " Interface FOVs:\n",pop)
                    summaryTitle <- paste0(s, " Interface All FOVs (total density):\n",pop) 

                    allPlots[[ms]][[pop]][[s]] <- makeSampleInfiltrationDensityPlots(sDen, popCfg, densityMarkers, bandWidth,
                                                                              sampleOrder=sampleOrder,
                                                                              plotTitle=plotTitle, yMax=yMax, sumYmax=popYmax, 
                                                                              absoluteDensity=absoluteDensity, 
                                                                              densityPercentage=densityPercentage,
                                                                              summarize=summarize, summaryTitle=summaryTitle,
                                                                              forceStack=forceStack)
                } #end each sample
            } #end each population
        } # end if indivPops

        ### if all negatives are set, making densities mutually exclusive, we can stack
        ### all cell types for a single marker set
        if(stacked && (length(unique(msCfg$Label)) == nrow(msCfg))){
            ctToStack <- stackable(unique(msCfg$CellType))
            if(is.null(ctToStack) || length(ctToStack) == 0){ 
                if(forceStack){ 
                    flog.warn("Cell types do not appear to be stackable, but forcing stacking. Plot MAY not be accurate.")
                    densityMarkers <- unique(msCfg$CellType)
                } else {
                    flog.warn("Cell types are not stackable. No stacked plot made")
                    next
                }
            } else {
                ## get stacked version
                densityMarkers <- ctToStack
            }
 
            tDen <- msDen %>% filter(CellType %in% densityMarkers) %>% ungroup()

            for(s in unique(tDen$Sample)){
                flog.debug(paste0("  ",s))
                sDen <- tDen %>% filter(Sample == s)
                #sDen <- getAllInfiltrationDensityValues(sDen, ia)
                plotTitle <- paste0(s, " Interface FOVs: ",unique(msCfg$MarkerSetAlias))
                summaryTitle <- paste0(s, " Interface All FOVs (total density): ",unique(msCfg$MarkerSetAlias)) 
                
                allPlots[[ms]][["stacked"]][[s]] <- makeSampleInfiltrationDensityPlots(sDen, msCfg, densityMarkers, bandWidth,
                                                                              sampleOrder=sampleOrder,
                                                                              plotTitle=plotTitle, yMax=yMax, sumYmax=msYmax, 
                                                                              absoluteDensity=absoluteDensity, 
                                                                              densityPercentage=densityPercentage,
                                                                              summarize=summarize, summaryTitle=summaryTitle, 
                                                                              forceStack=forceStack)
            } # end for each sample
        } # end if stacked
    } # end for each marker set
    return(allPlots)
}


#' Organize and print multiple sample summaries on a page or two
#' 
#' For each marker set in list of plots organized by markerSet, population, then sample,
#' gather all sample plots for a single population, organizing with each row including
#' a single patient and each column including a treatment
#'
#' @param den            density table including Sample_name, Patient, and Response in addition to standard density data
#'                       NOTE: currently Sample_name MUST be in the form "Pt[X]_[ResponseX]", e.g., "Pt1_UT1" or "Pt1_PR"
#' @param allPlots       list of plots organized by markerSet, population, sample, [summary], [absolute|percentage]
#' @param numPtsPerPage  number of cases to plot on each page
#' @return nothing
#' @export
sampleInfiltrationPlotsSingleFile <- function(den, allPlots, numPtsPerPage=4, ncol=NULL, nrow=NULL){

    ptInfo <- den %>%
           ungroup() %>% 
           select(Sample,Patient,`Lesion Response`) %>% 
           unique() %>% 
           arrange(Patient,desc(`Lesion Response`))

    nColEachResp <- ptInfo %>% 
                    group_by(Patient,`Lesion Response`) %>% 
                    summarize(NumEach=n()) %>% 
                    ungroup() %>% 
                    group_by(`Lesion Response`) %>% 
                    summarize(MaxEach=max(NumEach)) %>% 
                    ungroup() %>% 
                    arrange(desc(`Lesion Response`))
 
    ncol <- sum(nColEachResp$MaxEach) 
    nColEachResp <- nColEachResp %>% mutate(FirstCol=cumsum(MaxEach)+1 - MaxEach[1])

    if(is.null(nrow)){
        nrow <- length(unique(ptInfo$Patient))
    } 
    numPages <- ceiling(length(unique(ptInfo$Patient))/nrow)


    lgnd <- title <- xtitle <- ytitle <- NULL
 
    for(pltType in c("absolute","percentage")){

        tmp <- list()
        pgs <- list()

        for(x in 1:numPages){
            tmp[[x]] <- list()
            pgs[[x]] <- matrix(NA, ncol=ncol, nrow=nrow)
            cnames <- c()
            for(y in 1:length(nColEachResp$`Lesion Response`)){
                for(z in 1:nColEachResp$MaxEach[y]){
                    cnames <- c(cnames, paste0(nColEachResp$`Lesion Response`[y],z))
                }
            }
            colnames(pgs[[x]]) <- cnames
            rownames(pgs[[x]]) <- paste0("Pt_",sort(unique(ptInfo$Patient)))      
        }

        print(pltType)
        for(x in seq(names(allPlots))){
            sl <- names(allPlots)[x]
            print(sl)
 
            if(!pltType %in% names(allPlots[[sl]][["summary"]])){ next }

            plt <- allPlots[[sl]][["summary"]][[pltType]] +
                   theme(legend.text=element_text(size=rel(0.75)),
                         legend.key.size=unit(0.4,"cm"),
                         legend.background=element_rect(color="gray", size=0.5))

            if(is.null(plt)){ next }

            if(is.null(lgnd)){
                tryCatch({ lgnd <- g_legend(plt) },error = function(e){ warning(e) })
                  title <- paste0(gsub(".*? (.+)","\\1",plt$labels$title),"\n")
                  xtitle <- plt$labels$x
                  ytitle <- plt$labels$y
            }

            plt <- plt + theme(legend.position="none",
                               plot.title=element_text(size=10),
                               plot.margin=unit(c(0,-.25,-.25,.25),"cm"),
                               axis.text.x=element_blank(),
                               axis.text.y=element_text(size=8))
            plt$labels$title <- NULL
            plt$labels$x <- plt$labels$y <- ""

            sname <- unlist(strsplit(sl,"\\."))[1]
            resp <- ptInfo %>% filter(Sample == sl) %>% pull(`Lesion Response`)
            plRow <- sname

            ## find right layoutMatrix
            page <- 1
            for(pg in seq(pgs)){
                if(plRow %in% rownames(pgs[[pg]])){
                    page <- pg
                    break
                }
            }

            plCol <- paste0(resp, min(which(is.na(pgs[[pg]][sname,grep(resp, colnames(pgs[[pg]]))]))))

            tmp[[page]][[length(tmp[[page]])+1]] <- plt
            pgs[[pg]][plRow,plCol] <- length(tmp[[pg]])
        } # end for each samp

        if(!is.null(unlist(tmp))){
            tmp[[numPages]][[length(tmp[[numPages]])+1]] <- lgnd

            ### if there are no empty (NA) spaces in matrix, add another column
            ### for legend
            lastPg <- pgs[[numPages]]
            if(!(is.na(lastPg[nrow(lastPg),ncol(lastPg)]) || is.na(lastPg[1,ncol(lastPg)]))){
                lgndCol <- rep(length(tmp[[numPages]]),nrow(pgs[[numPages]]))
                pgs[[numPages]] <- cbind(pgs[[numPages]],lgndCol)
            } else {
                lastCol <- lastPg[,ncol(lastPg)]
                ## find longest stretch of NA at bottom or top of last column and fit legend there
                if(is.na(lastPg[nrow(lastPg),ncol(lastPg)])){
                    x = nrow(lastPg)
                    inc <- -1
                } else {
                    x <- 1
                    inc <- 1
                }
                lastCol[x] <- length(tmp[[numPages]])
                while((x + inc) %in% seq(1,nrow) && is.na(lastCol[x + inc])){
                    lastCol[x + inc] <- length(tmp[[numPages]])
                    x <- x + inc
                }
                pgs[[numPages]][,ncol(pgs[[numPages]])] <- lastCol 
           }
            for(pg in seq(pgs)){
                print(pg)
                grid.arrange(grobs=tmp[[pg]], 
                             layout_matrix=pgs[[pg]], 
                             bottom=xtitle, 
                             left=ytitle, 
                             top=textGrob(title, gp=gpar(fontsize=12, lineheight=1)))
            }
        }
    } # end each pltType
    return()
}

#' Print density plots 
#'
#' Print density plots
#'
#' @param den                tibble containing densities for all markers to be plotted
#' @param areaDat            tibble containing area data for all Samples, FOVs and Bands
#' @param config             tibble containing marker configureation (returned from getMarkerConfig())
#' @param bandWidth          if plotting by distance intervals from infiltration boundary, this is the size of 
#'                           each interval
#' @param yScaleConsistency  level on which y-scales should be the same; choices: ["byPopulation"|"all"|"markerSet"]
#' @param indivPopulations   logical; print a plot for each individual population; default=TRUE
#' @param absoluteDensity    logical; print plots of absolute density values; default=TRUE
#' @param densityPercentage  logical; print plots of density percentage values; default=TRUE
#' @param byFOV              logical; print plots of each FOV (all FOV on one panel); default=TRUE
#' @param summarize          logical; print summary plots, showing density of ALL FOV together; default=TRUE
#' @param stacked            logical; for plots with all exclusive markers, print stacked bar plots; default=TRUE
#' @param forceStack         force stacking of cell types even if they appear to be unstackable
#' @param sampleOrder        vector of samples in order that they should appear on plot
#' @param outDir             directory where all plots will be printed, organized by marker set
#' @param sampleSummaryRows  number of rows to include on each page of sample summary file; default=NULL
#' @param sampleSummaryCols  number of cols to include on each page of sample summary file; default=NULL
#' @param sampleSummaryPtsPerPage  number of patients to include on each page of sample summary file; can be used instead
#'                                 of sampleSummary[Rows|Cols]; default=NULL
#' @return nothing 
#' @export
printInfiltrationDensityPlots <- function(den, bandWidth=NULL, config, yScaleConsistency="population", 
                                     indivPopulations=TRUE, absoluteDensity=TRUE, densityPercentage=TRUE, byFOV=TRUE, 
                                     summarize=TRUE, stacked=TRUE, forceStack=FALSE, sampleOrder=NULL, separateLegend=TRUE, 
                                     infiltrationAreas=NULL, outDir=NULL, sampleSummaryRows=NULL, sampleSummaryCols=NULL, 
                                     sampleSummaryPtsPerPage=NULL){

    if(is.null(outDir)){
        outDir <- getwd()
    }

    flog.info("Getting all plots")
    allPlots <- getAllInfiltrationPlots(den, config, bandWidth=bandWidth, indivPopulations=indivPopulations, 
                                        absoluteDensity=absoluteDensity, densityPercentage=densityPercentage, 
                                        byFOV=byFOV, summarize=summarize, stacked=stacked, 
                                        infiltrationAreas=infiltrationAreas, forceStack=forceStack) 

    if(is.null(allPlots) || length(allPlots) == 0){
        flog.info("No plots made.")
        return()
    }

    for(ms in names(allPlots)){

        if(is.null(allPlots[[ms]]) || length(allPlots[[ms]]) == 0){
            flog.info(paste0("No plots made for marker set ",ms))
            next
        }
        ## create directories 
        msOutDir <- file.path(outDir,ms)
        dir.create(msOutDir, showWarnings=TRUE, recursive=TRUE)

        fovOutDir <- file.path(msOutDir, "individual_populations_individual_fovs")
        ctSumOutDir <- file.path(msOutDir, "individual_populations_whole_samples")
        stackedOutDir <- file.path(msOutDir, "stacked_populations_individual_fovs")
        stackedSumOutDir <- file.path(msOutDir, "stacked_population_whole_samples")
        sampleOutDir <- file.path(msOutDir, "individual_populations_sample_summary")
        sampleSumOutDir <- file.path(msOutDir, "stacked_populations_sample_summary")

        flog.info(paste0("Printing plots for marker set: ",ms))

        for(ct in names(allPlots[[ms]])){
            flog.info(paste0("  CellType/Population: ",ct))

            if(length(allPlots[[ms]][[ct]]) == 0){
                flog.info(paste0("No plots made for marker set ",ms," cell type ",ct))
                next
            }

            ## print sample plots separately (single page each)
            for(s in names(allPlots[[ms]][[ct]])){
                if(length(allPlots[[ms]][[ct]][[s]]) == 0){
                    flog.info(paste0("No plots made for marker set ",ms," cell type ",ct," sample ",s))
                    next
                }
                flog.info(paste0("    Sample: ",s))

                p1 <- allPlots[[ms]][[ct]][[s]][[1]][[1]]
                lgnd <- g_legend(p1) ## automatically get legend from first plot
  
                for(plotType in names(allPlots[[ms]][[ct]][[s]])){
                    if(ct == "stacked"){
                        if(plotType == "sep_fov"){
                            if(!file.exists(stackedOutDir)){ dir.create(stackedOutDir, showWarnings=FALSE, recursive=TRUE) }
                            fileName <- file.path(stackedOutDir, paste0(ms,"_",s,"_by_fov.pdf"))
                        } else {
                            if(!file.exists(stackedSumOutDir)){ dir.create(stackedSumOutDir, showWarnings=FALSE, recursive=TRUE) }
                            fileName <- file.path(stackedSumOutDir, paste0(ms,"_",s,".pdf"))
                        }
                   } else {
                        if(plotType == "sep_fov"){
                            if(!file.exists(fovOutDir)){ dir.create(fovOutDir, showWarnings=FALSE, recursive=TRUE) }
                            fileName <- file.path(fovOutDir, paste0(ct,"_",s,"_by_fov.pdf"))
                        } else {
                            if(!file.exists(ctSumOutDir)){ dir.create(ctSumOutDir, showWarnings=FALSE, recursive=TRUE) }
                            fileName <- file.path(ctSumOutDir, paste0(ct,"_",s,".pdf"))
                        }
                    }
                    plts <- allPlots[[ms]][[ct]][[s]][[plotType]]
                    if(is.null(lgnd) || is.null(plts)){ next }
                    pdf(fileName, height=8.5, width=11)
                    flog.debug("drawing legend")
                    if(!is.null(lgnd)) { grid.draw(lgnd) }
                    for(p in names(plts)){ 
                        p1 <- plts[[p]] + theme(legend.position="none")
                        if(!is.null(p1)){ print(p1) }
                    }
                    dev.off()
                }
            } 

            ## print summaries for each sample on one page
            if(!ct == "stacked"){
                if(!file.exists(sampleOutDir)){ dir.create(sampleOutDir, showWarnings=FALSE, recursive=TRUE) }
                flog.info(paste0("  Printing single-page summary for cell type: ",ct))
                fileName <- file.path(sampleOutDir, paste0(ct, "__combined_summary.pdf"))
            } else {
                if(!file.exists(sampleSumOutDir)){ dir.create(sampleSumOutDir, showWarnings=FALSE, recursive=TRUE) }
                flog.info(paste0("Printing single-page summary for marker set: ",ms))
                fileName <- file.path(sampleSumOutDir, paste0(ms, "__combined_summary.pdf"))
            }
            pdf(fileName, height=11, width=8.5)
            sampleInfiltrationPlotsSingleFile(den, allPlots[[ms]][[ct]], numPtsPerPage=sampleSummaryPtsPerPage, 
                                              ncol=sampleSummaryCols, nrow=sampleSummaryRows)
            dev.off()
        }
    }
}

#' Generate all infiltration density plots for a single sample
#' 
#' Generate all infiltration density plots for a single sample
#'
#' @param den                tibble with columns for Sample,SPOT,Band,Counts,Density
#' @param config             parsed marker config (return value from getMarkerConfig())
#' @param densityMarkers     vector of marker combinations to plot
#' @param sampleOrder        vector of sample names in order they should appear on plots
#' @param plotTitle          title to appear on plot
#' @param yMax               maximum y value for all plots
#' @param sumYmax            maximum y value for summary plots
#' @param byFOV              logical indicating whether to generate plots for each individual FOV (default=TRUE)
#' @param absoluteDensity    logical indicating whether to plot density values on absolute scale (default=TRUE)
#' @param densityPercentage  logical indicating whether to plot density percentages (default=TRUE)
#' @param summarize          logical indicating whether to generate summary plots of total band densities (ALL FOVs; default=TRUE)
#' @param stacked            logical indicating whether to generate plots of all cell types for a single marker set
#'                           in a stacked bar chart (default=TRUE)
#' @param forceStack         force stacking of cell types even if they appear to be unstackable
#' @param summaryTitle       title to appear on summary plots 
#' @return list of sample plots
#' @export
makeSampleTotalFOVDensityPlots <- function(den, config, densityMarkers, sampleOrder=NULL,
                                           plotTitle="", yMax=NULL, sumYmax=NULL, byFOV=TRUE,
                                           absoluteDensity=TRUE, densityPercentage=TRUE,
                                           summarize=TRUE, summaryTitle="", forceStack=FALSE, 
                                           annotationCols=NULL){

    plotList <- list()

    ## set up labels and colors based on config
    ctLabels <- config$Label
    names(ctLabels) <- config$CellType
    clrLbls <- unique(config[,c("Color","Label")])
    ctClrs <- pull(clrLbls[,"Color"])
    names(ctClrs) <- pull(clrLbls[,"Label"])

    ## first, plot individual FOV
    if(byFOV){
        tryCatch({
            if(absoluteDensity){
                flog.info("      plotting absolute total FOV density")
                plotList[["sep_fov"]][["absolute"]] <- plotTotalFOVMarkerDensity(den, densityMarkers, ctClrs,
                                                              sampleOrder=sampleOrder, xCol="SPOT", yCol="Total",
                                                              plotTitle=plotTitle, yMax=yMax, separateLegend=FALSE,
                                                              groupBy="FOV_type", cellTypeLabels=ctLabels, xAxisTitle="FOV")
            }
            if(densityPercentage && length(unique(den$CellType)) > 1){
                ctToStack <- stackable(unique(den$CellType))
                if(!is.null(ctToStack) && length(ctToStack) > 1){
                    den <- den %>% filter(CellType %in% ctToStack)
                    if(length(which(den$Total > 0)) > 1){
                        flog.info("      plotting total FOV density percentages")
                        sDenPct <- den %>% group_by(Sample,SPOT) %>% mutate(Pct=Counts/sum(Counts))                

                        sDenPct$Pct[is.na(sDenPct$Pct)] <- 0
                        plotList[["sep_fov"]][["percentage"]] <- plotTotalFOVMarkerDensity(sDenPct, markers,
                                                               ctClrs, sampleOrder=sampleOrder, xCol="SPOT",
                                                               plotTitle=plotTitle, yMax=yMax, separateLegend=FALSE,
                                                               printLegend=TRUE, yCol="Pct", groupBy="FOV_type", 
                                                               yAxisTitle="Percent Total FOV Density", xAxisTitle="FOV",
                                                               cellTypeLabels=ctLabels, pct=TRUE)
                    }
                }
            }
          }, error = function(){
                warning(paste0("Could not plot sample ",s,", for cell population ",pop))
          }
        )
    }

    ### summarize each band over all FOV
    if(summarize){
        plotTitle <- gsub(".*\\s","",summaryTitle) 
        tryCatch({
          if(!is.null(annotationCols)){
              ssum <- den %>% 
                      select_at(c("CellType","Sample","SampleCounts","CellTypeSampleDensity","Patient",annotationCols)) %>% 
                      unique()
          } else {
              ssum <- den %>% select(CellType,matches("Sample"),Patient) %>% unique()
          }
          if(!is.null(sampleOrder)){ ssum$Sample <- factor(ssum$Sample, levels=sampleOrder) }
          ssum[is.na(ssum)] <- 0
          if(is.null(sumYmax)){
              sumYmax <- max(sum(ssum$CellTypeSampleDensity)) * 1.25
          }
          if(length(which(!is.na(ssum$CellTypeSampleDensity))) == 0){
              flog.info(paste0("Could not plot summary for sample ",s,", cell population ",pop))
          } else {
              flog.info("      plotting total FOV density summary")
              if(absoluteDensity){
                  plotList[["summary"]][["absolute"]] <- plotTotalFOVMarkerDensity(ssum, densityMarkers,
                                                                    ctClrs, plotTitle=plotTitle, yMax=sumYmax,
                                                                    separateLegend=FALSE, printLegend=TRUE, groupBy="Patient",
                                                                    yCol="CellTypeSampleDensity", xCol="Sample",
                                                                    cellTypeLabels=ctLabels, sampleOrder=sampleOrder,
                                                                    annotationCols=annotationCols, xAxisTitle="")
              }
              if(densityPercentage & length(unique(ssum$CellType)) > 1){
                  ctToStack <- stackable(unique(den$CellType))
                  if(!is.null(ctToStack) && length(ctToStack) > 1){
                      ssum <- ssum %>% filter(CellType %in% ctToStack)
                      flog.info("      plotting total FOV summary percentages")
                      ssumPct <- ssum %>% group_by(Sample) %>% mutate(Pct=SampleCounts/sum(SampleCounts)) 
                      ssumPct$Pct[is.na(ssumPct$Pct)] <- 0
                      plotList[["summary"]][["percentage"]] <- plotTotalFOVMarkerDensity(ssumPct, densityMarkers,
                                                                          ctClrs, sampleOrder=sampleOrder, xCol="Sample",
                                                                          plotTitle=plotTitle, yMax=sumYmax, separateLegend=FALSE,
                                                                          printLegend=TRUE, yCol="Pct", groupBy="Patient",
                                                                          yAxisTitle="Percent Total Sample Density",
                                                                          cellTypeLabels=ctLabels, pct=TRUE) 
                  }
              }
          }
         }, error = function(){
               warning(paste0("Could not plot summary for sample ",s,", cell population ",pop))
        })
    } # end if summarize
    return(plotList)
}


#' Create and store all infiltration plots in a list to be printed later
#'
#' Create and store all infiltration plots in a list to be printed later
#'
#' @param den                tibble with columns for Sample,SPOT,Band,Counts,Density
#' @param config             parsed marker config (return value from getMarkerConfig())
#' @param indivPopulations   plot populations individually; default=TRUE; if set to FALSE, ONLY stacked plots
#'                           will be created
#' @param absoluteDensity    logical indicating whether to plot density values on absolute scale (default=TRUE)
#' @param densityPercentage  logical indicating whether to plot density percentages (default=TRUE)
#' @param summarize          logical indicating whether to generate summary plots of total band densities (ALL FOVs; default=TRUE)
#' @param stacked            logical indicating whether to generate plots of all cell types for a single marker set
#'                           in a stacked bar chart (default=TRUE)
#' @param forceStack         force stacking of cell types even if they appear to be unstackable
#' @param sampleOrder        vector of sample names in order they should appear on plots
#' @param fovAreas           area data for each sample/spot
#' @param yMax               maximum y value for all plots
#' @param annotationCols     vector of column names to use for annotation
#' @return list of all plots
#' @export
getAllTotalFOVPlots <- function(den, config, indivPopulations=TRUE, absoluteDensity=TRUE,
                                densityPercentage=TRUE, summarize=TRUE, stacked=TRUE,
                                fovAreas=NULL, forceStack=FALSE, yMax=NULL, annotationCols=NULL){

    allPlots <- list()

    ## get sample/band level density info
    den <- getAllFOVDensityValues(den, fovAreas)

    if(!"Total" %in% names(den) && "Density" %in% names(den)){
        den <- den %>% select(Total=Density, everything())
    }
    
    for(ms in unique(config$MarkerSet)){
        flog.info(paste0("MARKER SET: ",ms))

        msCfg <- config %>% filter(MarkerSet == ms)
        densityMarkers <- unique(msCfg$CellType)
        msDen <- den %>% filter(CellType %in% densityMarkers)
        msYmax <- msDen %>% 
                  select(Sample,CellType,CellTypeSampleDensity) %>% 
                  unique() %>% 
                  group_by(Sample) %>% 
                  summarise(TMP=sum(CellTypeSampleDensity)) %>% 
                  pull(TMP) %>% 
                  max()

        if(indivPopulations){
            for(pop in unique(msCfg$Population)){
                flog.info(paste0("  POPULATION: ",pop))

                popCfg <- msCfg %>% filter(Population == pop)
                ## extract density for this population
                densityMarkers <- unique(popCfg$CellType)
                popDen <- msDen %>% filter(CellType %in% densityMarkers)
                popYmax <- popDen %>% 
                           select(Sample,CellType,CellTypeSampleDensity) %>% 
                           unique() %>% 
                           group_by(Sample) %>% 
                           summarise(TMP=sum(CellTypeSampleDensity)) %>% 
                           pull(TMP) %>% 
                           max()

                for(s in unique(popDen$Sample)){
                    flog.info(paste0("    SAMPLE: ",s))
                    sDen <- popDen %>% filter(Sample == s) %>% ungroup()
                    plotTitle <- paste0(s, " All FOVs:\n",pop)
                    summaryTitle <- paste0(s, " All FOVs (total density):\n",pop)
                    dm <- intersect(densityMarkers, unique(sDen$CellType))
                    allPlots[[ms]][[pop]][[s]] <- makeSampleTotalFOVDensityPlots(sDen, popCfg, dm, byFOV=TRUE,
                                                                              sampleOrder=sampleOrder, plotTitle=plotTitle, 
                                                                              yMax=yMax, sumYmax=popYmax,
                                                                              absoluteDensity=absoluteDensity,
                                                                              densityPercentage=densityPercentage,
                                                                              summarize=FALSE, summaryTitle=summaryTitle,
                                                                              forceStack=FALSE, annotationCols=annotationCols)
                } #end each sample

                if(summarize){
                    #### complete summary (all samples one plot)
                    if(!is.null(annotationCols)){
                        sDen <- popDen %>% 
                                ungroup() %>% 
                                select_at(c("Sample","CellType","CellTypeSampleDensity","SampleCounts","Patient",annotationCols)) %>% 
                                unique()
                    } else {
                        sDen <- popDen %>% ungroup() %>% select(Sample,CellType,CellTypeSampleDensity,SampleCounts,Patient) %>% unique()
                    }
                    if(!is.null(sampleOrder)){
                        sDen$Sample <- factor(sDen$Sample, levels=sampleOrder)
                    }
                    sumMax <- sDen %>% group_by(Sample) %>% summarise(TMP=sum(CellTypeSampleDensity)) %>% pull(TMP) %>% max()
                    summaryTitle <- pop
                    allPlots[[ms]][[pop]][["full_summary"]] <- makeSampleTotalFOVDensityPlots(sDen, popCfg, densityMarkers,
                                                                                  sampleOrder=sampleOrder, plotTitle=plotTitle, 
                                                                                  yMax=yMax, sumYmax=sumMax, byFOV=FALSE,
                                                                                  absoluteDensity=absoluteDensity,
                                                                                  densityPercentage=densityPercentage,
                                                                                  summarize=summarize, summaryTitle=summaryTitle,
                                                                                  forceStack=TRUE, annotationCols=annotationCols)

                }
            } #end each population
        } # end if indivPops

        ### if all negatives are set, making densities mutually exclusive, we can stack
        ### all cell types for a single marker set
               ### if all negatives are set, making densities mutually exclusive, we can stack
        ### all cell types for a single marker set
        if(stacked && (length(unique(msCfg$Label)) == nrow(msCfg))){
            ctToStack <- stackable(unique(msCfg$CellType))
            if(is.null(ctToStack) || length(ctToStack) == 0){ 
                if(forceStack){
                    flog.warn("Cell types do not appear to be stackable, but forcing stacking. Plot MAY not be accurate.")
                    densityMarkers <- unique(msCfg$CellType)
                } else {
                    flog.warn("Cell types are not stackable. No stacked plot made")
                    next
                }
            } else {
                ## get stacked version
                densityMarkers <- ctToStack
            }
        
            tDen <- msDen %>% filter(CellType %in% densityMarkers) %>% ungroup()

            for(s in unique(tDen$Sample)){
                flog.debug(paste0("  ",s))
                sDen <- tDen %>% filter(Sample == s)
                plotTitle <- paste0(s, " All FOVs: ",unique(msCfg$MarkerSetAlias))
                summaryTitle <- paste0(s, " All FOVs (total density): ",unique(msCfg$MarkerSetAlias))

                allPlots[[ms]][["stacked"]][[s]] <- makeSampleTotalFOVDensityPlots(sDen, msCfg, densityMarkers,
                                                                              sampleOrder=sampleOrder, byFOV=FALSE,
                                                                              plotTitle=plotTitle, yMax=yMax, sumYmax=msYmax,
                                                                              absoluteDensity=absoluteDensity,
                                                                              densityPercentage=densityPercentage,
                                                                              summarize=FALSE, summaryTitle=summaryTitle,
                                                                              forceStack=forceStack, annotationCols=annotationCols)
            } # end for each sample

            if(!is.null(annotationCols)){
                tDen <- tDen %>% 
                        select_at(c("Sample","CellType","CellTypeSampleDensity","SampleCounts","Patient",annotationCols)) %>% 
                        unique()
            } else {
                tDen <- tDen %>%
                        select(Sample,CellType,CellTypeSampleDensity,SampleCounts,Patient) %>%
                        unique()
            }
            if(!is.null(sampleOrder)){
                tDen$Sample <- factor(tDen$Sample, levels=sampleOrder)
            }
            plotTitle <- paste0("All FOV")
            summaryTitle <- "All FOV"
            allPlots[[ms]][["stacked"]][["full_summary"]] <- makeSampleTotalFOVDensityPlots(tDen, msCfg, densityMarkers,
                                                                              sampleOrder=sampleOrder, byFOV=FALSE,
                                                                              plotTitle=plotTitle, yMax=yMax, sumYmax=msYmax,
                                                                              absoluteDensity=absoluteDensity,
                                                                              densityPercentage=densityPercentage,
                                                                              summarize=TRUE, summaryTitle=summaryTitle,
                                                                              forceStack=forceStack, annotationCols=annotationCols)
        } # end if stacked
    } # end for each marker set
    return(allPlots)

}





#' Print plots of total FOV marker density plots 
#'
#' Print plots of total FOV marker density plots
#'
#' @param den                tibble containing densities for all markers to be plotted
#' @param fovAreas           tibble containing area data for all Samples and FOVs
#' @param config             tibble containing marker configuration (returned from getMarkerConfig())
#' @param yScaleConsistency  level on which y-scales should be the same; choices: ["byPopulation"|"all"|"markerSet"]
#' @param indivPopulations   logical; print a plot for each individual population; default=TRUE
#' @param absoluteDensity    logical; print plots of absolute density values; default=TRUE
#' @param densityPercentage  logical; print plots of density percentage values; default=TRUE
#' @param summarize          logical; print summary plots, showing density of ALL FOV together; default=TRUE
#' @param stacked            logical; for plots with all exclusive markers, print stacked bar plots; default=TRUE
#' @param forceStack         force stacking of cell types even if they appear to be unstackable
#' @param sampleOrder        vector of samples in order that they should appear on plot
#' @param annotationCols     vector of column names to use for annotation
#' @param outDir             directory where all plots will be printed, organized by marker set
#' @return nothing 
#' @export
printTotalFOVDensityPlots <- function(den, fovAreas, config, yScaleConsistency="population",indivPopulations=TRUE, 
                                      absoluteDensity=TRUE, densityPercentage=TRUE, summarize=TRUE, 
                                      stacked=TRUE, forceStack=FALSE, sampleOrder=NULL, outDir=NULL, 
                                      annotationCols=NULL){

    if(is.null(outDir)){
        outDir <- getwd()
    }

    flog.info("Getting all plots")
    allPlots <- getAllTotalFOVPlots(den, config, indivPopulations=indivPopulations,
                                    absoluteDensity=absoluteDensity, densityPercentage=densityPercentage,
                                    summarize=summarize, stacked=stacked,
                                    fovAreas=fovAreas, forceStack=forceStack, annotationCols=annotationCols)

    if(is.null(allPlots) || length(allPlots) == 0){
        flog.info("No plots made.")
        return()
    }

    for(ms in names(allPlots)){

        if(is.null(allPlots[[ms]]) || length(allPlots[[ms]]) == 0){
            flog.info(paste0("No plots made for marker set ",ms))
            next
        }
        ## create directories 
        msOutDir <- file.path(outDir,ms)
        dir.create(msOutDir, showWarnings=TRUE, recursive=TRUE)

        fovOutDir <- file.path(msOutDir, "individual_populations_individual_fovs")
        ctSumOutDir <- file.path(msOutDir, "individual_populations_whole_samples")
        stackedOutDir <- file.path(msOutDir, "stacked_populations_individual_fovs")
        stackedSumOutDir <- file.path(msOutDir, "stacked_population_whole_samples")
        sampleOutDir <- file.path(msOutDir, "individual_populations_sample_summary")
        sampleSumOutDir <- file.path(msOutDir, "stacked_populations_sample_summary")

        flog.info(paste0("Printing plots for marker set: ",ms))

        for(ct in names(allPlots[[ms]])){
            flog.info(paste0("  CellType/Population: ",ct))

            if(length(allPlots[[ms]][[ct]]) == 0){
                flog.info(paste0("No plots made for marker set ",ms," cell type ",ct))
                next
            }


            ## print sample plots separately (single page each)
            for(s in names(allPlots[[ms]][[ct]])){
                if(length(allPlots[[ms]][[ct]][[s]]) == 0){
                    flog.info(paste0("No plots made for marker set ",ms," cell type ",ct," sample ",s))
                    next
                }
                flog.info(paste0("    Sample: ",s))

                p1 <- allPlots[[ms]][[ct]][[s]][[1]][[1]]
                lgnd <- g_legend(p1) ## automatically get legend from first plot

                for(plotType in names(allPlots[[ms]][[ct]][[s]])){
                    if(ct == "stacked"){
                        if(plotType == "sep_fov"){
                            if(!file.exists(stackedOutDir)){ dir.create(stackedOutDir, showWarnings=FALSE, recursive=TRUE) }
                            fileName <- file.path(stackedOutDir, paste0(ms,"_",s,"_by_fov.pdf"))
                        } else {
                            if(!file.exists(stackedSumOutDir)){ dir.create(stackedSumOutDir, showWarnings=FALSE, recursive=TRUE) }
                            fileName <- file.path(stackedSumOutDir, paste0(ms,"_",s,".pdf"))
                        }
                   } else {
                        if(plotType == "sep_fov"){
                            if(!file.exists(fovOutDir)){ dir.create(fovOutDir, showWarnings=FALSE, recursive=TRUE) }
                            fileName <- file.path(fovOutDir, paste0(ct,"_",s,"_by_fov.pdf"))
                        } else {
                            if(!file.exists(ctSumOutDir)){ dir.create(ctSumOutDir, showWarnings=FALSE, recursive=TRUE) }
                            fileName <- file.path(ctSumOutDir, paste0(ct,"_",s,".pdf"))
                        }
                    }
                    plts <- allPlots[[ms]][[ct]][[s]][[plotType]]
                    if(is.null(lgnd) || is.null(plts)){ next }
                    pdf(fileName, height=8.5, width=14, onefile=TRUE)
                    flog.debug("drawing legend")
                    if(!is.null(lgnd)) { grid.arrange(lgnd) }
                    for(p in names(plts)){
                        #plot.new()
                        p1 <- plts[[p]] 
                        if(!is.null(p1)){ 
                            p1$grobs[[which(sapply(p1$grobs, function(x) x$name) == "guide-box")]] <- zeroGrob()
                            grid.arrange(p1) 
                        }
                    }
                    dev.off()
                }
            }   
        }   
    }
}

#' Make annotated heatmap of cell type densities
#' 
#' Make annotated heatmap of cell type densities
#' 
#' @param den              table containing total FOV cell type densities with basic FOV/Sample annotations
#' @param annot            table of additional FOV annotation to include in annot tracks
#' @param annotColors      list of annotation legend colors, where keys match column names in annot, and
#'                         each value is a vector of colors named by the respective annotation values 
#' @param msToSummarize    a vector of marker set names to collapse into one "cell type"; if given, must 
#'                         provide a parsed cell type configuration and FOV area table; any mutually-exclusive 
#'                         (stackable) cell types belonging to a marker set in this vector will be converted to 
#'                         a single density value for each FOV.
#' @param ctConfig         if summarizing any marker sets, must provide a cell type config returned by
#'                         getMarkerConfig()
#' @param areas            if summarizing any marker sets, must provide table of FOV areas in order to calculate
#'                         new density. 
#' @export
makeDensityHeatmap <- function(den, annot=NULL, annotColors=NULL, msToSummarize=NULL, ctConfig=NULL, 
                               areas=NULL, separateLegend=FALSE, bySample=TRUE, byFOV=TRUE, clusterFOVs=TRUE,
                               clusterSamples=TRUE){

    den <- den %>% 
           filter(!is.na(Counts)) %>%
           select(-(SampleLabel)) %>%
           mutate(SampleLabel = paste0(Patient,".",Sample_number,".",SPOT)) %>%
           filter(CellType %in% ctConfig$CellType) %>%
           unique() 

    den <- getAllFOVDensityValues(den, areas)

    ## collapse any marker sets marked for summarization
    if(!is.null(msToSummarize) && length(msToSummarize) > 0){
        areas$SampleLabel <- paste0(areas$Patient,".",areas$Sample_number,".",fovAreas$SPOT)
        areas <- areas %>% select(SampleLabel, Area)
        areas$Sample <- areas$SampleLabel
        for(ms in msToSummarize){
            msCfg <- ctConfig %>% filter(MarkerSet == ms)
            msDen <- den %>% filter(CellType %in% stackable(unique(msCfg$CellType)))
            msCounts <- msDen %>% group_by(Sample,SPOT) %>% summarize(msCounts=sum(Counts), msSampCounts=sum(SampleCounts))
            msDen <- msDen %>% 
                     left_join(msCounts, by=c("Sample","SPOT")) %>%
                     mutate(msDensity=msCounts/Area, msSampDensity=msSampCounts/SampleArea) %>%
                     select(-c(CellType,Counts,Total,SampleCounts,CellTypeSampleDensity)) %>%
                     select(Counts=msCounts, Total=msDensity, SampleCounts=msSampCounts, CellTypeSampleDensity=msSampDensity, everything()) %>%
                     mutate(CellType=ms) %>%
                     select(names(den)) %>% 
                     unique()

            den <- den %>% 
                   filter(!CellType %in% stackable(unique(msCfg$CellType))) %>%
                   bind_rows(msDen) 
        }
    }

    den$Total <- log2(den$Total+1)
    den <- den %>% filter(!is.na(Total)) 

    if("CellTypeSampleDensity" %in% names(den)){
        den$CellTypeSampleDensity <- log2(den$CellTypeSampleDensity + 1)        
    }

    for(statLevel in c("byFOV","bySample")){

        if(!get(statLevel)){ next }

        if(statLevel == "bySample"){
            den2 <- den %>% 
                    mutate(SampleLabel=paste0(Patient,".",Sample_number)) %>%
                    filter(!is.na(CellTypeSampleDensity)) %>% 
                    select(-c(Counts,Total,SPOT,FOV_type,Sample_number,Sample,Area)) %>% 
                    select(Total=CellTypeSampleDensity, Counts=SampleCounts, everything()) %>% 
                    unique()
            clusterCols <- clusterSamples
        } else {
            den2 <- den
            clusterCols <- clusterFOVs
        }

        ## format data into numeric matrix for plotting
        htmpDen <- den2 %>% 
                   select(SampleLabel, CellType, Total) %>% 
                   unique()
        dups <- which(duplicated(htmpDen %>% select(SampleLabel,CellType)))
        if(length(dups) > 0){
            htmpDen <- htmpDen[-dups,] 
        }
        htmpDen <- htmpDen %>%
                   filter(!is.na(Total)) %>%
                   spread(CellType,Total) 

        tmp <- as.matrix(htmpDen)
        rownames(tmp) <- htmpDen$SampleLabel
        tmp <- tmp[,-1]
        tmp2 <- apply(tmp, 2, as.numeric)
        rownames(tmp2) <- rownames(tmp)
        idxs <- grep("-", colnames(tmp2))
        colnames(tmp2)[idxs] <- unlist(lapply(strsplit(colnames(tmp2)[idxs],","), function(x){ paste(x[-grep("-",x)],collapse=",") }))
        tmp2 <- t(tmp2) 

        ## prep annotation tracks
        if(!is.null(annot)){
            annot2 <- annot 
            if(!is.null(annotColors)){
                selVars <- unique(c(names(annotColors),"Patient","Sample_number","SampleLabel"))
                if(statLevel == "byFOV"){
                    selVars <- c(selVars, "SPOT")
                    annot2 <- annot2 %>% 
                              mutate(SampleLabel=paste0(Patient,".",Sample_number,".",SPOT))
                } else {
                    annot2 <- annot2 %>%
                              mutate(SampleLabel=paste0(Patient,".",Sample_number))
                }
            }
            annot2 <- annot2 %>% select(selVars) %>% unique()
            annot2 <- annot2 %>% 
                     filter(SampleLabel %in% colnames(tmp2)) %>%
                     unique()
            annot2[is.na(annot2)] <- "None"
            for(x in names(annot2)){
                lvls <- unique(annot2[[x]])
                addNA <- "NA" %in% lvls 
                if(addNA){
                    lvls <- lvls[-which(lvls == "NA")]
                }
                if(!any(is.na(as.numeric(lvls)))){
                    lvls <- as.numeric(lvls)
                } 
                lvls <- sort(lvls)
                if(addNA){ lvls <- c(lvls, "NA") }                    
                annot2[[x]] <- factor(annot2[[x]], levels=lvls)
            }

            rn <- as.vector(annot2$SampleLabel)
            annot2 <- as.data.frame(annot2)
            rownames(annot2) <- rn
            annot2 <- annot2[,-which(colnames(annot2) %in% c("SampleLabel","Sample_number","SPOT"))]
            annot2 <- annot2 %>% select(rev(names(annot2)))
        }

        breaks <- (0:ceiling(max(den2$Total)))
        col2 <- colorRampPalette(c("white", "navy"))(length(breaks))

        if(separateLegend){
            tmp4 <- matrix(0, ncol=ncol(tmp2), nrow=nrow(tmp2))
            rownames(tmp4) <- rownames(tmp2)
            colnames(tmp4) <- colnames(tmp2)
            pheatmap(mat = tmp4, 
                     annotation_col=annot2, 
                     fontsize=8, color=col2, 
                     breaks=breaks, 
                     show_colnames=FALSE, 
                     show_rownames=FALSE, 
                     cellwidth=.1)
        }  

        binned <- apply(tmp2, 2, function(x){ cut(x, breaks, include.lowest=TRUE, labels=breaks[-length(breaks)])  }) 
        rownames(binned) <- rownames(tmp2)
        binned <- apply(binned, 2, as.numeric)
        rownames(binned) <- rownames(tmp2)    

        cellwidth <- ifelse(statLevel == "bySample", 600/dim(binned)[2], 2)
        cellheight <- ifelse(statLevel == "bySample", cellwidth, 600/dim(binned)[1])
        pheatmap(mat = binned, 
                 color = col2,
                 main = paste0("Total FOV cell type density (log2)\n", tolower(gsub("by","by ",statLevel))),
                 legend = !separateLegend, 
                 legend_breaks = breaks[-length(breaks)],
                 #legend_labels = llabels, 
                 drop_levels=TRUE,
                 annotation_legend = !separateLegend, 
                 annotation_col = annot2, 
                 annotation_colors = annotColors, 
                 fontsize = 8, 
                 breaks = breaks, 
                 show_colnames = FALSE, 
                 cellwidth = cellwidth,
                 cellheight = cellheight, 
                 scale = "none",
                 cluster_cols=clusterCols)
    }
    return()
}
