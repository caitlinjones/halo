#' Get marker colors
#'
#' Get a named vector of predefined colors for certain markers
#'
#' @param markerColors  named vector of colors already in use, where names are 
#'                      marker names (e.g., CD3, CD4, MRC1, etc.) and values are
#'                      colors; only to be used when 'custom' is set to FALSE; default=NULL
#' @param markerCombos  vector of marker combinations to add to MarkerColors; only to be
#'                      used when 'custom' is set to FALSE
#' @export
getMarkerColors <- function(markerColors=NULL,markerCombos=NULL){
    qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual' | brewer.pal.info$category == 'div',]
    col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    col_vector <- col_vector[-which(col_vector == "#FFFF99")] ## remove ugly yellow

    if(is.null(markerColors) || length(markerColors) == 0){
        markerColors <- col_vector[1:length(markerCombos)]
        names(markerColors) <- markerCombos
    } else {
        newCombos <- setdiff(markerCombos,names(markerColors))
        numExtras = length(setdiff(markerCombos,names(markerColors)))
        if(numExtras > 0){
            first = length(markerColors)+1
            last = length(markerColors)+numExtras
            extraColors <- col_vector[first:last]
            names(extraColors) <- newCombos
            markerColors <- c(markerColors,extraColors)
        }
    }
    return(markerColors)
}

#' Get marker colors
#'
#' Get a named vector of predefined colors for certain markers
#'
#' @param markerColors  named vector of colors already in use, where names are 
#'                      marker names (e.g., CD3, CD4, MRC1, etc.) and values are
#'                      colors; only to be used when 'custom' is set to FALSE; default=NULL
#' @param markerCombos  vector of marker combinations to add to MarkerColors; only to be
#'                      used when 'custom' is set to FALSE
#' @export
getCustomColors <- function(markerColors=NULL,markerCombos=NULL){
    #CD3, CD4, CD8, FOXP3, AND/OR CD56 = different shades of blue
    #CD20 = green
    #SOX10 = yellow
    #CD4, CD68, CD14, CD163, MRC1, AND/OR TGM2 = different shades of red
    #any other combination of these markers = different shades of grey
    #any combination less than 2%/"other" = white

    bls <- colorRampPalette(c("lightblue", "darkblue"))(length(grep("CD3|CD4|CD8|FOXP3|CD56",markerCombos))*10)
    bls <- setdiff(bls,markerColors)
    reds <- colorRampPalette(c("red", "darkred"))(length(grep("CD4|CD68|CD14|CD163|MRC1|TGM2",markerCombos))*10)
    reds <- setdiff(reds,markerColors)
    grays <- colorRampPalette(c("white","black"))(256)   
    grays <- setdiff(grays,markerColors)
    #purples <- colorRampPalette(c("purple","black"))(246)
    browns <- colorRampPalette(c("brown","black"))(175)
    darkblues <- colorRampPalette(c("darkblue","black"))(140)
    darkgreens <- colorRampPalette(c("darkgreen","black"))(101)
    lightyellow <- colorRampPalette(c("lightyellow","white"))(32)
    pinks <- colorRampPalette(c("white","pink"))(103)
    other <- c(grays,browns,darkgreens,lightyellow)#,pinks)

    if(length(markerCombos) < 100){
        extras <- grays
    } else {
        extras <- other
    }

    for(m in markerCombos){
        mlst <- unlist(strsplit(m,","))
        if(m == "CD4"){
            newColor = "#094dba"
        } else if(m == "SOX10"){
            newColor <- "#ffff80"
        } else if(m == "CD20"){
            newColor <- "#28b463"
        } else if(all(mlst %in% c("CD3","CD4","CD8","FOXP3","CD56"))){
            n <- sample(1:length(bls),1)
            newColor <- bls[n]
            bls <- bls[-n]
        } else if(all(mlst %in% c("CD4","CD68","CD14","CD163","MRC1","TGM2"))){
            n <- sample(1:length(reds),1)
            newColor <- reds[n] 
            reds <- reds[-n]
        } else {
            n <- sample(1:length(extras),1)
            newColor <- extras[n]
            extras <- extras[-n]
        }
        markerColors <- c(markerColors,newColor)
        names(markerColors)[length(markerColors)] <- m
    }

    return(markerColors)
}

#' Extract a legend from a gplot
#' 
#' Extract a legend from a plot generated with ggplot
#' 
#' @param plt    a plot generated with ggplot 
getLegend <- function(plt){ 
    tmp <- ggplot_gtable(ggplot_build(plt)) 
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

#' Make a grid of pie charts, one for each given marker showing
#' the percentages of all the marker combinations in which that 
#' positive marker is included
#' 
#' @param   countsTable         tibble of counts generated by countMarkers()
#' @param   markerNames         vector of positive markers, for each of which a pie will be made
#' @param   type                type of chart ["bar"|"pie"]; default="bar"
#' @param   other_threshold     all markers making up less that this percentage will be put togther in 
#'                              a marker combination named 'other'; default=0.05
#' @param   custom_colors       use color scheme defined in getMarkerColors(); TO DO: CHANGE THIS
#' @param   exclude_sample_fov  a string of Sample FOVs to exclude from analysis in the form: 
#'                                  'Sample1:2+4+3,Sample2:1,Sample3:2+4' 
#' @param   pdfFile             output PDF file
#' @export
plotMarkerPercentages <- function(countsTable, markerNames, type="bar", other_threshold=0.05, 
                                    exclude_sample_fov=NULL, pdfFile=NULL, custom_colors=FALSE){ 

    ## set up colors
    markerColors <- c()

    if(!is.null(pdfFile)){
        pdf(pdfFile,width=10,height=8)#paper="a4r")
    }

    ## remove median rows from counts table
    countsTable <- filter(countsTable, !is.na(FOV) & !is.na(SLICE))

    ## remove any exclusions
    if(!is.null(exclude_sample_fov)){
        countsTable <- removeExclusions(countsTable, exclude_sample_fov)
    }

    allPlots <- c()
    for(samp in unique(countsTable$Sample)){
        flog.debug(samp)
        pieTbl <- tibble()
        for(m in markerNames){
            flog.debug(m)
            if(length(grep(m,names(countsTable))) == 0){
                flog.debug("skipping %s",m)
                next
            }

            tmp <- countsTable
            ## remove marker combos in counts table that do not have this POSITIVE marker
            posMarkerCombos <- names(tmp)[!grepl(paste0(m,"-"),names(tmp))]
            if(length(posMarkerCombos) == 3){
                flog.debug("No %s cells. Skipping.",m)
                next
            }
            ## get total marker counts and percentage of each marker for this sample
            tmp <- dplyr::select(tmp, posMarkerCombos) %>%
                 filter(Sample == samp) %>%
                 gather(-(Sample:SLICE), key="Marker", value="Count") %>%
                 group_by(Marker) %>%
                 summarise(TotalMarkerCount = sum(Count)) %>%
                 mutate(PercentAllSampleCells = TotalMarkerCount/sum(TotalMarkerCount))
 
            ## collapse markers that each make up less than X% of all counts
            other <- filter(tmp, PercentAllSampleCells <= other_threshold) %>%
                       summarise(TotalMarkerCount = sum(TotalMarkerCount),
                                 PercentAllSampleCells = sum(PercentAllSampleCells))
            other$Marker <- paste0("Other = <",other_threshold*100,"%")
            ## separate markers that make up more than X% of all counts
            all <- filter(tmp, PercentAllSampleCells > other_threshold) %>%
                     bind_rows(other) %>%
                     arrange(desc(PercentAllSampleCells)) %>%
                     dplyr::select(-PercentAllSampleCells)  %>%
                     spread(Marker,TotalMarkerCount)
            if(ncol(all) > 1){
                all$Marker <- paste0(m," : ",rowSums(all[1,2:ncol(all)]))
            } else {
                all$Marker <- paste0(m," : ",all$Other)
            } 
            all <- dplyr::select(all, Marker, everything())
            ## add data for this sample to the rest
            pieTbl <- bind_rows(pieTbl, all)
        }

        pie_data <- gather(pieTbl, names(pieTbl)[-1], key="Combo", value="Count") %>%
                    group_by(Marker) %>%
                    mutate(labels=paste0(round((Count/sum(Count,na.rm=TRUE)) * 100, 1),"%")) %>%
                    drop_na()
        pie_data <- filter(pie_data, Count > 0)
        ## remove negative markers from the combos for clarity 
        pie_data$Combo <- gsub(",[A-Z]{1,}\\d*-","",pie_data$Combo)

        flog.info("Sample %s has %s marker combinations",samp,length(unique(pie_data$Combo)))

        ## get colors for this sample to match any previous samples
        if(custom_colors){
            markerColors <- getCustomColors(markerColors=markerColors,markerCombos=unique(pie_data$Combo))
        } else {
            markerColors <- getMarkerColors(markerColors=markerColors,markerCombos=unique(pie_data$Combo))
        }
        sampleColors <- markerColors[unique(pie_data$Combo)] 

        ## plot data
        numRows <- NULL
        if(type == "bar"){ numRows <- 1 }

        ## TEMPORARY HACK until I can figure out how to print legend separately: 
        ## print pies without legend first in order to make sure they are bigger 
        ## and readable; then add legend if <= 70 combos 
        p <- ggplot(pie_data, aes(x = factor(1), y = Count, fill=Combo)) +
                   facet_wrap(~Marker,nrow = numRows) +
                   geom_bar(width=1,stat="identity",position="fill") +
                   labs(x = "", y = "", title = paste0(samp,"\n")) +
                   scale_fill_manual(values = sampleColors) +
                   theme(axis.text.x = element_blank(),
                   legend.position = "none",
                     strip.text.x = element_text(size=6)) +
                   guides(fill=guide_legend(title="",nrow=length(unique(pie_data$Combo))/4))
        if(length(unique(pie_data$Combo)) <= 70){
            p <- p + theme(legend.text = element_text(size=6),
                       legend.key.size = unit(0.3, "cm"),
                       legend.position = "right")
        }
        if(type == "pie"){
            p <- p + coord_polar("y")
        }

        print(p)

    }

    if(!is.null(pdfFile)){
        dev.off()
    }
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
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

#' Get plot theme for a given type of halo plot
#' 
#' Depending on what type of plot you're making, get a theme()
#' object with the elements standard to that plot type
#'
#' @param plotType  character string indicating which type of plot youre making;
#'                  options: ["infiltration density"]
#' @return theme() object with standard theme elements
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
                  plot.title = element_text(size=20, margin=margin(b=0.5, t=0.2, r=0, l=0.1, "cm"))
                  #panel.background = element_blank(),
                  #panel.grid.major = element_blank(),
                  #panel.grid.minor = element_blank()
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
              err = function(){ mDen$CellTypeLabels <- mDen$CellType
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
         geom_bar(stat="identity", position="stack", color="black", size=0.05, width=0.75) +
         plotTheme +
         scale_fill_manual(values=clrs, labels=names(clrs)) +
         scale_x_discrete(breaks=c(-355,-305,-205,-105,-5,95,195,295,355),
                          labels=c("-360","-300","-200","-100","0","100","200","300","360"),
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

    if(separateLegend | legendOnly){
        flog.debug("getting legend")
        legend <- g_legend(p1)
        if(printLegend & length(unique(mtd$CellType)) > 1){
            flog.debug("drawing legend")
            grid.draw(legend)
        }
        p1 <- p1 + theme(legend.position="none")
    }
    if(!legendOnly){
        flog.debug("printing plot")
        print(p1)
    }
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


#' Get a named vector of colors to use on total density bar chart
#' 
#' Using a tibble including at least columns for Sample,SPOT,Functional, 
#' create a named vector to be added to that tibble and used to fill
#' bars on a chart for total density
#'
#' @param dat              tibble including at least columns for Sample,SPOT,Functional
#' @param funcPosColor       color to represent density of cells positive for a functional marker
#' @param funcNegColor       color to represent density of cells negative for a functional marker
#' @param sampleColors     vector of colors, each named according to the sample it represents
#' @param sampleColorsFuncPos vector of colors for density of cells positive for functional marker, 
#'                            named according to the samples they represent
#' @return vector of colors, also named by those colors
getTotalDensityFillColors <- function(dat, funcPosColor, funcNegColor, sampleColors, sampleColorsFuncPos){
    ### color all bars grey (darker for funcMarker+)
    ### color median bars according to sample color (again, darker for funcMarker+)
    clrs <- rep(funcPosColor,len(dat$Sample))
    clrs[grep("\\-",dat$Functional)] <- funcNegColor
    for(s in unique(dat$Sample)){
        clrs[intersect(which(dat$Sample == s & dat$SPOT == "Median"),
                       grep("\\-",dat$Functional))] <- sampleColors[[s]]
        clrs[intersect(which(dat$Sample == s & dat$SPOT == "Median"),
                       grep("\\+",dat$Functional))] <- sampleColorsFuncPos[[s]]
    }

    names(clrs) <- clrs

    return(clrs)
}

#' Build bar plot for total cell type density
#'
#' Build bar plot for total cell type density
#'
#' @param  dat                 tibble containing all info needed for plotting [ TO DO: FILL IN EXACTLY WHAT'S NEEDED ]
#' @param  plotPer             generate plot per ['markerSet'|'marker'|'sample'|'sample+marker]; default='markerSet' 
#' @param  funcMarker          marker to be indicated as either positive or negative for each cell type
#' @param  funcNegColor        hex color indicating functional marker is negative in that cell type
#' @param  funcPosColor        hex color indicating functional marker is positive in that cell type
#' @param  sampleColor         vector of colors named by the sample which it will represent; this color will only
#'                             appear on the median bars
#' @param  sampleColorDark     vector of colors named by sample, to represent density of cells that are positive
#'                             for the functional marker, if one is given
#' @param  consistentYScaling  logical; when onePanel=FALSE, set this parameter to FALSE in order to scale separate 
#'                             plots individually; default=TRUE
#' @return list of ggplots to be printed
getTotalDensityPlots <- function(dat, markerList, plotPer='markerSet', funcMarker=NULL, funcPosColor="grey35",
                             funcNegColor="grey50", sampleColor="#f16913", sampleColorDark="#80380A",
                             consistentYScaling=TRUE){


    ## set up theme
    plotTheme = theme(legend.title=element_blank(),
                     axis.text=element_text(size=8),
                     axis.text.x = element_text(angle=45),
                     axis.ticks = element_blank(),
                     strip.text.x = element_text(face="bold", size=12, margin=margin(.4, .2, .4, .2, "cm")),
                     strip.text.y = element_text(margin = margin(.4, .2, .4, .2, "cm")),
                     strip.placement = "outside",
                     plot.background = element_blank(),
                     panel.background = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.line = element_line(color = "black", size=0.25),
                     legend.position = "right"
                     )

    ymax <- ifelse(consistentYScaling, max(dat$Density), NULL)

    #####
    ## generate plot(s)
    #####
    plotList <- list()

    dplot <- function(dat,funcMarker=NULL,ymax=NULL){
      clrs <- dat$Colors
      names(clrs) <- clrs

      ## remove sample name from SampleSpot for x axis labels
      xlbls <- gsub(".*_","",dat$SampleSpot)
      names(xlbls) <- dat$SampleSpot

      ## set up temporary legend labels in case there is no functional marker
      legendLabels <- c("TMP","TMP")
      if(!is.null(funcMarker)){
        legendLabels <- paste0(funcMarker,c("-","+"))
      }
      g <- ggplot(dat, aes(x=as.factor(SampleSpot2), y=Density,
                        fill=factor(Colors,levels=c(funcNegColor,funcPosColor,sampleColor,sampleColorDark)))) +
                        geom_bar(stat="identity",position="stack") +
                        plotTheme +
                        xlab("FOV") +
                        ylab("Density (counts/mm^2)") +
                        scale_fill_manual(values=clrs,
                                          breaks=c(funcNegColor,funcPosColor),
                                          labels=legendLabels) +
                        scale_x_discrete(labels=xlbls) +
                        geom_errorbar(aes(ymin=ymin,ymax=ymax))
      if(is.null(funcMarker)){
          g <- g + theme(legend.position="none")
      }
      if(!is.null(ymax)){
         g <- g + scale_y_continuous(limits = c(0,ymax))
      }
      return(g)
    }


    if(plotPer %in% c('markerSet','marker','sample')){
        flog.debug("building one plot per %s",plotPer)
        if(plotPer == 'markerSet'){
            for(mn in names(markerList)){
print(paste0('making plot for markerSet ',mn))
                mdat <- dat %>% filter(CellType %in% markerList[[mn]])
                g <- dplot(mdat, funcMarker=funcMarker, ymax=ymax)
                g <- g + facet_grid(CellTypeLabels ~ Sample2, scales="free_x", switch="y")
                plotList[[mn]] <- g
            }
        } else if(plotPer == 'sample') {
            for(s in unique(dat$Sample)){
                mdat <- dat %>% filter(Sample==s)
print(paste0('making plot for sample ',s))
                g <- dplot(mdat, funcMarker=funcMarker, ymax=ymax)
                g <- g + facet_wrap(~CellType, scales="free_x", strip.position=c("top")) +
                     labs(title=s) +
                     theme(strip.text.x = element_text(margin = margin(.4, .2, .4, .2, "cm")))
                plotList[[s]] <- g
            }
        } else { ## plotPer == 'marker'
            for(m in unlist(markerList)){
                if(len(grep(paste0(",",funcMarker),m) > 0)){ next }
print(paste0('making plot for marker ',m))
                mdat <- dat %>% filter(CellType==m)
                g <- dplot(mdat, funcMarker=funcMarker, ymax=ymax)
                g <- g + facet_wrap(~Sample, scales="free_x", strip.position=c("top")) +
                     labs(title=m)
                plotList[[m]] <- g
            }
        }
    } else {
        for(marker in unique(dat$CellType)){
            for(samp in unique(dat$Sample)){
print(paste0('making plot for marker ',marker,', sample ',samp))
                mdat <- dat %>% filter(CellType==marker,Sample==samp)
                flog.debug("printing each marker individually")
                g <- dplot(mdat,funcMarker=funcMarker,ymax=ymax) +
                     labs(title=paste0(samp, "\n", marker))
                if(is.null(funcMarker)){
                    g <- g + theme(legend.position="none")
                }
                plotList[[samp]][[marker]] <- g
            }
        }
    }
    return(plotList)
}


#' Print density plots 
#'
#' Print density plots
#'
#' @param den                tibble containing densities for all markers to be plotted
#' @param config             configuration of plots (already parsed)
#' @param bandWidth          if plotting by distance intervals from infiltration boundary, this is the size of 
#'                           each interval
#' @param yScaleConsistency  level on which y-scales should be the same; choices: ["byPopulation"|"all"|"markerSet"]
#' @param absoluteDensity    logical; print plots of absolute density values; default=TRUE
#' @param densityPercentage  logical; print plots of density percentage values; default=TRUE
#' @param byFOV              logical; print plots of each FOV (all FOV on one panel); default=TRUE
#' @param summarize          logical; print summary plots, showing density of ALL FOV together; default=TRUE
#' @param stacked            logical; for plots with all exclusive markers, print stacked bar plots; default=TRUE
#' @return nothing 
#' @export
printDensityPlots <- function(den, bandWidth=NULL, config, yScaleConsistency="population", absoluteDensity=TRUE,
                              densityPercentage=TRUE, byFOV=TRUE, summarize=TRUE, stacked=TRUE,
                              sampleOrder=NULL, separateLegend=TRUE){

    yMax <- NULL
    if(yScaleConsistency=="all"){
        yMax <- max(den$Density)
    }
    curdir <- getwd()
    for(ms in unique(config$MarkerSet)){
        idMarkerTag <- "+/-"
        if(grepl("_with_neg",ms)){
            idMarkerTag <- "-*"
        }
        print(paste0("MARKER SET: ",ms))

        ## create directories 
        msOutDir <- file.path(curdir,ms)
        dir.create(msOutDir, showWarnings=TRUE, recursive=TRUE)
        fovOutDir <- file.path(msOutDir,"by_fov")
        dir.create(fovOutDir, showWarnings=TRUE, recursive=TRUE)
        sumOutDir <- file.path(msOutDir,"summaries")
        dir.create(sumOutDir, showWarnings=TRUE, recursive=TRUE)

        msCfg <- config %>% filter(MarkerSet == ms)
        if(yScaleConsistency == "markerSet"){
            yMax <- max(den$Density[which(den$CellType %in% msCfg$CellType)])
        }

        for(pop in unique(msCfg$Population)){
            print(paste0("  POPULATION: ",pop))

            popCfg <- msCfg %>% filter(Population == pop)

            ## set up labels and colors for this population
            ctLabels <- popCfg$Label
            names(ctLabels) <- popCfg$CellType

            clrLbls <- unique(popCfg[,c("Color","Label")])
            ctClrs <- pull(clrLbls[,"Color"])
            names(ctClrs) <- pull(clrLbls[,"Label"])

            ## extract density for this population
            densityMarkers <- unique(popCfg$CellType)
            ctDen <- den %>% filter(CellType %in% densityMarkers)

            sumYmax <- NULL
            if(yScaleConsistency == "population"){
                yMax <- max(ctDen$Density)
                if(any(ctDen$Density > 0)){
                    mssum <- suppressMessages(getInfiltrationDensitySummary(densityMarkers, ia=ia, ba=ba))
                    sumYmax <- max(mssum$TotalDensity)
                }
            }
            for(s in unique(ctDen$Sample)){
                print(paste0("    SAMPLE: ",s))
                sDen <- ctDen %>% filter(Sample == s) %>% filter(complete.cases(.))

                tryCatch({
                    onefile=TRUE
                    if(length(unique(sDen$CellType)) == 1){ onefile <- FALSE }
                    pdfFile <- file.path(fovOutDir,
                                         paste0(pop,"_",s,"_infiltration_density_by_cell_type_and_FOV.pdf"))
                    pdf(pdfFile, height=11,width=8.5, onefile=onefile)

                    plotTitle <- paste0(s, " Interface FOVs: ",pop,",",idMarkerTag)

                    if(absoluteDensity){
                        print("      plotting absolute infiltration density")
                        ## print absolute density with legend on separate page
                        plotInfiltrationDensity(sDen, densityMarkers, bandWidth, ctClrs, sampleOrder=sampleOrder,
                                                plotTitle=plotTitle, yMax=yMax, separateLegend=TRUE,
                                                facetByFOV=TRUE, facetByCellType=FALSE, cellTypeLabels=ctLabels)
                    }
                    if(densityPercentage & length(unique(sDen$CellType)) > 1){
                        ## plot density percentage
                        if(length(which(sDen$Density > 0)) > 1){
                            print("      plotting intiltration density percentages")
                            sDenPct <- suppressMessages(getDensityPercentage(sDen))
                            plotInfiltrationDensity(sDenPct, densityMarkers, bandWidth, ctClrs, sampleOrder=sampleOrder,
                                                    plotTitle=plotTitle, yMax=yMax, separateLegend=TRUE,
                                                    printLegend=FALSE, yCol="Percent", facetByFOV=TRUE,
                                                    facetByCellType=FALSE, yAxisTitle="Percent Total Density",
                                                    cellTypeLabels=ctLabels,pct=TRUE)
                        }
                    }
                }, err = function(){
                    warning(paste0("Could not plot sample ",s,", for cell population ",pop))
                }, finally={ dev.off(); })

                if(summarize){
                    tryCatch({
                        ## plot summary over all FOV
                        sia <- ia %>% filter(Sample == s)
                        sba <- ba %>% filter(Sample == s)
                        ssum <- suppressMessages(getInfiltrationDensitySummary(densityMarkers, ia=sia, ba=sba))
                        if(length(which(!is.na(ssum$TotalDensity))) == 0){
                            print(paste0("Could not plot summary for sample ",s,", cell population ",pop))
                        } else {
                            onefile=TRUE
                            if(length(unique(ssum$CellType)) == 1){ onefile <- FALSE }
                            pdfFile <- file.path(sumOutDir,paste0(pop,"_",s,"_infiltration_density_summary.pdf"))
                            pdf(pdfFile, height=8.5, width=11, onefile=onefile)
                            print("      plotting infiltration summary")
                            ## plot summary over all FOV
                            plotTitle <- paste0(s, " Interface All FOV: ",pop,",",idMarkerTag)
                            plotInfiltrationDensity(ssum, densityMarkers, bandWidth, ctClrs, plotTitle=plotTitle, yMax=sumYmax,
                                                    separateLegend=TRUE, printLegend=TRUE, facetByFOV=FALSE,
                                                    facetByCellType=FALSE, yCol="TotalDensity",cellTypeLabels=ctLabels)
                            if(densityPercentage & length(unique(ssum$CellType)) > 1){
                                print("      plotting intiltration density percentages")
                                ssumPct <- suppressMessages(getDensityPercentage(ssum, by=c("Sample","Band")))
                                plotInfiltrationDensity(ssumPct, densityMarkers, bandWidth, ctClrs, sampleOrder=sampleOrder,
                                                        plotTitle=plotTitle, yMax=yMax, separateLegend=TRUE,
                                                        printLegend=FALSE, yCol="Percent", facetByFOV=FALSE,
                                                        facetByCellType=FALSE, yAxisTitle="Percent Total Density",
                                                        cellTypeLabels=ctLabels,pct=TRUE)
                            }
                        }
                    }, err = function(){
                        warning(paste0("Could not plot summary for sample ",s,", cell population ",pop))
                    },finally={ dev.off(); }) 
                }
            }
        }

        ### if all negatives are set, making densities mutually exclusive, we can stack
        if(stacked & idMarkerTag == "-*"){
            if(length(unique(msCfg$CellType)) == nrow(msCfg)){
                ## get stacked version
                densityMarkers <- unique(msCfg$CellType)
                tDen <- den %>% filter(CellType %in% densityMarkers)
   
                ctLabels <- msCfg$Label
                names(ctLabels) <- msCfg$CellType

                clrLbls <- unique(msCfg[,c("Color","Label")])
                ctClrs <- pull(clrLbls[,"Color"])
                names(ctClrs) <- pull(clrLbls[,"Label"]) 

                ## check for unique labels; if they are not unique as-is, revert back to just using cell types;
                ## TO DO: remove negatives, but only in population, not functional part
                numLabels <- msCfg %>% select(MarkerSet, Label) %>% group_by(MarkerSet) %>% summarize(nLbls=length(Label)) %>% pull(nLbls)
                unqLabels <- msCfg %>% select(MarkerSet, Label) %>% group_by(MarkerSet) %>% summarize(nLbls=length(Label)) %>% pull(nLbls)
                if(unqLabels != numLabels){
                    ctLabels <- msCfg$CellType
                    names(ctLabels) <- msCfg$CellType
                    names(ctClrs) <- ctLabels
                }

                for(s in unique(tDen$Sample)){
                    sDen <- tDen %>% filter(Sample == s)
                    tryCatch({
                        onefile=TRUE
                        #if(length(unique(sDen$CellType)) == 1){ onefile <- FALSE }
                        pdfFile <- file.path(paste0(ms,"_",s,"_infiltration_density_by_cell_type_and_FOV.pdf"))
                        pdf(pdfFile, height=11, width=8.5, onefile=onefile)

                        plotTitle <- paste0(s, " Interface FOVs: ",unique(msCfg$MarkerSetAlias)," (Cell ID Markers = ",idMarkerTag,")")

                        if(absoluteDensity){
                            print("      plotting absolute infiltration density (STACKED)")
                            ## print absolute density with legend on separate page
                            plotInfiltrationDensity(sDen, densityMarkers, bandWidth, ctClrs, sampleOrder=sampleOrder,
                                                    plotTitle=plotTitle, yMax=yMax, separateLegend=TRUE,
                                                    facetByFOV=TRUE, facetByCellType=FALSE, cellTypeLabels=ctLabels)
                        }
                        if(densityPercentage & length(unique(sDen$CellType)) > 1){
                            ## plot density percentage
                            if(length(which(sDen$Density > 0)) > 1){
                                print("      plotting intiltration density percentages (STACKED)")
                                sDenPct <- suppressMessages(getDensityPercentage(sDen))
                                plotInfiltrationDensity(sDenPct, densityMarkers, bandWidth, ctClrs, sampleOrder=sampleOrder,
                                                        plotTitle=plotTitle, yMax=yMax, separateLegend=TRUE,
                                                        printLegend=FALSE, yCol="Percent", facetByFOV=TRUE,
                                                        facetByCellType=FALSE, yAxisTitle="Percent Total Density",
                                                        cellTypeLabels=ctLabels, pct=TRUE)
                            }
                        }
                    }, err = function(){
                        warning(paste0("Could not plot sample ",s,", for cell population ",pop))
                    }, finally = { dev.off(); }) 


                    if(summarize){
                        tryCatch({
                            print("      plotting infiltration summary (STACKED)")
                            ## plot summary over all FOV
                            sia <- ia %>% filter(Sample == s)
                            sba <- ba %>% filter(Sample == s)
                            ssum <- suppressMessages(getInfiltrationDensitySummary(densityMarkers, ia=sia, ba=sba))
                            if(length(which(!is.na(ssum$TotalDensity))) == 0){
                                print(paste0("Could not plot summary for sample ",s,", marker set ",ms))
                            } else {
                                onefile=TRUE
                                if(length(unique(ssum$CellType)) == 1){ onefile <- FALSE }
                                pdfFile <- file.path(sumOutDir,paste0(ms,"_",s,"_infiltration_density_summary.pdf"))
                                pdf(pdfFile, height=8.5, width=11, onefile=onefile)
                                print("      plotting infiltration summary (STACKED)")
                                ## plot summary over all FOV
                                plotTitle <- paste0(s, " Interface All FOV: ",unique(msCfg$MarkerSetAlias)," (Cell ID Markers = ",idMarkerTag,")")
                                plotInfiltrationDensity(ssum, densityMarkers, bandWidth, ctClrs, plotTitle=plotTitle, yMax=yMax,
                                                        separateLegend=TRUE, printLegend=TRUE, facetByFOV=FALSE,
                                                        facetByCellType=FALSE, yCol="TotalDensity",cellTypeLabels=ctLabels)
                                if(length(unique(ssum$CellType)) > 1){
                                    ssumPct <- suppressMessages(getDensityPercentage(ssum, by=c("Sample","Band")))
                                    print("      plotting intiltration summary percentages (STACKED)")
                                    plotInfiltrationDensity(ssumPct, densityMarkers, bandWidth, ctClrs, sampleOrder=sampleOrder,
                                                            plotTitle=plotTitle, yMax=yMax, separateLegend=TRUE,
                                                            printLegend=FALSE, yCol="Percent", facetByFOV=FALSE,
                                                            facetByCellType=FALSE, yAxisTitle="Percent Total Density",
                                                            cellTypeLabels=ctLabels, pct=TRUE)
                                    }
                            }
                        }, err = function(){
                            warning(paste0("Could not plot summary for sample ",s,", marker set ",ms))
                        },finally = { dev.off(); }) 
                    }
                }
            }
        } ## end for each population
    } ## end for each marker set

}

