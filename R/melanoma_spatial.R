#' Generate data frame of random points (is it random???)
#' 
#' Generate a data frame of random points that fall within the plot area
#' 
#' @param   nGrid   number of points to generate
#' @param   bbG     a list containing X0,X1,Y0,Y1 representing the boundary
#'                  of the plot area
#' @return data frame of X and Y values of the random points
generateSobolGrid<-function(nGrid,bbG) {
    gg=sobol(nGrid,2,scrambling=2)
    data.frame(
        X=gg[,1]*(bbG$X1-bbG$X0)+bbG$X0,
        Y=gg[,2]*(bbG$Y1-bbG$Y0)+bbG$Y0
        )
}



#' Compute band areas for each interface bin
#' 
#' Given a list of distances from a tumor boundary, or 'interface bins', calculate area
#' for each bin using a sobol grid
#' 
#' @param tumorBoundariesList   a list where each element is a tibble containing tumor boundaries
#'                              returned from readHaloAnnotations()
#' @param bbData                a list containing X0,X1,Y0,Y1 representing the boundary
#'                              of the trimmed FOV
#' @param interfaceBins         vector of distances (integers) that define each band; default = (-20:20)*10 
#' @param nGrid                 number of random points to use for area calculations; default=100
#' @param maxSeg                maximum segments to divide long boundary segments into; default=1000
#' @param exclusionsB           a list where each element is a tibble containing exclusion boundaries
#'                              returned from readHaloAnnotations()
computeBandAreasMC <- function(tumorBoundariesList,bbData,interfaceBins,nGrid=100,
                               maxSeg=1000,exclusionB) {

    ## generate random points all across the plot
    flog.debug("generating sobol grid")
    gg <- generateSobolGrid(nGrid,bbData)
    if(len(exclusionB)>0) {
        flog.debug("There are %s exclusion boundaries. Removing excluded points",len(exclusionB))
        excludedPts <- pointsInsidePolygonList(gg,exclusionB)
        gg <- gg[!excludedPts,]
    }

    flog.debug("trimming boundaries")
    tumBTrim <- tumorBoundariesList %>%
        map(trimDFToBB,bbData) %>%
        map(subDivideLongSegments,maxSeg)

    flog.debug("getting X,Y points for tumor boundaries")
    tumorBoundariesMergeXY <- tumBTrim %>% bind_rows %>% select(X,Y) %>% as.data.frame

    flog.debug("finding distance from points to interface")
    gg$Z <- findDistanceFromPointsToInterfacePoint(gg,tumorBoundariesMergeXY)

    flog.debug("assign negative distance to points inside tumor")
    insideTumor <- pointsInsidePolygonList(gg,tumorBoundariesList)
    gg$Z[insideTumor] <- -gg$Z[insideTumor]

    flog.debug("getting area table")
    areaTable <- table(cut(gg$Z*pixel2um,interfaceBins))

    flog.debug("making bandAreas dataframe")
    bandAreas <- as.data.frame.table(p2tomm*areaBB(bbData)*areaTable/nGrid)
    colnames(bandAreas) <- c("Band","Area")
    return(bandAreas)
}



#' Calculate area of interface bands, each defined as the collection of points that are located
#' between X and Y microns of the tumor boundary [ TO DO: REWORD THIS? ]
#' 
#' Calculate area of each band, defined as the collection of cells that are located
#' between X and Y microns of the tumor boundary [ TO DO: REWORD THIS? ]
#' **Do this for 10^maxG and for 10^maxG-1
#' 
#' @param allBoundaries  halo boundaries after removal of contained boundaries (list returned by 
#'                       cleanBoundaries())
#' @param maxG           maximum factor of 10 to use for generating random points on grid [ REWORD??? ]
#' @param bbData         a list containing X0,X1,Y0,Y1 representing the boundary
#'                       of the trimmed FOV
#' @param interfaceBins  vector of distances (integers) that define each band; default = (-20:20)*10
#' @return a dataframe with a row for each interface bin and a column for each area calculation (one for 
#'         10^maxG and one for 10^maxG-1)
calculateBandAreas <- function(allBoundaries,maxG=5,bbData,interfaceBins) {

    if(is.null(interfaceBins)){
        interfaceBins <- (-20:20)*10
    }
    tumB <- allBoundaries$tumB
    excB <- allBoundaries$excB

    bandArea <- list()
    for(nGrid in 10^((maxG-1):maxG)) {
        maxSeg <- max(unlist(bbData))/sqrt(10*nGrid)
        flog.debug("nGrid = %s  maxSeg = %s ",nGrid,maxSeg)
        bandArea[[as.character(nGrid)]] <- computeBandAreasMC(tumB,bbData,interfaceBins,
                                                               nGrid,maxSeg,excB)
    }

    aa <- bind_cols(bandArea) %>% select(matches("Area"))
    rownames(aa) <- bandArea[[1]][,1]

    aa

}



#' Compute lattice of boundary bins
#'
#' Given an annotation file and halo infiltration data, compute tibble of Sample, FOV, 
#' Halo Version, MaxG, Band, Area, Halo
#'
#' *NOTE: formerly computeBoundaryBinsLattice
#' 
#' @param aFile                Halo boundary annotation file in XML format
#' @param sample               sample name
#' @param fov                  FOV
#' @param maxG                 maximum factor of 10 to use for generating random points on grid [ REWORD??? ]
#' @param writeCSVfiles        write lattice band area file in csv format; default=FALSE
#' @param outDir               if writeCSVfiles=TRUE, files will be written here
#' @param bbData               a list containing X0,X1,Y0,Y1 representing the boundary
#'                             of the trimmed FOV
#' @param interfaceBins        vector of distances (integers) that define each band; default = (-20:20)*10 
#' @return tibble 
getAreaPerBand <- function(allBoundaries, sample, fov, maxG, outDir=NULL,
                                       bbData=NULL, writeCSVfiles=FALSE, interfaceBins=NULL){ 

    if(is.null(interfaceBins)){
        interfaceBins <- (-20:20)*10
    }

    SAMPLE <- sample
    sampleSPOT <- as.numeric(fov)
    flog.debug("%s SPOT %s",SAMPLE,sampleSPOT)
    aa=calculateBandAreas(allBoundaries,maxG=maxG,bbData,interfaceBins)

    flog.debug("creating tibble")
    xx <- as.tibble(data.frame(
            Sample=SAMPLE,
            Spot=sampleSPOT,
            Version="v3.4",
            MaxG=maxG,
            Band=rownames(aa), 
            Area=aa[,1] 
            )
          )

    if(writeCSVfiles){
        flog.debug("writing csv files")
        if(!is.null(bbData)){
            maxSeg=max(unlist(bbData))/sqrt(10^(maxG+1))
            csvfile <- file.path(outDir,cc("latticeBandAreaV4Exc",SAMPLE,sampleSPOT,"maxG",maxG,round(maxSeg,2),".csv"))
        } else {
            csvfile <- file.path(outDir,cc("latticeBandAreaV4Exc",tag,"maxG",maxG,".csv"))
        }
        write_csv(xx,csvfile)
    }
    xx
}


#' Add columns for band area to existing tibble with a column of bands
#'
#' Add columns for band area to existing tibble with a column of bands
#'
#' @param bdat      tibble of band data, including at least a column named Band
#' @param bandArea  a data frame or tibble with an area value for each band in bdat
#' @return tibble containing all data in bdat and bandArea
joinBandArea<-function(bdat,bandArea) {

    bdat$Band=as.character(bdat$Band)
    flog.debug("adding band area to band counts")
    xx=full_join(bdat,bandArea)
    flog.debug("filling in NAs with 0s")
    xx$n[is.na(xx$n)]=0
    xx$Band=factor(xx$Band,levels=bandArea$Band)
    flog.debug("arranging by Band")
    xx %<>% arrange(Band)
    xx

}

#' Calculate total area and density of entire FOVs
#' 
#' Calculate total area and density of FOVs that do not contain 
#' tumors
#' 
#' @param aFiles              a vector of Halo boundary annotation files, with 
#'                            '.annotations' extension in XML format
#' @param data                Halo data for one sample loaded from a *.rda 
#' @param markerSet           a vector markers to calculate density for
#' @param cellTypeName        a name to represent the cell types in marker file
#' @param funcMarker          this marker will be added to all markers in markerSet in order to
#'                            compare cells that are +/- for this functional marker [ TO DO: REWORD THIS ]
#' @param fovBB               a list with four elements: X0, X1, Y0, Y1, representing the minimum
#'                            and maximum coordinates of the FOV to be plotted
#' @param pad                 amount to trim from FOV before calculating
#' @param maxG                the maximum factor of 10 to use for generating random points for area calculation
#' @param writeCSVfiles       logical indicating whether to write density and area tables to file; 
#' @param outDir              if writeCSVfiles=T, directory to which these files will be written
#' @return a list containing two tables: 'density' and 'area' 
#' @export
calculateFOVStats <- function(aFiles, data, markerSet, cellTypeName,
                            pad, maxG, writeCSVfiles=TRUE,
                            funcMarker=NULL,sortByMarker=NULL,outDir){
  pad <- as.numeric(pad)
  dd <- data
  sampleName <- gsub("_ObjectAnalysisData","",unique(dd$Sample)) 

  rho <- list()
  areas <- list()

  flog.debug("getting spots")
  spots <- dd %>% dplyr::select(SPOT) %>% distinct(SPOT) %>% arrange(SPOT) %$% as.vector(SPOT)
  for(spot in spots){
    flog.debug("working on spot %s",spot)

    flog.debug("checking for annotations file...")
    sTag <- cc(sampleName, paste0("Spot",spot,".annotations"))
    aai <- grep(sTag, aFiles)
    tumB <- NULL
    excB <- NULL

    if(len(aai) > 0){
      flog.debug("File found. Getting tumor and/or exclusion boundaries.")
      aFile <- aFiles[aai]
      boundaries <- readHaloAnnotations(aFile)

      flog.debug("removing contained boundaries")
      ## remove boundaries that are completely contained in another one
      allBoundaries <- cleanBoundaries(boundaries)
      tumB <- allBoundaries$tumB
      excB <- allBoundaries$excB
    }

    if(len(tumB) > 0){
      flog.debug("SPOT %s contains tumor boundaries. skipping.",spot)
      next ## FOV with tumor boundaries are to be handled by calculateInterfaceArea()      
      ## QUESTION: do we ever want to calculate total density of FOV even if there ARE
      ## tumor boundaries? 
    }

    ds <- dd %>% filter(SPOT==spot & ValueType=="Positive")
    bbFOV <- getBoundingBox(ds)
    bbData <- bbFOV
    if(pad > 0){
      bbData <- padBoundingBox(bbFOV,-pad/pixel2um)
    }

    ds <- convertCellMinMaxToMidpoints(ds) ## this can be done in removeExcludedPoints too 
                                           ## but we need it done even if there arent any
    ds <- trimDFToBB(ds,bbData) ## needs X and Y coords

    ## calculate area of entire FOV
    gg <- generateSobolGrid(10^maxG,bbData)
    if(len(excB) > 0){
      ## remove excluded points
      ds <- removeExcludedPoints(ds,excB=excB)
      excludedPts <- pointsInsidePolygonList(gg,excB)
      gg <- gg[!excludedPts,]
    }

    markers <- ds %>% distinct(Marker) %$% as.character(Marker)
    ds <- ds %>% spread(Marker,Value)
    for(mi in markers) {
      ds[[paste0(mi,"-")]] <- ifelse(ds[[mi]]==0,1,0)
    }

    ## get total number of cells
    mt <- computeMultiMarkerTable(ds, markerSet) %>% filter(CD3==1) ## <- figure out how to make this dynamic
    mtCounts <- mt %>% select(markerSet) %>% colSums
    totalCounts <- as.matrix(mtCounts)

    totalArea <- p2tomm * areaBB(bbData) * nrow(gg)/10^maxG
    MM <- 1/totalArea

    rho[[len(rho)+1]] <- data.frame(CellType=rownames(totalCounts),
                                    Counts=totalCounts,
                                    Total=totalCounts/totalArea,
                                    Sample=sampleName,
                                    SPOT=spot) %>%
                         as.tibble

    areas[[len(areas)+1]] <- data.frame(Area=totalArea) %>%
                                as_tibble %>%
                                mutate(Sample=sampleName,SPOT=spot)
  }
  rhos <- NULL
  areass <- NULL
  flog.debug("Num elements in list rho: %s",length(rho))
  if(length(rho) > 0){
    rhos <- rho %>% bind_rows
  }
  flog.debug("Num elements in list areas: %s",length(areas))
  if(length(areas) > 0){
    areass <- areas %>% bind_rows 
  }

  if(writeCSVfiles){
    write_csv(as.data.frame(areass), file.path(outDir,cc(sampleName, "_area.csv")))
    write_csv(as.data.frame(rhos), file.path(outDir,cc(sampleName, "_density.csv")))
  }

  return(list(density=rhos, area=areass))

}      

#' Calculate total area and density of tumor interface 
#' 
#' Calculate area and density of cells within a given distance to a tumor boundary
#'
#' @param aFiles              a vector of Halo boundary annotation files, with 
#'                            '.annotations' extension in XML format
#' @param data                Halo data for one sample, loaded from *.rda file  
#' @param pad                 amount to trim from FOV before calculating
#' @param writeCSVfiles       logical indicating whether to write density and area tables to file; 
#'                            default=TRUE
#' @param maxG                the maximum factor of 10 to use for generating random points for area calculation
#' @param outDir              if writeCSVfiles=T, directory to which these files will be written
#' @param interfaceBins       a vector of distances from tumor interface in which cells will be binned; 
#'                            default=(-20:20)*10
#' @param statsByBand         logical; return stats broken down by interface bands (bins) as opposed to total
#'                            stats; default=FALSE 
#' @return a list containing two tables: 'density' and 'area' 
#' @export
calculateInterfaceArea <- function(aFiles, dat, pad, writeCSVfiles=TRUE,
                                   maxG=5,outDir=NULL,statsByBand=FALSE,
                                   interfaceBins=(-20:20)*10){
  pad <- as.numeric(pad)

  dd <- dat 
  sampleName <- gsub("_ObjectAnalysisData","",unique(dd$Sample)) 

  areas <- list()
  bandAssignments <- list()

  flog.debug("getting spots")
  spots <- dd %>% dplyr::select(SPOT) %>% distinct(SPOT) %>% arrange(SPOT) %$% as.vector(SPOT)
  for(spot in spots){
    flog.debug("working on spot %s",spot)

    flog.debug("filtering dataset for spot..")
    ds <- dd %>% filter(SPOT==spot & ValueType=="Positive")

    bbFOV <- getBoundingBox(ds)
    bbData <- bbFOV
    if(pad > 0){
      bbData <- padBoundingBox(bbFOV,-pad/pixel2um)
    }     

    sTag <- cc(sampleName, paste0("Spot",spot,".annotations"))
    aai <- grep(sTag, aFiles)
    if(len(aai) > 0){
      aFile <- aFiles[aai]

      flog.debug("getting halo boundaries")
      boundaries <- readHaloAnnotations(aFile)

      flog.debug("removing contained boundaries and separating tumors and exclusions")
      ## remove boundaries that are completely contained in another one
      allBoundaries <- cleanBoundaries(boundaries)
      tumB <- allBoundaries$tumB
      excB <- allBoundaries$excB

      if(len(tumB) == 0){
        next
      }
      ##
      ## get Distance and Band assignments for all points that fall inside a tumor
      ##
      flog.debug("adding to data distance from tumor boundary and band assignment")
      ba <- getBandAssignments(ds, bbData, allBoundaries, interfaceBins)
      bandAssignments[[len(bandAssignments)+1]] <- ba
     
      flog.debug("filtering band count data for CD3+ and counting cells in each band")
      bdat <- ba %>% filter(CD3==1) %>% filter(!is.na(Band)) %>% count(Band) %$% data.frame(Band,n)
    
      flog.debug("calculating area of each band") 
      bandArea <- getAreaPerBand(allBoundaries,sampleName,spot,maxG=maxG, 
                                 outDir=outDir, bbData=bbData,
                                 interfaceBins=interfaceBins)
      if(statsByBand){
        areas[[len(areas)+1]] <- bandArea %>% select(Band,Area,Sample,SPOT=Spot)
      } else {
        flog.debug("calculating area per band and then summing for total area of space 200um around boundaries")
        bdat <- joinBandArea(bdat, bandArea)

        flog.debug("calculating total area for bands each inside and outside tumor")
        bandAreaCollapse <- tapply(bdat$Area,interfaceSide(bdat$Band),sum)

        flog.debug("getting total area for ALL bands")
        totalArea <- sum(bandAreaCollapse)

        flog.debug("storing area")
        areas[[len(areas)+1]] <- data.frame(Area=bandAreaCollapse) %>%
                                      rownames_to_column("Side") %>%
                                      as_tibble %>%
                                      mutate(Sample=sampleName,SPOT=spot)
      }
    }
  }
  areass <- NULL
  allBandAssignments <- NULL
  flog.debug("Num elements in list areas: %s",length(areas))
  if(length(areas) > 0){
    if(statsByBand){ 
      areass <- areas %>% bind_rows %>% spread(Band,Area) 
    } else {
      areass <- areas %>% bind_rows %>% spread(Side,Area) %>% mutate(Area=Inside+Outside)
    }
  }
  if(length(bandAssignments)>0){
    allBandAssignments <- bandAssignments %>% bind_rows  
  }

  if(writeCSVfiles){
    write_csv(as.data.frame(areass), file.path(outDir,cc(sampleName, "Interface_Area.csv")))
    write_csv(as.data.frame(allBandAssignments), file.path(outDir,cc(sampleName, "Band_Assignments.csv")))
  }

  return(list(area=areass, bandAssignments=allBandAssignments))
}

#' Remove exclusion boundaries that are contained in another one
#' 
#' Remove data from boundaries table that represent exclusion boundaries
#' that are completely surrounded by another exclusion boundary
#' 
#' @param boundaries            table generated by readHaloAnnotations
#'                             
cleanBoundaries <- function(boundaries){

    regionTable <- boundaries %>% bind_rows %>% count(RegionNum,RegionCode)
    excB <- boundaries[regionTable$RegionCode!="Tum"]
    tumB <- boundaries[regionTable$RegionCode=="Tum"]

    if(len(excB) > 1){
        ## make list of spatial polygons, one element for each exclusion boundary
        spExcB <- seq(excB) %>% map(function(b){boundaryToSpatialPolygon(excB[[b]],b)})
        
        ## are any contained within another?
        containedBoundary <- rep(FALSE,len(excB))
        for(i in seq(len(excB))){
            containedBoundary[i] <- spExcB[-i] %>% map(gContains,spExcB[[i]]) %>% unlist %>% any
        }
        ## remove those
        excB <- excB[!containedBoundary]
    }
    return(list(excB=excB,tumB=tumB))
}

#' Remove points from data set that fall inside exclusion boundaries
#' 
#' Read Halo boundaries annotations file and remove from data set any cells
#' that fall inside those regions
#' 
#' @param ds     tibble containing at minimum, Sample, SPOT, UUID, X, Y, Marker, Value
#' @param excB   a list of excluded region tables as returned by cleanBoundaries;
#'               if NULL, must provide the original Halo annotations file in XML format
#' @param aFile  Halo boundaries annotation file in XML format; if NULL, must provide
#'               the list of excluded region tables as returned by cleanBoundaries()
#' @return tibble filtered to remove points that fall inside exclusion boundaries
removeExcludedPoints<-function(ds,aFile=NULL,excB=NULL) {

    if(!all(c('X','Y') %in% colnames(ds))){
        ds <- convertCellMinMaxToMidpoints(ds)
    }

    if(is.null(excB)){
        boundaries <- readHaloAnnotations(aFile)
        allBoundaries <- cleanBoundaries(boundaries)
        excB <- allBoundaries$excB
    }

    spCells <- ds %>% dplyr::select(UUID,X,Y)
    coordinates(spCells) <- ~X+Y

    if(len(excB)>0){
        spExcB <- seq(excB) %>% map(function(b){boundaryToSpatialPolygon(excB[[b]],b)})
        excludedCellIdx <- NULL
        for(jj in len(spExcB)) {
            excCellJJ <- unlist(over(spExcB[[jj]],geometry(spCells),returnList=T))
            excludedCellIdx <- union(excludedCellIdx,excCellJJ)
        }
        if(len(excludedCellIdx)>0){
            ds <- ds[-excludedCellIdx,]
        }
    }

    return(ds)

}

#' Get counts of cells that fall within given distances of
#' tumor boundaries for a single FOV
#'
#' Given a set of interface bands, or set distances from a tumor 
#' boundary, count the numbers of cells that fall within those bands
#' 
#' @param  ds             a tibble containing data for one FOV, including columns
#'                        for X, Y, and one for each marker
#' @param bbData          a list containing X0,X1,Y0,Y1 representing the boundary
#'                        of the trimmed FOV
#' @param allBoundaries   a list returned from readHaloAnnotations
#' @param interfaceBins   vector of distances (integers) that define each band; default = (-20:20)*10 
getBandAssignments<-function(ds,bbData,allBoundaries,interfaceBins=(-20:20)*10){

    flog.debug("getting boundaries")
    tumB <- allBoundaries$tumB
    excB <- allBoundaries$excB
    if(is.null(tumB) || len(tumB) == 0){
        flog.info("No TUMOR boundaries found in file %s. Skipping.")
        return()
    }

    markers <- ds %>% distinct(Marker) %$% as.character(Marker)

    flog.debug("getting midpoint of each cell")
    ## get midpoint of each cell
    ds <- ds %>%
          mutate(X=(XMax+XMin)/2,Y=-(YMax+YMin)/2) %>%
          dplyr::select(Sample,SPOT,UUID,X,Y,Marker,Value) %>%
          spread(Marker,Value)

    flog.debug("filtering out cells that fall outside bounding box")
    ## filter out cells that fall outside bounding box
    ds <- trimDFToBB(ds,bbData)

    ## add values for negative markers
    for(mi in markers) {
      ds[[paste0(mi,"-")]]=ifelse(ds[[mi]]==0,1,0)
    }

    ## create a SpatialPointsDataFrame
    spCells <- ds %>% dplyr::select(UUID,X,Y)
    coordinates(spCells) <- ~X+Y

    flog.debug("removing cells in exclusion boundaries")
    ## if there are still any exclusion boundaries, remove the cells that fall
    ## within those boundaries?
    if(len(excB)>0){
        spExcB <- seq(excB) %>% map(function(b){boundaryToSpatialPolygon(excB[[b]],b)})
        excludedCellIdx <- NULL
        for(jj in len(spExcB)) {
            excCellJJ <- unlist(over(spExcB[[jj]],geometry(spCells),returnList=T))
            excludedCellIdx <- union(excludedCellIdx,excCellJJ)
        }
        #points(spCells[excludedCellIdx,],col=3,cex=.7)
        if(len(excludedCellIdx)>0){
            ds <- ds[-excludedCellIdx,]
        }
        flog.debug("updating spatialpointsdataframe")
        ## update the SpatialPointsDataFrame
        spCells <- ds %>% dplyr::select(UUID,X,Y)
        coordinates(spCells) <- ~X+Y
    }

    flog.debug("getting tumor contour")
    ## filter data for cells that are near a tumor boundary AND inside
    ## the bounding box 
    tumorContour <- tumB %>% map(trimDFToBB,bbData)
    flog.debug("merging X and Y")
    ## get x and y coords of all remaining tumor cells
    tumorBoundariesMergeXY <- tumorContour %>% bind_rows %>% dplyr::select(X,Y) %>% as.data.frame

    flog.debug("finding distances closest to boundary")
    ## for each point, find distance to closes boundary point
    dsPt <- ds %>% dplyr::select(X,Y) %>% as.data.frame
    dTumor <- findDistanceFromPointsToInterfacePoint(
        dsPt,
        tumorBoundariesMergeXY
    )
    ## list of true/false indicating whether point is in a tumor
    tumorPoints <- pointsInsidePolygonList(dsPt,tumB)
    ## assign neg to distances for points inside tumor
    dTumor <- ifelse(tumorPoints,-dTumor,dTumor)
    ## convert to um
    dTumor <- dTumor*pixel2um
    ds$Distance <- dTumor
    ds$Band <- cut(dTumor,interfaceBins) ## assign each point to a bin

    return(ds)
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

#' Add to data MAD min/max y-values for plotting median error bars
#' 
#' For each SPOT=="Median" row in tibble, add MAD of densities for
#' all SPOTs for a specific Sample+CellType. Then add ymin and ymax values
#' representing total median density +/- MAD
#'
#' @param dat  tibble including at least columns for Sample,CellType & Density;
#'             Also must include rows with SPOT=="Median" representing median of
#'             all spots for each Sample+CellType 
#' @return original tibbles with added columns for TotalDensity, mad, ymin, ymax
addMedianErrorBarValues <- function(dat){
    ## get table of MAD values for Medians only
    md <-  dat %>%
              filter(SPOT != "Median") %>%
              group_by(Sample, CellType) %>%
              summarise(mad = mad(Density)) 

    ## Density in dat is now divided by FuncPos and FuncNeg, 
    ## so need to get total density to calc ymin and ymax 
    ## of error bars
    tmp <- filter(dat,SPOT=="Median") %>%
            group_by(Sample,SPOT,CellType) %>%
            summarise(TotalDensity = sum(Density)) %>%
            full_join(md) %>%
            mutate(ymin=TotalDensity-mad, ymax=TotalDensity+mad)

    ## add values to original data
    full_join(dat, tmp)

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
                     axis.text.x = element_text(angle=45, hjust=1),
                     axis.ticks = element_blank(),
                     strip.text.x = element_text(face="bold", size=12, margin=margin(.4, .2, .4, .2, "cm")),
                     strip.text.y = element_text(margin = margin(.4, .2, .4, .2, "cm")),
                     strip.placement = "outside",
                     plot.background = element_blank(),
                     panel.background = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.line = element_line(color = "black", size=0.5),
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

#' Plot density for each marker, broken down by FOV
#' 
#' Bar plot(s) showing density of a marker in each FOV of one sample
#'
#' @param markerList          vector of markers/marker combinations to be plotted; there will
#'                            be one plot per marker
#' @param density             tibble containing columns for Sample, SPOT, CellType, Total
#' @param funcMarker          a functional marker that is included in the marker combos in markerTable;
#'                            for each marker, there should be a combo that includes the functional marker
#'                            (indicating that cell type is positive for this marker) and one that does
#'                            not (indicating that the cell type may or may not be positive for this marker);
#'                            default=NULL;
#' @param sampleColor         all bars will be colored gray except for the Median bar, which will be
#'                            this color; default=#f16913
#' @param sampleColorDark     a darker version of sampleColor (or really any color of your choosing)
#'                            that will represent the median of cells positive for the funcMarker;
#'                            default=#80380A
#' @param spotOrder           a vector of FOV numbers indicating the order in which FOVs should be in 
#'                            the plot(s); default=NULL
#' @param sampleOrder         a vector indicating the order in which samples should appear
#' @param logPlot             logical indicating that density should be plotted on log scale; default=FALSE
#' @param boxPlots            logical indicating whether to create a box plot for each cell type showing 
#'                            the mean difference between samples; default=TRUE
#' @param consistentYScaling  logical; when onePanel=FALSE, set this parameter to FALSE in order to scale separate 
#' @param pdfFile             if given, plot(s) will be written to this file; default=NULL
plotTotalDensity <- function(markerList, rho, funcMarker=NULL, funcPosColor="grey35",
                             funcNegColor="grey50", sampleColor="#f16913", sampleColorDark="#80380A", 
                             spotOrder=NULL, sampleOrder=NULL, plotPer='markerSet', pdfFile=NULL,
                             sortByMarker=NULL, includeSortByMarker=FALSE, boxPlots=TRUE,
                             consistentYScaling=TRUE){

    #if(boxPlots){
    #    if(!is.null(pdfFile)){
    #        bpdfFile <- gsub(".pdf","_boxplots.pdf",pdfFile)
    #        pdf(bpdfFile,height=8,width=8)
    #    }

    #    densityBoxPlots(markers, rho, sampleOrder)

     #   if(!is.null(pdfFile)){
     #      dev.off()
     #   }
    #}


    #########
    ## prep data for plotting
    #########

    ## spread density table to get to get a column for each cell type
    mt <- rho %>%
          dplyr::select(CellType,Sample,SPOT,Total) %>%
          spread(CellType,Total)
    mt$SPOT <- as.character(mt$SPOT)

    ## get median of all spots for each sample+cellType
    meds <- gather(mt, 3:ncol(mt), key="CellType", value="Density") %>% 
               group_by(Sample,CellType) %>% 
               summarize(Median=median(Density)) %>% 
               spread(CellType, Median) %>%
               mutate(SPOT="Median")
    mt <- mt %>% bind_rows(meds)

    ## add column to indicate whether sample/spot/celltype is
    ## positive or negative for the functional marker, and calculate
    ## density of the negative ones
    if(!is.null(funcMarker)){
        fm <- paste0(",",funcMarker)
        fmp <- paste0(funcMarker,"+")
        fmn <- paste0(funcMarker,"-")
        dat <- gather(mt, 3:ncol(mt),key="CellType",value="TotalDensity") %>% 
                  mutate(FunctionalMarker=ifelse(grepl(fm,CellType),fmp,fmn),
                         CellType=gsub(fm,"",CellType),
                         SampleSpot=paste(Sample,SPOT,sep="_")
                        ) %>% 
                  as.tibble %>% 
                  spread(FunctionalMarker, TotalDensity) %>% 
                  replace(is.na(.),0)
        dat[[fmn]] <- dat[[fmn]] - dat[[fmp]]
        dat <- gather(dat, -(Sample:SampleSpot), key="Functional",value="Density") %>%
               arrange(desc(Functional))
    } else {
        ## manipulate mt to have Functional col that's all the same??
        dat <- gather(mt, 3:ncol(mt), key="CellType", value="Density") %>%
                mutate(Functional="X")
    }

    ## sort samples
    if(!is.null(sampleOrder)){
        flog.debug("Sorting samples in this order: %s", paste(sampleOrder, collapse=", "))
        dat$Sample2 <- factor(dat$Sample, levels=sampleOrder)
    }

    ## sort FOVs by density of a particular marker, if given
    if(!is.null(sortByMarker)){
        flog.debug("Sorting FOVs by density of %s",sortByMarker)
        spotOrder <- c(filter(dat, SPOT!="Median", 
                                 CellType==sortByMarker, 
                                 grepl("-",Functional)) %>% 
                                 arrange(Sample,Density) %>% 
                                 pull(SampleSpot),
                       unique(filter(dat,SPOT=="Median") %>%
                          pull(SampleSpot)))
        dat$SampleSpot2 <- factor(dat$SampleSpot, levels=spotOrder) 

        if(!includeSortByMarker){
            flog.debug("removing ",includeSortByMarker," from data")
            dat <- dat %>% filter(CellType != sortByMarker)
        }
    }

    ## add columns for error bar values
    flog.debug("adding median error bars")
    dat <- addMedianErrorBarValues(dat)

    ## add column for color fill
    dat$Colors <- getTotalDensityFillColors(dat,funcPosColor,funcNegColor,sampleColor,sampleColorDark)

    ## add new column for cell type labels in order to make extra long ones wrap
    cellTypeLabels <- str_wrap(gsub(",",", ",unique(dat$CellType)), width=25)
    dat$CellTypeLabels <- gsub(", ",",",cellTypeLabels)


    ##########
    #### FINALLY, MAKE PLOTS!
    ##########

    plotList <- getTotalDensityPlots(dat, markerList, plotPer, funcMarker=funcMarker, funcPosColor=funcPosColor,
                     funcNegColor=funcNegColor, sampleColor=sampleColor, sampleColorDark=sampleColorDark,
                     consistentYScaling=consistentYScaling)

    return(plotList)
}


#' Plot density broken down by distance from tumor interface 
#' 
#' Given a vector of markers and a tibble containing columns for Density and 
#' Band, generate a bar plot for each FOV, with one bar per Band, colored by 
#' marker
#'
#' @param density    tibble with columns for Density, Band and CellType
#' @param sampleName name of sample from which the FOVs came from
#' @param markers    vector of markers/marker combos
#' @param pdfFile    pdfFile to which plot will be printed
plotDensityByBand <- function(density,sampleName,markers,pdfFile=NULL){
    if(!is.null(pdfFile)){
        pdf(pdfFile,height=8.5,width=11)
    }

    plotTheme <- theme(legend.position="bottom",
                    legend.title=element_blank(),
                     axis.title.x=element_blank(),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())

    g <- ggplot(density, aes(x=Band,y=Density)) +
              geom_bar(aes(fill = CellType),
                       stat="identity",
                       position=position_stack()) +
              plotTheme +
              facet_wrap(~SPOT) +
              ggtitle(sampleName) +
              ylab("Density (per mm^2)") +
              scale_fill_manual(values = getMarkerColors(markerColors=NULL,markerCombos=markers)) +
              geom_vline(aes(xintercept = 20.5), linetype=3)
    print(g)
    if(!is.null(pdfFile)){
        dev.off()
    }
}

calculateDensity <- function(dat, areas, bandAssignments, markerSet, cellTypeName,
                              pad, funcMarker=NULL, sortByMarker=NULL,writeCSVfiles=TRUE,
                              outDir=NULL,statsByBand=FALSE){
  pad <- as.numeric(pad)

  dd <- dat
  sampleName <- gsub("_ObjectAnalysisData","",unique(dd$Sample))

  fillNA <- list()
  dummy <- lapply(markerSet, function(x){ fillNA[x] <<- 0 })

  rho <- list()

  flog.debug("computing multimarker table and summarizing each band")
  mt <- computeMultiMarkerTable(bandAssignments,markerSet)
  mtCounts <- mt %>%
              group_by(SPOT,Band) %>%
              summarize_at(markerSet,sum) %>%
              complete(Band,fill=fillNA)

  for(spot in unique(areas$SPOT)){
    flog.debug("working on spot %s",spot)

    if(statsByBand){
      ## bandArea = ???
      MM <- 1/bandArea$Area
      densityByMarker <- as.data.frame(t(mtB * MM))
      rho[[len(rho)+1]] <- densityByMarker %>%
                            cbind(data.frame(Counts=t(mtB),check.names=FALSE)) %>%
                            rownames_to_column("CellType") %>%
                            mutate(Sample=sampleName,SPOT=spot)

    } else {

      flog.debug("counting number of cells inside and outside the tumor for each marker")
      
      mtB <- mtCounts %>% filter(SPOT==spot,!is.na(Band))
      rn <- mtB$Band
      mtB <- as.matrix(mtB[3:ncol(mtB)],check.names=F)
      rownames(mtB) <- rn
      mtBCollapse <- apply(mtB,2,function(x){tapply(x,interfaceSide(rownames(mtB)),sum)})

      interfaceArea <- areas %>% filter(SPOT==spot)

      ## divide into 1, to later be multiplied by counts to calculate density
      inoutArea <- interfaceArea %>% select(Inside,Outside)
      MM <- diag(1/inoutArea)
      rownames(MM) <- c("Inside","Outside")
      colnames(MM) <- c("Inside","Outside")

      flog.debug("getting total counts and area for ALL bands")
      totalCounts <- colSums(mtBCollapse)
      flog.debug("getting total area for ALL bands")
      totalArea <- interfaceArea$Area

      flog.debug("calculating density")
      rho[[len(rho)+1]] <- (t(mtBCollapse) %*% MM) %>%
                              as.data.frame %>%
                              cbind(data.frame(Counts=t(mtBCollapse))) %>%
                              rownames_to_column("CellType") %>%
                              mutate(Total=totalCounts/totalArea,Sample=sampleName,SPOT=spot)
    }
  }

  rhos <- NULL
  flog.debug("Num elements in list rho: %s",length(rho))
  if(length(rho) > 0){
    rhos <- rho %>% bind_rows
  }

  if(writeCSVfiles){
    write_csv(as.data.frame(rhos), file.path(outDir,cc(sampleName, "interface_density.csv")))
  }

  return(rhos)
}

#' Internal utility function to print total density plots to pdf
#' 
#' Print plots either all in one file, with one file per marker set, 
#' or one file per marker. One file per marker can only be done if there
#' is no facetting involved (i.e., singlePanelPlots must be set to FALSE)
#' 
#' @param plotList     a list of ggplots to be printed; may or may not be a nested list,
#'                     depending on how the plots were built (specifically on, 'plotPer') 
#' @param markerList   named list of marker sets
#' @param filePer      print one file per this value; choices are ['plot'|'sample'|'marker'|NULL]; 
#'                     set to NULL to print all plots in one file, named by pdfFile (required when
#'                     printing to a single file)
#' @param plotPer      one plot was built for each of these units; choies are ['markerSet','marker','sample',
#'                     'sample+marker']
#' @param pdfFile      if printing all plots to a single file (filePer is NULL), file name must be specified here 
#' @param height       height of PDF image; default=8.5
#' @param width        width of PDF image; default=11
#' @return nothing 
printTotalDensityPlots <- function(plotList, markerList, filePer=NULL, plotPer='markerSet', 
                                   pdfFile=NULL, height=8.5, width=11){
    h=height
    w=width

print(names(plotList))

    ## validate choices
    if(!is.null(filePer)){
        if(filePer %in% c('sample','marker') & !(plotPer %in% c('sample+marker','marker+sample'))){
            stop(paste0("Can not print files per ",filePer," unless plotPer is 'sample+marker' or 'marker+sample'"))
        }
    }

    if(is.null(filePer)){ 
        ## print all plots in a single file; pdfFile must be specified
        if(is.null(pdfFile)){
            stop("Must specify pdfFile when printing all plots to a single file (filePer==NULL)")
        }
        pdf(pdfFile,height=h,width=w)
        if(plotPer=='sample+marker'){
            ## plots are in a nested list, with the outer list named by sample and inner lists
            ## by markers
            for(s in seq(names(plotList))){
                for(m in seq(names(plotList[[s]]))){
                    print(plotList[[s]][[m]])
                }
            }
        } else {
            ## if plots were build in any other way, they are just in a flat list named by 
            ## the way they were wrapped (by sample, or marker, markerSet, etc.
            for(p in names(plotList)){
                print(plotList[[p]])
            }
        }
        dev.off()
    } else if(filePer == 'plot'){
        ## print one file per plot, regardless of how the plots were build
        if(plotPer %in% c('sample+marker','marker+sample')){
            for(s in names(plotList)){
                for(m in names(plotList[[s]])){
                    pdf(paste0(s,"_",m,".pdf"),height=h,width=w)
                        print(plotList[[s]][[m]])
                    dev.off()
                }
            }
        } else {
            for(pn in names(plotList)){
                pdf(paste0(pn,"_totalDensity.pdf"),height=h,width=w)
                print(plotList[[pn]])
                dev.off()
            }
        }
    } else if(filePer == 'sample'){
        ## if each plot is one marker in one sample, print all
        ## plots for one sample in one file
        for(s in names(plotList)){
            pdf(paste0(s,"_totalDensity.pdf"),height=h,width=w)
            for(m in names(plotList[[s]])){
                print(plotList[[s]][[m]])
            }
            dev.off()
        } 
    } else if(filePer == 'marker'){
        ## if each plot is one marker in one sample, print all
        ## plots for one marker in one file
        allMarkers <- c()
        for(s in names(plotList)){
            allMarkers <- c(allMarkers, names(plotList[[s]]))
        }
        allMarkers <- unique(allMarkers)
        for(m in allMarkers){
            pdf(paste0(m,"_totalDensity.pdf"),height=h,width=w)
            for(s in names(plotList)){
                if(m %in% names(plotList[[s]])){
                    print(plotList[[s]][[m]])
                }
            }
            dev.off()
        }
    } else {
        stop(paste0("Unrecognized print option: filePer='",filePer,"'"))
    }
}

#' Plot density of each FOV for one sample
#'
#' @param dataFiles                 vector of *.rda files containing Halo data for a single sample
#' @param annotationsDirs           vector of directories containing *.annotations files, XML files of boundary
#'                                  annotations from Halo; order must correspond to order of dataFiles
#' @param markerSetFiles             a vector of text files with each line containing a cell type marker; may be
#'                                  marker combinations delimited by commas
#' @param cellTypeName              a character string naming the group of cell types (e.g., T Cells)
#' @param fovBB                     a list of four elements: X0,X1,Y0,Y1, representing the min/max coords
#'                                  for FOVs
#' @param pad                       amount to trim from each FOV
#' @param pdfFile                   output PDF file
#' @param logPlot                   logical; when set to TRUE, plot density per mm^2; default=FALSE
#' @param funcMarker                character string of functional marker to plot on top of every other marker
#' @param sampleColor               character string indicating HTML color to represent this particular sample
#' @param sampleColorDark           a slightly darker shade of sampleColor, to be used for density of funcMarker+
#' @param sortByMarker              sort FOVs by density of this marker (cell type)in ascending order;
#'                                  default=CD3
#' @param exclude_sample_fov        list where names are sample names and values are vectors of FOV to be  
#'                                  removed from data before any calculations 
#' @param exclude_sample_marker     list where names are sample names and values are vectors of markers to be 
#'                                  should be removed from data before any calculations
#' @param writeCSVfiles             logical; write CSV files, one of density values and one of area values; 
#'                                  default=TRUE
#' @param maxG                      maximum factor of 10 to use for generating random points on grid [ REWORD??? ]
#' @param outDir                    if writeCSVfiles=T, directory to which these files will be written
#' @param byBand                    logical indicating whether to break down density by distance from tumor interface; 
#'                                  default=FALSE
#' @param bandWidth                 width of each band around tumor interface, to be used when byBand is TRUE; 
#'                                  default=10
#' @param maxDistanceFromInterface  include only those cells within this distance from tumor interface
#' @param boxPlots                  logical indicating whether to generate a box plot for each cell type
#'                                  comparing mean densities between samples; default=TRUE
#' @param plotPer                   generate plot per ['markerSet'|'marker'|'sample'|'sample+marker]; default='markerSet' 
#' @param fileSetup                 ['single'|'per-marker'|'per-markerSetFile']; when set to 'single', ALL plots will be in
#'                                  a single PDF file; must provide pdfFile 
#' @param consistentYScaling        logical; when onePanel=FALSE, set this parameter to FALSE in order to scale separate 
#'                                  plots individually; default=TRUE
#' @export
plotDensity <- function(dataFiles, annotationsDirs, markerSetFiles, cellTypeName, fovBB, pad, pdfFile, 
                         logPlot=F, funcMarker=NULL, funcPosColor="grey50", funcNegColor="grey35",
                         sampleColor="#f16913", sampleColorDark="#80380A",sortByMarker="CD3", 
                         sampleOrder=NULL, exclude_sample_fov=NULL, exclude_sample_marker=NULL,
                         writeCSVfiles=TRUE, maxG=5, outDir=getwd(), byBand=FALSE, bandWidth=10,
                         maxDistanceFromInterface=200, densityFiles=NULL, boxPlots=TRUE, 
                         plotPer='markerSet', includeSortByMarker=TRUE, consistentYScaling=TRUE, 
                         filePer=NULL,pdfHeight=8.5,pdfWidth=11){

    ### create list of markers from all cell type files
    markerList <- list()
    for(ctf in markerSetFiles){
        msetName <- gsub(".txt","",basename(ctf))
        markerList[[msetName]] <- scan(ctf,"")
        if(!is.null(funcMarker)){
            markerList[[msetName]] <- c(markerList[[msetName]],paste0(markerList[[msetName]],",",funcMarker))
        }
    }
    ## collapse marker list into vector for calculation of area and density
    markers <- unique(unlist(markerList))
    ## add marker to sort by if it's not in the list
    if(!is.null(sortByMarker) & !sortByMarker %in% markers){
        markers <- c(markers,sortByMarker)
    }

    rho <- NULL
   
    if(!is.null(densityFiles)){
        for(df in densityFiles){
            den <- readRDS(df)
            if("density" %in% names(den)){
                rho <- bind_rows(rho, den$density)
            } else {
                rho <- bind_rows(rho, den)
            }
        }
        rho <- select(rho, c(Sample,SPOT,CellType,Total)) ## this shouldn't be needed after running again
    }
    if(!is.null(dataFiles) & !is.null(annotationsDirs)){
        for(n in 1:len(dataFiles)){
            dataFile <- dataFiles[n]
            dd <- readRDS(dataFile)
            print(paste0("getting IS stats for ",dataFile))
            flog.debug("getting IS stats for %s",dataFile)
            annotationsDir <- annotationsDirs[n] 
            print(paste0("using annotations in ",annotationsDir))
            flog.debug("using annotations in %s",annotationsDir)
            if(!byBand){
                ## for plotting TOTAL density around tumor interface, create one band/bin
                interfaceBins = c(-maxDistanceFromInterface,0,maxDistanceFromInterface)
            } else {
                mdfi <- maxDistanceFromInterface/10
                interfaceBins = (-mdfi:mdfi)*bandWidth
            }
            flog.debug("interface bins: ",paste(interfaceBins,collapse="  "))
            aFiles <- file.path(annotationsDir,dir(annotationsDir)[grep("\\.annotations$",
                                                                   dir(annotationsDir))])

            ## get density and area of interface stats
            flog.debug("calculating interface stats")
            samp_ia <- calculateInterfaceArea(aFiles, dd, pad, writeCSVfiles=writeCSVfiles,
                                         maxG=maxG,outDir=outDir, interfaceBins=interfaceBins,
                                         statsByBand=byBand)
            saveRDS(samp_ia, file=file.path(outDir,
                                              paste(getSampleFromFileName(dataFile),
                                                "area.rda",sep="_")))
            if(!is.null(samp_ia$area) & !is.null(samp_ia$bandAssignments)){
                ## calculate density
                density <- calculateDensity(dd, samp_ia$area, samp_ia$bandAssignments,
                                            markers, cellTypeName, pad, funcMarker=funcMarker,
                                            sortByMarker=sortByMarker,writeCSVfiles=writeCSVfiles,
                                            outDir=outDir, statsByBand=byBand)
                if(!is.null(density)){
                    if(byBand){
                        density <- density %>% select(-one_of(names(density)[grep("Count",names(density))]))
                    } else {
                        density <- density %>% select(Sample,SPOT,CellType,Total)
                    }
                    rho <- rho %>% bind_rows(density)
                    saveRDS(density, file=file.path(outDir, paste(getSampleFromFileName(dataFile),
                                                    "density.rda",sep="_")))
                }
            } else {
                if(!byBand){
                    ## calculate densities for entire FOV 
                    flog.debug("calculating stats for FOVs without tumor boundaries")
                    fs <- calculateFOVStats(aFiles, dd, markers, cellTypeName,
                                            pad, maxG, writeCSVfiles=writeCSVfiles,
                                            funcMarker=funcMarker,sortByMarker=sortByMarker,outDir)
                    if(!is.null(fs$density)){
                        rho <- fs$density %>% select(Sample,SPOT,CellType,Total) %>%
                                 bind_rows(rho)
                        saveRDS(fs, file=file.path(outDir,
                                              paste(getSampleFromFileName(dataFile),
                                                "FOV_area_density.rda",sep="_")))
                    }    
                }
            }
            flog.debug("saving density to %s",file.path(outDir,"all_samples_density.rda")) 
            ## save as you go
            saveRDS(rho,file=file.path(outDir,"all_samples_density.rda"))
        }
    }

    if(is.null(rho)){
        return()
    }

    ## remove any exclusions
    flog.debug("removing sample_fov/sample_marker exclusions")
    if(!is.null(exclude_sample_fov)){
         rho <- removeExclusions(rho, exclude_sample_fov=exclude_sample_fov)
    }

    #######
    ## BUILD PLOT(S)
    #######
    if(!byBand){
        flog.debug("generating total density plots")
        tdPlots <- plotTotalDensity(markerList,rho,funcMarker=funcMarker,funcPosColor=funcPosColor,
                         funcNegColor=funcNegColor,sampleColor=sampleColor,
                         sampleColorDark=sampleColorDark,spotOrder=spotOrder,
                         plotPer=plotPer,pdfFile=pdfFile,sampleOrder=sampleOrder,
                         sortByMarker=sortByMarker,includeSortByMarker=includeSortByMarker,
                         boxPlots=boxPlots,consistentYScaling=consistentYScaling)
    } else {
        mt <- rho %>%
                gather(key="Band",value="Density",names(rho)[2:(ncol(rho)-2)]) %>%
                spread(CellType,Density)
        dat <- rho %>% 
                gather(key="Band",value="Density",names(rho)[2:(ncol(rho)-2)])
        ## to do: change this to print inside plotDensityByBand instead of returning plot
        plotDensityByBand(dat, sampleName, c(markers,paste(markers,funcMarker,sep=",")))
    }


    ######
    ## PRINT PLOTS(S)
    ######
    if(!is.null(filePer)){
        flog.debug("printing one file per %s",filePer)
    } else {
        flog.debug("printing all plots to ONE file: %s",pdfFile)
    }
    printTotalDensityPlots(tdPlots, markerList, filePer=filePer, plotPer=plotPer, 
                           pdfFile=pdfFile, height=pdfHeight, width=pdfWidth)
}

