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
    tumorBoundariesMergeXY <- tumBTrim %>% bind_rows %>% dplyr::select(X,Y) %>% as.data.frame

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

    aa <- bind_cols(bandArea) %>% dplyr::select(matches("Area"))
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
  spots <- dd %>% dplyr::select(SPOT) %>% distinct(SPOT) %>% arrange(SPOT) %>% pull(SPOT)
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
        areas[[len(areas)+1]] <- bandArea %>% dplyr::select(Band,Area,Sample,SPOT=Spot)
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
    allBandAssignments <- allBandAssignments[!is.na(allBandAssignments$Band),]
  }

  if(writeCSVfiles){
    write_csv(as.data.frame(areass), file.path(outDir,cc(sampleName, "Interface_Area.csv")))
    write_csv(as.data.frame(allBandAssignments), file.path(outDir,cc(sampleName, "Band_Assignments.csv")))
  }

  return(list(area=areass, bandAssignments=allBandAssignments))
}




#' Calculate density of set of markers (marker combinations)
#'
#' Based on pre-computed areas and band assignments (may be one band, depending on
#' how area was calculated), calculate densities for a set of markers in each band
#' 
#' @param dat
#' @param areas
#' @param bandAssignments
#' @param markerSet
#' @param cellTypeName
#' @param pad
#' @param funcMarker
#' @param sortByMarker
#' @param writeCSVfiles
#' @param outDir
#' @param statsByBand
#' @return 
#' @export
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
      mtB <- mtCounts %>% filter(SPOT == spot, !is.na(Band))
      a <- areas %>% filter(SPOT == spot) %>%
                     gather(key="Band", value="Area", names(areas)[3:ncol(areas)])
      MM <- 1/a$Area

      densityByMarker <- as.tibble(lapply(mtB[,3:ncol(mtB)], function(x){ x * MM })) %>%
                         bind_cols(mtB[,c("SPOT","Band")]) %>%
                         dplyr::select(SPOT, Band, everything()) #%>%
      densityByMarker <- densityByMarker %>% gather(key="CellType", value="Total", 
                                                    names(densityByMarker)[3:ncol(densityByMarker)])
      mtB <- mtB %>% gather(key="CellType", value="Counts", names(mtB)[3:ncol(mtB)])
      densityWithCounts <- densityByMarker
      densityWithCounts$Counts <- mtB$Counts
      densityWithCounts$Sample <- sampleName

      rho[[len(rho)+1]] <- densityWithCounts

    } else {

      flog.debug("counting number of cells inside and outside the tumor for each marker")

      mtB <- mtCounts %>% filter(SPOT==spot,!is.na(Band))
      rn <- mtB$Band
      mtB <- as.matrix(mtB[3:ncol(mtB)],check.names=F)
      rownames(mtB) <- rn
      mtBCollapse <- apply(mtB,2,function(x){tapply(x,interfaceSide(rownames(mtB)),sum)})

      interfaceArea <- areas %>% filter(SPOT==spot)

      ## divide into 1, to later be multiplied by counts to calculate density
      inoutArea <- interfaceArea %>% dplyr::select(Inside,Outside)
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
    mtCounts <- mt %>% dplyr::select(markerSet) %>% colSums
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


getDensityPercentage <- function(den, by=c("Sample","SPOT","Band")){
    if("TotalDensity" %in% names(den)){
        den <- den %>% select(Density=TotalDensity, everything())
    }
    tmp <- den %>% group_by_(.dots = by) %>%
           summarise(Tot=sum(Density))
    den <- den %>% left_join(tmp) %>%
            mutate(Percent=Density/Tot) %>%
            select(-(Tot))
    return(den)
}


#' Get a table of marker combination densities
#' 
#' Generate a table of cell type densities
#'
#' @param dat                       tibble containing halo data, specifically Sample, UUID, SPOT, Marker,
#'                                  Xmin/max, Ymin/max, ValueType, Value
#' @param dataFiles                 vector of *.rda files containing Halo data for a single sample
#' @param annotationsDirs           vector of directories containing *.annotations files, XML files of boundary
#'                                  annotations from Halo; order must correspond to order of dataFiles
#' @param pad                       amount to trim from each FOV
#' @param funcMarker                character string of functional marker to plot on top of every other marker
#' @param sortByMarker              sort FOVs by density of this marker (cell type)in ascending order;
#'                                  default=NULL
#' @param writeCSVfiles             logical; write CSV files, one of density values and one of area values; 
#'                                  default=TRUE
#' @param maxG                      maximum factor of 10 to use for generating random points on grid [ REWORD??? ]
#' @param outDir                    if writeCSVfiles=T, directory to which these files will be written
#' @param byBand                    logical indicating whether to break down density by distance from tumor interface; 
#'                                  default=FALSE
#' @param bandWidth                 width of each band around tumor interface, to be used when byBand is TRUE; 
#'                                  default=10
#' @param maxDistanceFromInterface  include only those cells within this distance from tumor interface
#' @return tibble containing density for all markers in all samples
#' @export
getMarkerDensityTable <- function(markers, dat=NULL, outDir=getwd(), dataFiles=NULL, annotationsDirs=NULL, 
                                  areaFiles=NULL, densityFiles=NULL, sampleOrder=NULL, 
                                  writeCSVfiles=FALSE, pad=30, byBand=TRUE, maxDistanceFromInterface=200, 
                                  bandWidth=10, markerSetName=NULL, sortByMarker=NULL, funcMarker=NULL, 
                                  maxG=5, calcFOVarea=FALSE){

    
    interfaceBins <- (-(maxDistanceFromInterface/bandWidth):(maxDistanceFromInterface/bandWidth))*10

    ia  <- NULL  ## interface area
    ba  <- NULL  ## band assignments
    rho <- NULL  ## density

    if(!is.null(densityFiles)){
        print("reading density files")
        ## read density files
        for(df in densityFiles){
            rho <- rho %>% bind_rows(readRDS(df))
        }
    } else {

        if(is.null(dat)){
            print("reading data files")
            for(df in dataFiles){
                ddat <- readRDS(df)
                ddat <- ddat %>% select(Sample, UUID, SPOT, Marker, matches("XM|YM"), ValueType, Value)
                dat <- dat %>% bind_rows(ddat)
            }
        }
        dat$Sample <- gsub("_ObjectAnalysisData","",dat$Sample)

        if(!is.null(areaFiles)){
            print("reading area files")
            ## read area files
            for(af in areaFiles){
                a <- readRDS(af)
                ia <- ia %>% bind_rows(a$area)
                ba <- ba %>% bind_rows(a$bandAssignments)
            }
            ia$Sample <- gsub("_ObjectAnalysisData","",ia$Sample)
            ba$Sample <- gsub("_ObjectAnalysisData","",ba$Sample)
        } else {
            print("reading data and annotations files")
            ## calculate area
            if(is.null(dat) | is.null(annotationsDirs)){
                stop(paste0("Need either density files, area files + dat, or dat +",
                            " annotations directories"))
            }

            ### calculate area then add it to the area tibble for entire analysis
            for(annotationsDir in annotationsDirs){
                samp <- gsub("_.*","",basename(annotationsDir))
                print(paste0("calculating area for ",samp))
                dd <- dat %>% filter(Sample == samp)
                aFiles <- file.path(annotationsDir,dir(annotationsDir)[grep("\\.annotations$",
                                                                   dir(annotationsDir))])

                flog.debug("calculating interface area")
                samp_ia <- calculateInterfaceArea(aFiles, dd, pad, writeCSVfiles=writeCSVfiles,
                                                  maxG=maxG, outDir=outDir, statsByBand=byBand,
                                                  interfaceBins=interfaceBins)

                if(!is.null(samp_ia$area) & !is.null(samp_ia$bandAssignments)){
                    ia <- ia %>% bind_rows(samp_ia$area)
                    ba <- ba %>% bind_rows(samp_ia$bandAssignments)
                    saveRDS(samp_ia, file=file.path(outDir,
                                                    paste(samp,
                                                          "area.rda",sep="_")))
                } else {
                    if(!calcFOVarea){
                        print(paste0("WARNING: Interface area could not be calculated for sample ",samp))
                        flog.warn("Interface area could not be calculated for sample %s",samp)
                    } else if(byBand){
                        ## calculate densities for entire FOV 
                        flog.debug("calculating stats for FOVs without tumor boundaries")
                        fs <- calculateFOVStats(aFiles, dat, markers, cellTypeName,
                                        pad, maxG, writeCSVfiles=writeCSVfiles,
                                        funcMarker=funcMarker,sortByMarker=sortByMarker,outDir)
                        if(!is.null(fs$density)){
                            rho <- fs$density %>% 
                                   dplyr::select(Sample,SPOT,CellType,Total) %>%
                                   bind_rows(rho)
                            saveRDS(fs, file=file.path(outDir,
                                                       paste(getSampleFromFileName(dataFile),
                                                      "FOV_area_density.rda",sep="_")))
                        }
                    }
                }
            }
        }
        ## calculate density
        if(!is.null(ia) & !is.null(ba)){
            for(s in unique(ia$Sample)){
                print(paste0("calculating density for ",s))
                sdat <- dat %>% filter(Sample == s)
                sia <- ia %>% filter(Sample == s)
                sba <- ba %>% filter(Sample == s)
                den <- calculateDensity(sdat, sia, sba,
                                    markers, cellTypeName, pad, funcMarker=funcMarker,
                                    sortByMarker=sortByMarker,writeCSVfiles=writeCSVfiles,
                                    outDir=outDir, statsByBand=byBand)
                if(!is.null(den)){
                    if(byBand){
                        den <- den %>% dplyr::select(-one_of(names(den)[grep("Count",names(den))]))
                    } else {
                        den <- den %>% dplyr::select(Sample,SPOT,CellType,Total)
                    }
                    rho <- rho %>% bind_rows(den)
                    saveRDS(rho, file=file.path(outDir, "All_density.rda"))
                }
            }
        } else {
            print("WARNING: Interface area could not be calculated.")
            flog.warn("Interface area could not be calculated.")
        }
    }

    if(!"Density" %in% names(rho)){
        rho <- rho %>% dplyr::select(Density=Total, everything())
    }

    return(rho)
}

getInfiltrationDensitySummary <- function(markers, ia=NULL, ba=NULL, areaFiles=NULL){

    if(is.null(ia) | is.null(ba)){
        ## read areaFiles
        for(af in areaFiles){
            a <- readRDS(af)
            ia <- ia %>% bind_rows(a$area)
            ba <- ba %>% bind_rows(a$bandAssignments)
        }
    }
    ia$Sample <- gsub("_ObjectAnalysisData","",ia$Sample)
    ba$Sample <- gsub("_ObjectAnalysisData","",ba$Sample)
    bandAreas <- ia %>% gather(3:ncol(ia), key="Band", value="Area")

    markerCounts <-computeMultiMarkerTable(ba,markers) %>%
                   select(Sample,SPOT,Band,markers) %>%
                   filter(!is.na(Band))

    markerCounts <- markerCounts %>%
                    gather(4:ncol(markerCounts), key="CellType", value="Positive") %>%
                    group_by(Sample,SPOT,Band,CellType) %>%
                    summarise(Count=sum(Positive))

    summaryDat <- left_join(markerCounts, bandAreas) %>%
               group_by(Sample,Band,CellType) %>%
               summarise(TotalCountAllFOV=sum(Count), TotalAreaAllFOV=sum(Area)) %>%
               mutate(TotalDensity=TotalCountAllFOV/TotalAreaAllFOV)

    return(summaryDat)
}
