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
computeBandAreasMC <- function(tumorBoundariesList,bbData,interfaceBins,nGrid=100,maxSeg=1000,
    exclusionB) {

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
#' @param aFile          Halo boundary annotation file in XML format
#' @param maxG           maximum factor of 10 to use for generating random points on grid [ REWORD??? ]
#' @param bbData         a list containing X0,X1,Y0,Y1 representing the boundary
#'                       of the trimmed FOV
#' @param interfaceBins  vector of distances (integers) that define each band; default = (-20:20)*10
#' @return a dataframe with a row for each interface bin and a column for each area calculation (one for 
#'         10^maxG and one for 10^maxG-1)
calculateBandAreas <- function(aFile,maxG=5,bbData,interfaceBins) {

    if(is.null(interfaceBins)){
        interfaceBins <- (-20:20)*10
    }

    flog.debug("getting halo boundaries")
    boundaries <- readHaloAnnotations(aFile)

    flog.debug("removing contained boundaries")
    ## remove boundaries that are completely contained in another one
    sepB <- removeContainedBoundaries(boundaries)
    tumB <- sepB$tumB
    excB <- sepB$excB

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
#' @param aFile                Halo boundary annotation file in XML format
#' @param sample               sample name
#' @param fov                  FOV
#' @param maxG                 [TO DO: FIND OUT WHAT THIS IS]
#' @param haloInfiltrationDir  directory containing halo infiltration files with suffix 
#'                             "_infiltration_ex_data"
#' @param writeCSVfiles        write lattice band area file in csv format; default=FALSE
#' @param outDir               if writeCSVfiles=TRUE, files will be written here
#' @param bbData               a list containing X0,X1,Y0,Y1 representing the boundary
#'                             of the trimmed FOV
#' @param interfaceBins        vector of distances (integers) that define each band; default = (-20:20)*10 
#' @return tibble 
computeBoundaryBinsLattice <- function(aFile, sample, fov, maxG, haloInfiltrationDir, outDir=NULL,
                                       bbData=NULL, writeCSVfiles=FALSE, interfaceBins=NULL){ 

    if(is.null(interfaceBins)){
        interfaceBins <- (-20:20)*10
    }

    SAMPLE <- sample
    sampleSPOT <- as.numeric(fov)
    flog.debug("%s SPOT %s",SAMPLE,sampleSPOT)
    aa=calculateBandAreas(aFile,maxG=maxG,bbData,interfaceBins)

    flog.debug("reading infiltration data file %s", paste0(SAMPLE,sampleSPOT,"_infiltration_ex_data.csv"))
    mp <- read_csv(file.path(haloInfiltrationDir,paste0(SAMPLE,sampleSPOT,"_infiltration_ex_data.csv")))
    imp <- grep("Band",colnames(mp))
    a2 <- data.frame(aa[,ncol(aa),drop=F],AreaMP=as.numeric(mp[,imp]))

    rmsConvergence1 <- NA
    if(ncol(aa)>2){
        rmsConvergence1 <- sqrt(mean((aa[,ncol(aa)-1]-aa[,ncol(aa)-2])^2))
    }
    rmsConvergence <- sqrt(mean((aa[,ncol(aa)]-aa[,ncol(aa)-1])^2))
    rmsHalo <- sqrt(mean((a2[,ncol(a2)]-a2[,ncol(a2)-1])^2))

    amin <- min(a2)
    amax <- max(a2)
    amax <- max(amax,2.5*amin)

    flog.debug("creating tibble")
    xx <- as.tibble(data.frame(
            Sample=SAMPLE,
            Spot=sampleSPOT,
            Version="v3.4",
            MaxG=maxG,
            Band=rownames(a2),
            Area=a2[,1],
            Halo=a2[,2])
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



joinBandArea<-function(bdat,sampleName,spot, aFile, maxG, haloInfiltrationDir, outDir, bbData) {

    bdat$Band=as.character(bdat$Band)
    flog.debug("computing boundary bins lattice")
    bandArea <- computeBoundaryBinsLattice(aFile, sampleName, spot, maxG=maxG, haloInfiltrationDir, outDir=outDir, bbData=bbData)
    flog.debug("joining bdat and bandArea")
    xx=full_join(bdat,bandArea)
    flog.debug("filling in NAs with 0s")
    xx$n[is.na(xx$n)]=0
    xx$Band=factor(xx$Band,levels=bandArea$Band)
    flog.debug("arranging by Band")
    xx %<>% arrange(Band)
    xx

}
      

#' Calculate total area and density of interface 
#' 
#' Calculate area and density of cells within a given distance to a tumor boundary
#'
#' @param aFiles              a vector of Halo boundary annotation files, with 
#'                            '.annotations' extension in XML format
#' @param dataFile            a *.rda file with Halo data for one sample
#' @param cellTypeFile        file containing cell types to be calculated; each line
#'                            is one cell type and each cell type can be either a single
#'                            marker name or a comma-separated list of markers
#' @param cellTypeName        a name to represent the cell types in marker file
#' @param funcMarker          this marker will be added to all markers in cellTypeFile in order to
#'                            compare cells that are +/- for this functional marker [ TO DO: REWORD THIS ]
#' @param fovBoundingBox      a list with four elements: X0, X1, Y0, Y1, representing the minimum
#'                            and maximum coordinates of the FOV to be plotted
#' @param plotBoundingBox     a list with four elements: X0, X1, Y0, Y1, representing the minimum
#'                            and maximum coordinates of the outermost plot box [ TO DO: UPDATE THIS ] 
#' @param pad                 amount to trim from FOV before calculating
#' @param writeCSVfiles       logical indicating whether to write density and area tables to file; 
#'                            default=TRUE
#' @param haloInfiltrationDir directory containing Halo Infiltration data (csv, plot files)
#' @param maxG                the maximum factor of 10 to use for generating random points for area calculation
#' @param outDir              if writeCSVfiles=T, directory to which these files will be written
#' @return a list containing two tables: 'density' and 'area' 
#' @export
calculateInterfaceStats <- function(aFiles, dataFile, cellTypeFile, cellTypeName,
                                    fovBoundingBox, pad, plotBoundingBox, funcMarker=NULL,
                                    writeCSVfiles=TRUE,haloInfiltrationDir,maxG,outDir){
  pad <- as.numeric(pad)
  bbFOV <- fovBoundingBox
  bbPlot <- plotBoundingBox

  flog.debug("reading data file %s",dataFile)
  dd <- readRDS(dataFile)
  sampleName <- getSampleFromFileName(dataFile)

  fovTag <- "ORIG"
  bbData <- bbFOV
  if(pad > 0){
      fovTag <- paste0("clip",pad)
      bbData <- padBoundingBox(bbFOV,-pad/pixel2um)
  }

  #interfaceBins <- (-20:20)*10

  rho <- list()
  areas <- list()

  flog.debug("getting spots")
  spots <- dd %>% dplyr::select(SPOT) %>% distinct(SPOT) %>% arrange(SPOT) %$% as.vector(SPOT)
  for(spot in spots){
    flog.debug("working on spot %s",spot)
    sTag <- cc(sampleName, paste0("Spot",spot,".annotations"))
    aai <- grep(sTag, aFiles)
    if(len(aai) > 0){
      flog.debug("filtering dataset for spot..")
      ds <- dd %>% filter(SPOT==spot & ValueType=="Positive")
      aFile <- aFiles[aai]

      ## get Distance and Band assignments for all points that fall inside a tumor
      flog.debug("getting band counts")
      ds <- getBandCounts(ds, bbData,bbPlot,bbFOV, aFile)
      if(is.null(ds)){ 
        next
      }
      
      flog.debug("filtering band count data for CD3+ and counting cells in each band")
      bdat <- ds %>% filter(CD3==1) %>% filter(!is.na(Band)) %>% count(Band) %$% data.frame(Band,n)
      
      flog.debug("joining band area")
      bdat <- joinBandArea(bdat, sampleName, spot, aFile, maxG, haloInfiltrationDir, outDir, bbData)
      
      markerSet <- scan(cellTypeFile,"")
      if(!is.null(funcMarker)){
        markerSet <- c(markerSet,paste0(markerSet,",",funcMarker))
      }
      fillNA <- list()
      dummy <- lapply(markerSet, function(x){ fillNA[x] <<- 0 })
      flog.debug("computing multimarker table and summarizing each band")
      mt <- computeMultiMarkerTable(ds,markerSet)
      mtCounts <- mt %>%
                  group_by(Band) %>%
                  summarize_at(markerSet,sum) %>%
                  complete(Band,fill=fillNA)
      mtB <- as.matrix(mtCounts[!is.na(mtCounts$Band),2:ncol(mtCounts)],check.names=F)
      rownames(mtB) <- mtCounts$Band[!is.na(mtCounts$Band)]

      flog.debug("counting number of cells inside and outside the tumor for each marker")
      mtBCollapse <- apply(mtB,2,function(x){tapply(x,interfaceSide(rownames(mtB)),sum)})

      flog.debug("calculating total area for bands each inside and outside tumor")
      bandAreaCollapse <- tapply(bdat$Area,interfaceSide(bdat$Band),sum)
      MM <- diag(1/bandAreaCollapse)
      rownames(MM) <- c("Inside","Outside")
      colnames(MM) <- c("Inside","Outside")

      flog.debug("getting total counts and area for ALL bands")
      totalCounts <- colSums(mtBCollapse)
      totalArea <- sum(bandAreaCollapse)

      flog.debug("calculating density")
      rho[[len(rho)+1]] <- (t(mtBCollapse) %*% MM) %>%
                            as.data.frame %>%
                            cbind(data.frame(Counts=t(mtBCollapse))) %>%
                            rownames_to_column("CellType") %>%
                            mutate(Total=totalCounts/totalArea,Sample=sampleName,SPOT=spot)

      flog.debug("getting area")
      areas[[len(areas)+1]] <- data.frame(Area=bandAreaCollapse) %>%
                                rownames_to_column("Side") %>%
                                as_tibble %>%
                                mutate(Sample=sampleName,SPOT=spot)
    }
  }
  rhos <- NULL
  areass <- NULL
  flog.debug("Num elements in list rho: %s",length(rho))
  if(length(rho) > 0){
    rhos <- rho %>% bind_rows
  }
  flog.debug("Num elements in list areas: %s",length(areas))
  if(length(areas) > 0){
    areass <- areas %>% bind_rows %>% spread(Side,Area)
  }

  if(writeCSVfiles){
    ctf <- basename(gsub("\\.txt","",cellTypeFile))
    write_csv(as.data.frame(areass), cc(sampleName, ctf, "interfaceStatsV3_Area.csv"))
    write_csv(as.data.frame(rhos), cc(sampleName, ctf, "interfaceStatsV3_Density.csv"))
  }

  return(list(density=rhos, area=areass))
}

#' Remove exclusion boundaries that are contained in another one
#' 
#' Remove data from boundaries table that represent exclusion boundaries
#' that are completely surrounded by another exclusion boundary
#' 
#' @param boundaries            table generated by readHaloAnnotations
#'                             
removeContainedBoundaries <- function(boundaries){

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
#' @param aFile  Halo boundaries annotation file in XML format
removeExcludedPoints<-function(ds,aFile) {

    ds=ds %>%
        mutate(X=(XMax+XMin)/2,Y=-(YMax+YMin)/2) %>%
        dplyr::select(Sample,SPOT,UUID,X,Y,Marker,Value)

    boundaries <- readHaloAnnotations(aFile)
    regionTable <- boundaries %>% bind_rows %>% count(RegionNum,RegionCode)
    excB <- boundaries[regionTable$RegionCode!="Tum"]

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
#' @param bbPlot          a list containing X0,X1,Y0,Y1 to be used as boundary of
#'                        entire plot
#' @param bbFOV0          a list containing X0,X1,Y0,Y1 showing the actual boundary
#'                        of the FOV (untrimmed)
#' @param aFile           Halo boundaries annotations file in XML format 
#' @param interfaceBins   vector of distances (integers) that define each band; default = (-20:20)*10 
getBandCounts<-function(ds,bbData,bbPlot,bbFOV0,aFile,interfaceBins=NULL){

    flog.debug("getting boundaries")
    ## boundary XML data to table
    boundaries <- readHaloAnnotations(aFile)
    sepB <- removeContainedBoundaries(boundaries)
    tumB <- sepB$tumB
    excB <- sepB$excB
    if(is.null(tumB) || len(tumB) == 0){
        #flog.info("No TUMOR boundaries found in file %s. Skipping.")
        #return()
    }

    if(is.null(interfaceBins)) {
        interfaceBins <- (-20:20)*10
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

#' Plot total density of each FOV for one sample
#'
#' @param dataFile               *.rda file containing Halo data for a single sample
#' @param annotationsDir         directory containing *.annotations files, XML files of boundary
#'                               annotations from Halo
#' @param cellTypeFile           a text file with each line containing a cell type marker; may be
#'                               marker combinations delimited by commas
#' @param cellTypeName           a character string naming the group of cell types (e.g., T Cells)
#' @param fovBB                  a list of four elements: X0,X1,Y0,Y1, representing the min/max coords
#'                               for FOVs
#' @param pad                    amount to trim from each FOV
#' @param plotBB                 a list of four elements: X0,X1,Y0,Y1, representing the min/max coords
#'                               for whole plots
#' @param pdfFile                output PDF file
#' @param logPlot                logical; when set to TRUE, plot density per mm^2; default=FALSE
#' @param funcMarker             character string of functional marker to plot on top of every other marker
#' @param sampleColor            character string indicating HTML color to represent this particular sample
#' @param sampleColorDark        a slightly darker shade of sampleColor, to be used for density of funcMarker+
#' @param sortByMarker           sort FOVs by density of this marker (cell type)in ascending order;
#'                               default=CD3
#' @param exclude_sample_fov     list where names are sample names and values are vectors of FOV to be  
#'                               removed from data before any calculations 
#' @param exclude_sample_marker  list where names are sample names and values are vectors of markers to be 
#'                               should be removed from data before any calculations
#' @param writeCSVfiles          logical; write CSV files, one of density values and one of area values; 
#'                               default=TRUE
#' @param ymax                   maximum y value; by default, this value is dictated by the data
#' @param haloInfiltrationDir    directory containing Halo Infiltration data (csv, plot files)
#' @param maxG                   ???
#' @param outDir                 if writeCSVfiles=T, directory to which these files will be written#' @return nothing
#' @export
plotTotalDensity <- function(dataFile, annotationsDir, cellTypeFile, cellTypeName,
                             fovBB, pad, plotBB, pdfFile, ymax=NULL,
                             logPlot=F, funcMarker="KI67", sortByMarker="CD3", sampleColor="#f16913",
                             sampleColorDark="#80380A", exclude_sample_fov=NULL, exclude_sample_marker=NULL,
                             writeCSVfiles=TRUE, haloInfiltrationDir,maxG,outDir){

    aFiles <- file.path(annotationsDir,dir(annotationsDir)[grep("\\.annotations$",dir(annotationsDir))])

    ## get density and area
    flog.debug("calculating interface stats")
    is <- calculateInterfaceStats(aFiles, dataFile, cellTypeFile, cellTypeName,
                                 fovBB, pad, plotBB, writeCSVfiles=writeCSVfiles,
                                 funcMarker=funcMarker,haloInfiltrationDir,maxG,outDir)

    rho <- is$density
    if(is.null(rho)){
        return()
    }

    ## remove any exclusions
    flog.debug("removing exclusions")
    if(!is.null(exclude_sample_fov) || !is.null(exclude_sample_marker)){
        rho <- removeExclusions(rho, exclude_sample_fov=exclude_sample_fov,
                             exclude_sample_marker=exclude_sample_marker)
    }

    ## get max density
    mt     <- rho %>% dplyr::select(CellType,Sample,SPOT,Total) %>% spread(CellType,Total)
    if(is.null(ymax)){
        maxRho <- mt %>% dplyr::select(-Sample,-SPOT) %>% summarize_all(funs(max)) %>% max
    } else {
        maxRho <- ymax
    }
    flog.debug("maxRho = %s",maxRho)
    sampleName <- mt$Sample[1]
    flog.debug("sorting by %s",sortByMarker)
    sortedByMarker <- cc(sampleName,mt %>% arrange(get(sortByMarker)) %>% pull(SPOT))
    markers    <- scan(cellTypeFile, "")

    pdf(file=pdfFile, width=11,height=8.5)

    ddt <- NULL
    for(mi in markers){
        miF <- NULL

        selectCols <- c("SampleID",mi)
        if(!is.null(funcMarker)){
            flog.debug("adding %s to marker list",funcMarker)
            miF <- paste0(mi,",",funcMarker)
            flog.debug(paste0("selecting cols: ",paste(selectCols,collapse=", ")))
            selectCols <- c(selectCols, miF)
        }

        dt <- mt %>%
              unite(SampleID,Sample,SPOT) %>%
              dplyr::select(selectCols) %>%
              remove_rownames %>%
              as.data.frame %>%
              column_to_rownames("SampleID")

        ## transpose so markers are rows
        dt=as.matrix(t(dt[sortedByMarker,,drop=F]))

        ## if functional marker is given, calculate values of markers "neg" for
        ## funcMarker by subtracting the value of markers "pos" from it (??)
        if(nrow(dt) > 1){
            ddt <- rbind(dt[2,],dt[1,]-dt[2,])
            rownames(ddt) <- c(paste0(funcMarker,"+"),paste0(funcMarker,"-"))
            ddt <- cbind(ddt,Median=apply(ddt,1,median) )
        } else {
            ddt <- dt
        }

        eb <- mad(apply(dt,2,sum))
        cols <- c(rep("#dddddd",len(ddt)-1),sampleColor)
        nP <- ncol(ddt)

        if(!logPlot) {
            bb=barplot(ddt,ylim=c(0,maxRho),xaxt='n',beside=F,main=mi,ylab="Density (per mm^2)")

            mdt=ddt
            mdt[,-ncol(mdt)]=0
            barplot(mdt,ylim=c(0,maxRho),xaxt='n',beside=F,add=T,
                col=c(sampleColorDark,sampleColor))

            text(bb+.5,-par()$usr[4]/40,colnames(ddt),xpd=T,srt=45,pos=2)
            dMed=sum(ddt[,nP])
            arrows(bb[nP],dMed-eb,bb[nP],dMed+eb,length=.1,angle=90,code=3)
        } else {
            bb=barplot(ddt+1,ylim=c(1,1e4),las=2,beside=F,main=mi,
                    ylab="Density (per mm^2)",log='y')
            dMed=sum(ddt[,nP])
            arrows(bb[nP],dMed-eb,bb[nP],dMed+eb,length=.1,angle=90,code=3)
        }
    }
    if(!is.null(funcMarker)){
        ## draw key
        td=ddt
        td[T]=0
        bb=barplot(td,ylim=c(0,maxRho),xaxt='n',yaxt='n',beside=F,main="",
                      ylab="",legend=T,args.legend=list(x=10,y=maxRho/2,cex=2))
    }

    dev.off()
}

