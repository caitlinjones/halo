#' Generate data frame of random points (is it random???)
#' 
#' Generate a data frame of random points that fall within the plot area
#' 
#' @param   nGrid   number of points to generate
#' @param   bbG     a list containing X0,X1,Y0,Y1 representing the boundary
#'                  of the plot area
#' @return data frame of X and Y values of the random points
generateSobolGrid <- function(nGrid,bbG) {
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
#' @param bb                a list containing X0,X1,Y0,Y1 representing the boundary
#'                              of the trimmed FOV
#' @param interfaceBins         vector of distances (integers) that define each band; default = (-20:20)*10 
#' @param nGrid                 number of random points to use for area calculations; default=100
#' @param maxSeg                maximum segments to divide long boundary segments into; default=1000
#' @param exclusionsB           a list where each element is a tibble containing exclusion boundaries
#'                              returned from readHaloAnnotations()
computeBandAreasMC <- function(tumorBoundariesList,bb,interfaceBins,nGrid=100,
                               maxSeg=1000,exclusionB) {

    ## generate random points all across the plot
    flog.debug("generating sobol grid")
    gg <- generateSobolGrid(nGrid,bb)
    if(len(exclusionB)>0) {
        flog.debug("There are %s exclusion boundaries. Removing excluded points",len(exclusionB))
        excludedPts <- pointsInsidePolygonList(gg,exclusionB)
        gg <- gg[!excludedPts,]
    }

    flog.debug("trimming boundaries")
    tumBTrim <- tumorBoundariesList %>%
        map(trimDFToBB,bb) %>%
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
    bandAreas <- as.data.frame.table(p2tomm*areaBB(bb)*areaTable/nGrid)
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
#' @param bb         a list containing X0,X1,Y0,Y1 representing the boundary
#'                       of the trimmed FOV
#' @param interfaceBins  vector of distances (integers) that define each band; default = (-20:20)*10
#' @return a dataframe with a row for each interface bin and a column for each area calculation (one for 
#'         10^maxG and one for 10^maxG-1)
calculateBandAreas <- function(allBoundaries,maxG=5,bb,interfaceBins) {

    if(is.null(interfaceBins)){
        interfaceBins <- (-20:20)*10
    }
    tumB <- allBoundaries$tumB
    excB <- c(allBoundaries$excB,allBoundaries$epiB,allBoundaries$glsB)

    bandArea <- list()
    for(nGrid in 10^((maxG-1):maxG)) {
        maxSeg <- max(unlist(bb))/sqrt(10*nGrid)
        flog.debug("nGrid = %s  maxSeg = %s ",nGrid,maxSeg)
        bandArea[[as.character(nGrid)]] <- computeBandAreasMC(tumB,bb,interfaceBins,
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
#' @param bb               a list containing X0,X1,Y0,Y1 representing the boundary
#'                             of the trimmed FOV
#' @param interfaceBins        vector of distances (integers) that define each band; default = (-20:20)*10 
#' @return tibble 
getAreaPerBand <- function(allBoundaries, sample, fov, maxG, outDir=NULL,
                                       bb=NULL, writeCSVfiles=FALSE, interfaceBins=NULL){

    if(is.null(interfaceBins)){
        interfaceBins <- (-20:20)*10
    }

    SAMPLE <- sample
    sampleSPOT <- as.numeric(fov)
    flog.debug("%s SPOT %s",SAMPLE,sampleSPOT)
    aa=calculateBandAreas(allBoundaries,maxG=maxG,bb,interfaceBins)

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
        if(!is.null(bb)){
            maxSeg=max(unlist(bb))/sqrt(10^(maxG+1))
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
joinBandArea <- function(bdat,bandArea) {

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
removeExcludedPoints <- function(ds,aFile=NULL,excB=NULL) {

    if(!all(c('X','Y') %in% colnames(ds))){
        ds <- convertCellMinMaxToMidpoints(ds)
    }

    if(is.null(excB)){
        boundaries <- readHaloAnnotations(aFile)
        allBoundaries <- cleanBoundaries(boundaries)
        excB <- c(allBoundaries$excB,allBoundaries$glsB,allBoundaries$epiB)
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
#' @param  ds             a tibble containing data for one FOV
#' @param bb              a list containing X0,X1,Y0,Y1 representing the boundary
#'                        of the trimmed FOV
#' @param allBoundaries   a list returned from readHaloAnnotations
#' @param interfaceBins   vector of distances (integers) that define each band; default = (-20:20)*10 
getBandAssignments <- function(ds,bb,allBoundaries,interfaceBins=(-36:36)*10){

    flog.debug("removing cells in exclusion boundaries")
    ## in case excluded cells have not been removed
    if("EXCLUDE" %in% names(ds)){
        ds <- ds %>% filter(EXCLUDE == "")
    }

    flog.debug("getting boundaries")
    tumB <- allBoundaries$tumB
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
    ds <- trimDFToBB(ds,bb)

    ## add values for negative markers
    for(mi in markers) {
      ds[[paste0(mi,"-")]]=ifelse(ds[[mi]]==0,1,0)
    }

    ## create a SpatialPointsDataFrame
    spCells <- ds %>% dplyr::select(UUID,X,Y)
    coordinates(spCells) <- ~X+Y

    flog.debug("getting tumor contour")
    ## filter data for cells that are near a tumor boundary AND inside
    ## the bounding box 
    tumorContour <- tumB %>% map(trimDFToBB,bb)
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
calculateInterfaceArea <- function(dat, aFiles=NULL, haloAnnotations=NULL, writeCSVfiles=TRUE,
                                   maxG=5,outDir=NULL, interfaceBins=(-36:36)*10){

    areas <- list()
    #bandAssignments <- list()

    for(sampleName in unique(dat$Sample)){
        flog.debug(sampleName)
        dd <- dat %>% filter(Sample == sampleName)
        flog.debug("getting spots")
        spots <- dd %>% dplyr::select(SPOT) %>% distinct(SPOT) %>% arrange(SPOT) %>% pull(SPOT)
        for(spot in spots){
            flog.debug(spot)
            allBoundaries <- NULL
            flog.debug("working on spot %s",spot)

            flog.debug("filtering dataset for spot..")
            ds <- dd %>% filter(SPOT==spot)

            if(is.null(haloAnnotations) && !is.null(aFiles)){
                sTag <- cc(sampleName, paste0("Spot",spot,".annotations"))
                aai <- grep(sTag, aFiles)
                if(len(aai) > 0){
                    aFile <- aFiles[aai]
                    flog.debug("getting halo boundaries")
                    allBoundaries <- cleanBoundaries(readHaloAnnotations(aFile)) 
                }
            } else if(!is.null(haloAnnotations)){
                allBoundaries <- haloAnnotations[[sampleName]][[as.character(spot)]]
            }
            ## skip if no tumors
            if(is.null(allBoundaries) || length(allBoundaries) == 0 || len(allBoundaries$tumB) == 0){
                flog.debug("no tumor boundaries found. skipping")
                next
            }
            allBoundaries$tumB <- allBoundaries$tumB[!is.na(names(allBoundaries$tumB))] ### TO DO: invesigate why
                                                                                        ### NAs are showing up (PR,23)
            ##
            ## get Distance and Band assignments for all points that fall inside a tumor
            ##
            #flog.debug("adding to data distance from tumor boundary and band assignment")
            #ba <- getBandAssignments(ds, bb, allBoundaries, interfaceBins)
            #bandAssignments[[len(bandAssignments)+1]] <- ba

            ## remove cells that do not fall within a band
            #bdat <- ba %>% filter(!is.na(Band)) %>% count(Band) %$% data.frame(Band,n)

            flog.debug("calculating area of each band")
            bandArea <- getAreaPerBand(allBoundaries,sampleName,spot,maxG=maxG,
                                       outDir=outDir, bb=bb,
                                       interfaceBins=interfaceBins)

            ## assume calcuating by band (so that when it comes time to do total, in theory 
            ## we can just use one large band
            areas[[len(areas)+1]] <- bandArea %>% dplyr::select(Band,Area,Sample,SPOT=Spot)
        }
    }
    areass <- NULL
    #allBandAssignments <- NULL
    flog.debug("Num elements in list 'areas': %s",length(areas))
    if(length(areas) > 0){
        areass <- areas %>% bind_rows %>% spread(Band,Area)
    }
    #if(length(bandAssignments)>0){
    #    allBandAssignments <- bandAssignments %>% bind_rows
    #    allBandAssignments <- allBandAssignments[!is.na(allBandAssignments$Band),]
    #}

    if(writeCSVfiles){
        write_csv(as.data.frame(areass), file.path(outDir, "all_samples_interface_area.csv"))
        #write_csv(as.data.frame(allBandAssignments), file.path(outDir, "all_samples_band_assignments.csv"))
    }

    return(list(area=areass, bandAssignments=allBandAssignments))
}

#' Calculate distances from each cell to nearest point on tumor boundary 
#' 
#' Calculate distances from each cell to nearest point on tumor boundary
#'
#' @param aFiles              a vector of Halo boundary annotation files, with 
#'                            '.annotations' extension in XML format
#' @param dat                 Halo data for one sample or more samples, loaded from *.rda file(s)  
#' @param writeCSVfiles       logical indicating whether to write density and area tables to file; 
#'                            default=TRUE
#' @param outDir              if writeCSVfiles=T, directory to which these files will be written
#' @param interfaceBins       a vector of distances from tumor interface in which cells will be binned; 
#'                            default=(-36:36)*10
#' @return table of all cells including distance to nearest tumor boundary and band assignments 
#' @export
calcDistancesFromTumorBoundary <- function(dat, aFiles=NULL, haloAnnotations=NULL, 
                                           writeCSVfiles=TRUE, outDir=NULL, interfaceBins=(-36:36)*10){

    bandAssignments <- list()

    for(sampleName in unique(dat$Sample)){
        flog.debug(sampleName)
        dd <- dat %>% filter(Sample == sampleName)
        flog.debug("getting spots")
        spots <- dd %>% dplyr::select(SPOT) %>% distinct(SPOT) %>% arrange(SPOT) %>% pull(SPOT)
        for(spot in spots){
            flog.debug(spot)
            allBoundaries <- NULL
            flog.debug("working on spot %s",spot)

            flog.debug("filtering dataset for spot..")
            ds <- dd %>% filter(SPOT==spot)

            if(is.null(haloAnnotations) && !is.null(aFiles)){
                sTag <- cc(sampleName, paste0("Spot",spot,".annotations"))
                aai <- grep(sTag, aFiles)
                if(len(aai) > 0){
                    aFile <- aFiles[aai]
                    flog.debug("getting halo boundaries")
                    allBoundaries <- cleanBoundaries(readHaloAnnotations(aFile))
                }
            } else if(!is.null(haloAnnotations)){
                allBoundaries <- haloAnnotations[[sampleName]][[as.character(spot)]]
            }

            ## skip if no tumors
            if(is.null(allBoundaries) || length(allBoundaries) == 0 || len(allBoundaries$tumB) == 0){
                flog.debug("no tumor boundaries found. skipping")
                next
            }
            allBoundaries$tumB <- allBoundaries$tumB[!is.na(names(allBoundaries$tumB))] 

            ##
            ## get Distance and Band assignments for all points that fall inside a tumor
            ##
            flog.debug("adding to data distance from tumor boundary and band assignment")
            ba <- getBandAssignments(ds, bb, allBoundaries, interfaceBins)
            bandAssignments[[len(bandAssignments)+1]] <- ba
        }
    }

    allBandAssignments <- NULL

    if(length(bandAssignments)>0){
        allBandAssignments <- bandAssignments %>% bind_rows
        allBandAssignments <- allBandAssignments[!is.na(allBandAssignments$Band),]
    }

    if(writeCSVfiles){
        write_csv(as.data.frame(allBandAssignments), file.path(outDir, "all_samples_band_assignments.csv"))
    }

    return(allBandAssignments)
}



#' Calculate density of set of markers (marker combinations)
#'
#' Based on pre-computed areas and band assignments (may be one band, depending on
#' how area was calculated), calculate densities for a set of markers in each band
#' 
#' @param areas              tibble containing area by band
#' @param bandAssignments    tibble containing band assignments for each cell
#' @param markerSet          vector of markers for which density should be calculated
#' @param funcMarker         functional marker for which to calculate density; default=NULL
#' @param sortByMarker       sort data by density of this marker; default=NULL
#' @param writeCSVfiles      logical; write density to csv files; default=FALSE
#' @param outDir             logical; if writeCSVfiles is TRUE, write them to this directory
#' @param statsByBand        logical; calculate density by interval bands
#' @return tibble containing densities for each marker type
#' @export
OLD_calculateInfiltrationDensity <- function(areas, bandAssignments, markerSet, 
                              funcMarker=NULL, sortByMarker=NULL,writeCSVfiles=TRUE,
                              outDir=NULL,statsByBand=FALSE){

  sampleName <- gsub("_ObjectAnalysisData","",unique(bandAssignments$Sample))

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
      bnds <- unique(mtB$Band)
      a <- areas %>% filter(SPOT == spot) %>%
                     select(Sample,SPOT,which(names(areas) %in% bnds)) #%>%
                     #gather(3:length(bnds)+2, key="Band", value="Area")
      a <- a %>% gather(3:ncol(a), key="Band", value="Area")
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


#' Calculate density of set of markers (marker combinations)
#'
#' Based on pre-computed areas and band assignments (may be one band, depending on
#' how area was calculated), calculate densities for a set of markers in each band
#' 
#' @param areas              tibble containing area by band
#' @param bandAssignments    tibble containing band assignments for each cell
#' @param markerSet          vector of markers for which density should be calculated
#' @param writeCSVfiles      logical; write density to csv files; default=FALSE
#' @param outDir             logical; if writeCSVfiles is TRUE, write them to this directory
#' @return tibble containing densities for each marker type
#' @export
calculateInfiltrationDensity <- function(areas, bandAssignments, markerSet, outFile, writeCSVfiles=TRUE){

    fillNA <- list()
    dummy <- lapply(markerSet, function(x){ fillNA[x] <<- 0 })

    flog.debug("computing multimarker table and summarizing each band")
    mt <- computeMultiMarkerTable(bandAssignments,markerSet)

    bandidx <- which(names(mt) == "Band")
    rhos <- mt %>% select(Sample,SPOT,Band,bandidx:ncol(mt))
    
    rhos <- rhos %>% 
           filter(!is.na(Band)) %>% 
           gather(4:ncol(rhos), key="CellType", value="Positive") %>% 
           group_by(Sample, SPOT, Band, CellType) %>% 
           summarise(Counts=sum(Positive)) %>%
           left_join(areas, by=c("Sample","SPOT","Band")) %>%
           mutate(Density=Counts/BandArea)

    if(writeCSVfiles){
        write_csv(as.data.frame(rhos), outFile)
    }

    return(rhos)
}

#' Get all absolute and percentage density values for complete FOVs
#' 
#' Given FOV-level density values, add sample-level values to data frame by  
#' calculating total densities of each marker in ALL FOV for a single sample
#'
#' @param den
#' @param areaDat
#' @return table of all density values
#' @export
getAllFOVDensityValues <- function(den, areaDat){

    totalSpotCounts <- den %>% 
                       group_by(Sample,SPOT) %>%
                       summarise(TotalSpotCounts=sum(Counts)) %>%
                       ungroup()

    #allVals <- den %>% 
    #           left_join(areaDat, by=c("Sample","SPOT")) 

    allVals <- den %>%
               left_join(areaDat, by=intersect(names(areaDat), names(den)))

    areaPerSample <- areaDat %>% 
                     group_by(Sample) %>%
                     summarize(SampleArea=sum(Area)) %>%
                     ungroup()

    datPerSample <- allVals %>% 
                    group_by(Sample,CellType) %>% 
                    summarize(SampleCounts=sum(Counts)) %>%
                    left_join(areaPerSample, by=c("Sample")) %>%
                    mutate(CellTypeSampleDensity=SampleCounts/SampleArea)

    allVals <- allVals %>%
               left_join(datPerSample, by=c("Sample","CellType"))

    return(allVals)
}


#' Calculate absolute and percentage density values by Sample, FOV, Band
#' 
#' Calculate absolute and percentage density values by Sample, FOV, Band
#' 
#' @param den
#' @param areaDat
#' @return table of all density values
#' @export
getAllInfiltrationDensityValues <- function(den, areaDat){
    if(!"BandArea" %in% names(areaDat)){
        interfaceArea <- areaDat %>% 
                         gather(3:ncol(areaDat), key="Band", value="BandArea")
    } else {
        interfaceArea <- areaDat
    }

    sbArea <- interfaceArea %>% 
              group_by(Sample,Band) %>% 
              summarise(SampleBandArea=sum(BandArea)) %>%
              ungroup()

    den <- den %>% 
           left_join(sbArea, by=c("Sample","Band")) %>%
           group_by(Sample,Band,CellType) %>%
           mutate(SampleBandCellTypeCounts=sum(Counts), 
                  SampleBandCellTypeDensity=SampleBandCellTypeCounts/SampleBandArea) %>% 
           ungroup()

    return(den)
}

#' Summarize positive/negative marker matrix by cell type definition
#'
#' Get cell counts and densities by cell definition rather than marker combination
#' 
#' @param assignments     table of 0s and 1s, where columns are markers (may be pos or neg)
#'                        and each row is a unique cell; zero indicates cell is NOT positive for that marker
#'                        (note, if marker name is negative a zero means the opposite) and one 
#'                        indicates positivity
#' @param bandArea        table indicating calculated area for each sample, fov and band
#' @param config          parsed marker config of info for each cell type definition including marker 
#'                        combinations (result of getMarkerConfig())
#' @param plotConfigFile  full path to file containing plot options including colors for cell types
#' @return list of two items: den=density table, config=updated config to represent summary info
#' @export
summarizeByCellTypeDefinition <- function(assignments, bandArea, config, plotConfigFile){
    
    pltCfg <- read_yaml(plotConfigFile)

    if(!"BandArea" %in% names(bandArea)){
        bandArea <- bandArea %>% gather(3:ncol(bandArea), key="Band", value="BandArea")
    }

    assignments$CellType <- "Unknown"
    for(ct in unique(config$CellType)){
        def <- unique(config$MarkerSet[which(config$CellType == ct)])
        if(length(def) != 1){
            flog.warn(paste0("Cell type ",ct," belongs to multiple marker sets. Can not choose for defining cells."))
            next
        }
        tmp <- assignments %>% select(unlist(strsplit(ct,",")))
        idx <- which(rowSums(tmp) == ncol(tmp))
        assignments$CellType[idx] <- def
    }

    den <- assignments %>% 
           filter(CellType != "Unknown") %>%
           select(Sample,SPOT,Band,CellType) %>% 
           group_by(Sample,SPOT,Band,CellType) %>% 
           summarise(Counts=n()) %>%
           left_join(bandArea, by=c("Sample","SPOT","Band")) %>%
           mutate(Density=Counts/BandArea)

    ## revamp config
    config$CellType <- config$MarkerSet
    config$Label <- config$MarkerSet
    config$MarkerSet <- "Cell_Identity"
    config$MarkerSetAlias <- "Cell_Identity"
    config$Population <- "Cell_Identity"
    config$Color <- "gray"
    config <- unique(config)
    for(lbl in config$Label){
        if(!lbl %in% names(pltCfg[["plot_options"]][["marker_colors"]])){ next }
        config$Color[which(config$Label == lbl)] <- pltCfg[["plot_options"]][["marker_colors"]][[lbl]]
    }        

    return(list(den=den, config=config))        
}


#' Summarize density data by cell type definition
#' 
#' Summarize density data by cell type definition
#' 
#' @param den            density table
#' @param config         parsed marker config of info for each cell type definition including marker 
#'                         combinations (result of getMarkerConfig())
#' @param areas          fovAreas
#' @param plotConfigFile full path to file containing plot options including colors for cell types
#' @return list of two items: den=density table, config=updated config to represent summary info
#' @export
summarizeFOVDataByCellTypeDefinition <- function(den, config, areas, plotConfigFile){

    pltCfg <- read_yaml(plotConfigFile)

    den$Definition <- "Unknown"
    for(ct in unique(config$CellType)){
        def <- unique(config$MarkerSet[which(config$CellType == ct)])
        if(length(def) != 1){
            flog.warn(paste0("Cell type ",ct," belongs to multiple marker sets. Can not choose for defining cells."))
            next
        }
        den$Definition[which(den$CellType == ct)] <- def
    }
  
    den$CellType <- den$Definition
    den <- den %>% 
           filter(CellType != "Unknown") %>%
           group_by(CellType,Sample,SPOT) %>%
           summarize(Counts=sum(Counts))

    den <- den %>%
           left_join(areas, by=c("Sample","SPOT")) %>%
           mutate(Total=Counts/Area) 

    config$CellType <- config$MarkerSet
    config$Label <- config$MarkerSet
    config$MarkerSet <- "Cell_Identity"
    config$MarkerSetAlias <- "Cell_Identity"
    config$Population <- "Cell_Identity"
    config$Color <- "gray"
    config <- unique(config)
    for(lbl in config$Label){
        if(!lbl %in% names(pltCfg[["plot_options"]][["marker_colors"]])){ next }
        config$Color[which(config$Label == lbl)] <- pltCfg[["plot_options"]][["marker_colors"]][[lbl]]
    }

    return(list(den=den, config=config))
}


#' Calculate total area in each FOV
#'
#' Calculate total area in each FOV
#' 
#' @param dat
#' @param metaFiles
#' @param areaDir
#' @param haloAnnotations
#' @param annotationsDirs
#' @param writeCSVfiles
#' @param maxG
#' @return table of area values for all FOVs
#' @export
calculateAreaTotalFOV <- function(dat, metaFiles, areaDir, haloAnnotations=NULL, annotationsDirs=NULL, writeCSVfiles=TRUE, maxG=5){

    if(is.null(haloAnnotations) && is.null(annotationsDirs)){
        stop("Can not calculate total FOV areas. Must provide either pre-parsed halo annotations (haloAnnotations) or vector of annotations directories (annotationsDirs.")
    }

    sampAnnFile <- metaFiles[grep("Sample",metaFiles)]
    sampAnn <- read.xlsx(sampAnnFile,1)

    areas <- list()

    for(s in unique(dat$Sample)){
        flog.debug(paste0("Sample: ",s))
        dd <- dat %>% filter(Sample == s)

        annFiles <- NULL
        if(is.null(haloAnnotations) && !is.null(annotationsDirs)){
            annFiles <- NULL
            for(annDir in annotationsDirs){
                flog.debug(paste0(" Annotation dir: ",annDir))
                if(length(grep(paste0(s,"_"),dir(annDir))) > 0){
                    annFiles <- c(annFiles, file.path(annDir,dir(annDir)[grep(paste0(s,"_"),dir(annDir))]))
                    flog.debug("  Annotation files: ")
                    cat(c(paste0(annFiles,collapse="\n    "),"\n"))
                } else {
                    flog.debug("  no annotation files found for this sample.")
                    next
                }
            }
        } 
        for(fov in unique(dd$SPOT)){
            flog.debug(paste0("  FOV: ",fov))
            excB <- NULL
            if(is.null(haloAnnotations) && !is.null(annotationsDirs)){
                aFiles <- NULL
                afi <- grep(paste0("Spot",fov,".annotations"),annFiles)
                if(length(afi) > 0){
                    aFiles <- annFiles[afi]
                    for(af in aFiles){
                        boundaries <- cleanBoundaries(readHaloAnnotations(af))
                        excB <- c(excB, boundaries$excB, boundaries$epiB, boundaries$glsB)
                    }
                }
            } else {
                boundaries <- haloAnnotations[[s]][[as.character(fov)]]
                excB <- c(boundaries$excB, boundaries$epiB, boundaries$glsB)
            }
            ds <- dd %>% filter(SPOT==fov)
            ds <- convertCellMinMaxToMidpoints(ds) ## this can be done in removeExcludedPoints too 
                                                   ## but we need it done even if there arent any
            ds <- trimDFToBB(ds,bb) ## needs X and Y coords
            gg <- generateSobolGrid(10^maxG,bb)

            if(len(excB) > 0){
                ## remove excluded points
                ds <- removeExcludedPoints(ds,excB=excB)
                excludedPts <- pointsInsidePolygonList(gg,excB)
                gg <- gg[!excludedPts,]
            }
            totalArea <- p2tomm * areaBB(bb) * nrow(gg)/10^maxG
         
            areas[[len(areas)+1]] <- data.frame(Area=totalArea) %>%
                                as_tibble %>%
                                mutate(Sample=s,SPOT=fov)
        }
        areass <- NULL

        flog.debug("Num elements in list areas: %s",length(areas))
        if(length(areas) > 0){
            areass <- areas %>% bind_rows
        }

        if(writeCSVfiles){
            write_csv(as.data.frame(areass), file.path(areaDir,"All_samples_area.csv"))
        }
    }
    return(areass)
}


#' Calculate marker densities for full FOVs
#'
#' Calculate marker densities for full FOVs
#' 
#' @param dat
#' @param areas
#' @param markers
#' @param writeCSVfiles
#' @param densityDir
#' @return table of marker densities per FOV
#' @export
calculateMarkerDensityTotalFOV <- function(dat, areas, markers, writeCSVfiles=TRUE, densityDir=NULL){
    rho <- NULL
    for(s in unique(dat$Sample)){
        dd <- dat %>% filter(Sample == s)
        for(fov in unique(dd$SPOT)){
            
            ds <- dd %>% filter(SPOT==fov)
            ds <- ds %>% spread(Marker,Value)
            indivMarkers <- unique(gsub("-","",unique(unlist(strsplit(markers,",")))))
            
            for(mi in indivMarkers) {
                if(mi %in% names(ds)){
                    ds[[paste0(mi,"-")]] <- ifelse(ds[[mi]]==0,1,0)
                } else {
                    ##### FOR EXCLUDED MARKERS, SET BOTH POS AND NEG TO 0 ???
                    ds[[mi]] <- NA
                    ds[[paste0(mi,"-")]] <- NA
                }
            }

            ## get total number of cells
            #mt <- computeMultiMarkerTable(ds, markers) %>% filter(CD3==1) ### NOT SURE ABOUT THIS, OMITTING FOR NOW
            mt <- computeMultiMarkerTable(ds, markers)
            mtCounts <- mt %>% dplyr::select(markers) %>% colSums
            totalCounts <- as.matrix(mtCounts)

            totalArea <- areas %>% filter(Sample == s, SPOT == fov) %>% pull(Area)
  
            rho[[len(rho)+1]] <- data.frame(CellType=rownames(totalCounts),
                                        Counts=totalCounts,
                                        Total=totalCounts/totalArea,
                                        Sample=s,
                                        SPOT=fov) %>%
                                        as.tibble

        }
        rhos <- NULL
        flog.debug("Num elements in list rho: %s",length(rho))
        if(length(rho) > 0){
            rhos <- rho %>% bind_rows
        }

        if(writeCSVfiles){
            write_csv(as.data.frame(rhos), file.path(densityDir,"All_samples_density.csv"))
        }
    }

    return(rhos)

}

