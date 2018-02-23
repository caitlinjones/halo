#' Calculate area and density of interface 
#' 
#' calculate interface stats
#'
#' @param aFiles             a vector of Halo boundary annotation files, with 
#'                           '.annotations' extension in XML format
#' @param dataFile           a *.rda file with Halo data for one sample
#' @param cellTypeFile       file containing cell types to be calculated; each line
#'                           is one cell type and each cell type can be either a single
#'                           marker name or a comma-separated list of markers
#' @param cellTypeName       a name to represent the cell types in marker file
#' @param funcMarker         this marker will be added to all markers in cellTypeFile in order to
#'                           compare cells that are +/- for this functional marker [ TO DO: REWORD THIS ]
#' @param fovBoundingBox     a list with four elements: X0, X1, Y0, Y1, representing the minimum
#'                           and maximum coordinates of the FOV to be plotted
#' @param plotBoundingBox    a list with four elements: X0, X1, Y0, Y1, representing the minimum
#'                           and maximum coordinates of the outermost plot box [ TO DO: UPDATE THIS ] 
#' @param pad                amount to trim from FOV before calculating
#' @param writeCSVfiles      logical indicating whether to write density and area tables to file; 
#'                           default=TRUE
#' @return a list containing two tables: 'density' and 'area' 
#' @export
calculateInterfaceStats <- function(aFiles, dataFile, cellTypeFile, cellTypeName,
                                    fovBoundingBox, pad, plotBoundingBox, funcMarker=NULL,
                                    writeCSVfiles=TRUE){
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

  interfaceBins <- (-20:20)*10

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
      ds <- getBandCounts(ds, bbData,bbPlot,bbFOV, aFile, doPlot=F)
      if(is.null(ds)){ 
        next
      }

      flog.debug("filtering band count data for CD3+ % counting bands")
      bdat <- ds %>% filter(CD3==1) %>% filter(!is.na(Band)) %>% count(Band) %$% data.frame(Band,n)

      #### TO DO: LOOK AT THE CODE THAT GENERATES THE FILES
      #### HARDCODED IN JOINBANDAREA (interfaceTools.R)
      flog.debug("joining band area")
      bdat <- joinBandArea(bdat,sampleName,spot)

      markerSet <- scan(cellTypeFile,"")
      if(!is.null(funcMarker)){
        markerSet <- c(markerSet,paste0(markerSet,",",funcMarker))
      }
      fillNA <- list()
      dummy <- lapply(markerSet, function(x){ fillNA[x] <<- 0 })
      flog.debug("computing multimarker table")
      mt <- computeMultiMarkerTable(ds,markerSet)
      mtCounts <- mt %>%
                  group_by(Band) %>%
                  summarize_at(markerSet,sum) %>%
                  complete(Band,fill=fillNA)
      mtB <- as.matrix(mtCounts[!is.na(mtCounts$Band),2:ncol(mtCounts)],check.names=F)
      rownames(mtB) <- mtCounts$Band[!is.na(mtCounts$Band)]

      mtBCollapse <- apply(mtB,2,function(x){tapply(x,interfaceSide(rownames(mtB)),sum)})
      bandAreaCollapse <- tapply(bdat$Area,interfaceSide(bdat$Band),sum)
      MM <- diag(1/bandAreaCollapse)
      rownames(MM) <- c("Inside","Outside")
      colnames(MM) <- c("Inside","Outside")

      totalCounts <- colSums(mtBCollapse)
      totalArea <- sum(bandAreaCollapse)

      flog.debug("getting density")
      ## density
      rho[[len(rho)+1]] <- (t(mtBCollapse) %*% MM) %>%
                            as.data.frame %>%
                            cbind(data.frame(Counts=t(mtBCollapse))) %>%
                            rownames_to_column("CellType") %>%
                            mutate(Total=totalCounts/totalArea,Sample=sampleName,SPOT=spot)
      flog.debug("getting area")
      ## area
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
    ctf <- gsub("\\.txt","",cellTypeFile)
    write_csv(as.data.frame(areass), cc(sampleName, ctf, "interfaceStatsV3_Area.csv"))
    write_csv(as.data.frame(rhos), cc(sampleName, ctf, "interfaceStatsV3_Density.csv"))
  }

  return(list(density=rhos, area=areass))
}





#' Recursively process an XML node for one FOV
#'
#' Visit each node in an object of class XMLNode and gather information
#' to print one line for each vertex which represents the location of a cell
#' on the slide [NOTE: DOUBLE CHECK THIS FOR ACCURACY] 
#' 
#' Output is a tab-delimited file containing: 
#'     1)  FOV name
#'     2)  Annotation number
#'     3)  Line color
#'     4)  Name
#'     5)  Visible
#'     6)  Region number
#'     7)  Region type
#'     8)  HasEndcaps
#'     9)  NegativeROA
#'     10) X coordinate
#'     11) Y coordinate ***(negative of what is output by Halo)***
#'     12) Region according to color (Tumor, Exclusion, Epidermis)
#' 
#' @param node           object of class XMLNode containing Halo annotations for one FOV
#' @param fovName        unique identifier of FOV
#' @param outFile        output file
#' @param LineColor      code for line color of the annotation ["655280"|"65535"|"255"]
#' @param Name           name of the annotation
#' @param Visible        logical indicating whether annotation is visible [????]
#' @param RegionType     region type *Note: any region with type other than Polygon is skipped
#' @param HasEndcaps     [???]
#' @param NegativeROA    [???]
#' @return nothing
parseAnnotationXML <- function(node, fovName=NULL, outFile=NULL, LineColor="", 
                               Name="", Visible="", RegionType="", HasEndcaps=0,
                               NegativeROA=0){
  ## use global variables in order to keep track of annotation and region nums
  ## while traversing xml doc; remove them later
  annotationNum <<- 0
  regionNum <<- 0

  processNode <- function(node, fovName=NULL, outFile=NULL, 
                        LineColor="", Name="", Visible="", RegionType="", HasEndcaps=0,
                        NegativeROA=0){

    if(is.null(node) | length(node) == 0){
        return()
    }
    colorToRegionMap <- list("65280"="Tum", "65535"="Exc", "255"="Epi")

    nodeName <- xmlName(node)
    nodeAttrs <- xmlAttrs(node)
    
    if(nodeName == "Annotation"){
        annotationNum <<- annotationNum + 1
        flog.debug("AnnotNum: %s",annotationNum)
        Name <- xmlGetAttr(node, "Name")
        Visible <- xmlGetAttr(node, "Visible")
        LineColor <- xmlGetAttr(node, "LineColor")
    } else if(nodeName == "Region"){
        regionNum <<- regionNum + 1
        flog.debug("RegionNum: %s",regionNum)
        HasEndcaps <- xmlGetAttr(node, "HasEndcaps")
        NegativeROA <- xmlGetAttr(node, "NegativeROA")
        RegionType <- xmlGetAttr(node, "Type")
    } else if(nodeName == "V"){
        X = xmlGetAttr(node, "X")
        Y = -as.numeric(xmlGetAttr(node, "Y"))
        if(RegionType == "Polygon"){
            line <- paste(fovName, annotationNum, LineColor, Name, Visible, 
                          regionNum, RegionType, HasEndcaps, NegativeROA, X, Y, 
                          colorToRegionMap[[LineColor]], sep="\t")
            write(line, file=outFile, append=TRUE)
        }
        return()
    }

    if(xmlSize(node) > 0){
       for(i in 1:xmlSize(node)){
           processNode(node[[i]], fovName=fovName, LineColor=LineColor, 
                       Name=Name, Visible=Visible, RegionType=RegionType, 
                       HasEndcaps=HasEndcaps, NegativeROA=NegativeROA,
                       outFile=outFile)
       }
    }
  }

  res <- processNode(node, fovName=fovName, outFile=outFile, 
           LineColor=LineColor, Name=Name, Visible=Visible, RegionType=RegionType, 
           HasEndcaps=HasEndcaps, NegativeROA=NegativeROA)

  suppressWarnings(rm(annotationNum))
  suppressWarnings(rm(regionNum))
}

#' Parse boundary annotations of a single FOV in a single XML file
#' 
#' Parse XML file containing annotations for one FOV. Write tab-delimited file
#' of select information (see ?parseAnnotationXML for details)
#' 
#' @param annoteFile    XML file containing boundary annotations for a single FOV
#' @export
readHaloAnnotations <- function(annoteFile){
    tfile <- tempfile(tmpdir=".")
    header = paste("File","AnnotationNum","LineColor","Name","Visible","RegionNum","Type",
                   "HasEndcaps","NegativeROA","X","Y","RegionCode",sep="\t")
    write(header, file=tfile, append=FALSE)

    rootNode <- xmlRoot(xmlParse(annoteFile))
    fovName <- gsub("\\.annotations","",basename(annoteFile))
    parseAnnotationXML(rootNode, fovName=fovName, outFile=tfile) 

    aa <- read_tsv(tfile, col_types="ciicciciiiic")
    fixFile <- file.path(dirname(annoteFile), "haloBoundaryRegionReassignments.csv")
    if(file.exists(fixFile)){
        flog.debug("Reassign file=%s",fixFile)
        fix <- read_csv(fixFile, col_types="ccic")
        ff <- fix %>% filter(Spot==aa$File[1])
        if(nrow(ff)>0){
            for(kk in seq(nrow(ff))){
                aa$RegionCode[ aa$RegionNum==ff$RegionNum[kk] ] <- ff$RegionCode[kk]
            }
        }
    }
    file.remove(tfile)
    split(aa,aa$RegionNum)
}


#' Get min and max X and Y coordinates of a FOV
#' 
#' @param  spObj  spot object(?); data frame containing X and Y coordinates
#'                of all cells in a FOV
#' @return list where X0 = minimum X, X1 = maximum X, and same for Y 
getBoundingBoxL<-function(spObj) {
    list(X0=min(spObj$X),Y0=min(spObj$Y),X1=max(spObj$X),Y1=max(spObj$Y))
}

#' Get X and Y coordinates of box enclosing two existing boxes
#'
joinBoundingBoxes<-function(b1,b2){
    list(
        X0=min(b1$X0,b2$X0),
        Y0=min(b1$Y0,b2$Y0),
        X1=max(b1$X1,b2$X1),
        Y1=max(b1$Y1,b2$Y1)
        )
}

interfaceSide<-function(bnd) {
    ifelse(grepl("-",bnd),"Inside","Outside")
}

drawBoundaries<-function(bnd,textLabels=T) {
    for(jj in seq(bnd)) {
        drawBoundary(bnd[[jj]],textLabels=textLabels)
    }
}

drawBoundary<-function(b,bCols=bndColors,col=NULL,lwd=2,textLabels=T,...) {
    if(is.null(col)) {
        bCol=bCols[[b$RegionCode[1]]]
    } else {
        bCol=col
    }
    lines(b$X,b$Y,lwd=lwd,col=bCol,...)
    if(textLabels) {
        text(mean(b$X),mean(b$Y),b$RegionNum[1],cex=1.4,col=bCol)
    }
}

drawBoundingBox<-function(b,...){
    rect(b$X0,b$Y0,b$X1,b$Y1,...)
}

padBoundingBox<-function(bb,pad){
    bn=list()
    bn$X0=bb$X0-pad
    bn$Y0=bb$Y0-pad
    bn$X1=bb$X1+pad
    bn$Y1=bb$Y1+pad
    bn
}

#' Distance between points (L2-norm)
distance<-function(p1,p2){sqrt(sum((p1-p2)^2))}

#' Compute theta
theta<-function(p1,p2){vv=p1-p2;atan(vv[2]/vv[1])*(180/pi)}

# Value:

#      integer array; values are: 0: point is strictly exterior to pol;
#      1: point is strictly interior to pol; 2: point lies on the
#      relative interior of an edge of pol; 3: point is a vertex of pol.


pointsInPolygon<-function(pts,poly){
    point.in.polygon(pts$X,pts$Y,poly$X,poly$Y)==1
}

pointsInsidePolygonList<-function(pts,bList) {

    inside=bList %>%
        map(function(bn){point.in.polygon(pts$X,pts$Y,bn$X,bn$Y)>0}) %>%
        as.data.frame %>%
        apply(.,1,any)

    inside

}

intersectWithBoundingBox<-function(x,bb) {
     as.tibble(x) %>%
        mutate(X=ifelse(X>bb$X1,bb$X1,ifelse(X<bb$X0,bb$X0,X))) %>%
        mutate(Y=ifelse(Y>bb$Y1,bb$Y1,ifelse(Y<bb$Y0,bb$Y0,Y))) %>%
        distinct
}

expandClippedPathToOuterBox<-function(x,bin,bout){
    x %>%
        mutate(X=ifelse(X>bin$X1,bout$X1,ifelse(X<bin$X0,bout$X0,X))) %>%
        mutate(Y=ifelse(Y>bin$Y1,bout$Y1,ifelse(Y<bin$Y0,bout$Y0,Y))) %>%
        distinct
}

boundingBoxToRect<-function(bb){
    rr=rbind(c(bb$X0,bb$Y0),c(bb$X1,bb$Y0),c(bb$X1,bb$Y1),c(bb$X0,bb$Y1))
    colnames(rr)=c("X","Y")
    as.data.frame(rr)
}

boundingBoxToPath<-function(bb,Z=10){
    poly=rbind(
        cbind(bb$X0:bb$X1,bb$Y0),
        cbind(bb$X1,bb$Y0:bb$Y1),
        cbind(bb$X1:bb$X0,bb$Y1),
        cbind(bb$X0,bb$Y1:bb$Y0))
    colnames(poly)=c("X","Y")
    as.data.frame(poly)
}

areaBB<-function(bb){
    (bb$X1-bb$X0)*(bb$Y1-bb$Y0)
}

boundaryToSpatialPolygon<-function(bnd,tag) {
    SpatialPolygons(list(
         Polygons(list(Polygon(bnd[,c("X","Y")])),tag)))
}

clipSPBoundaryToBB<-function(boundarySP,bb) {

    box=boundingBoxToRect(bb)
    spClip=SpatialPolygons(list(
        Polygons(list(Polygon(box)),"Clip")))

    res=try({clipBnd=gIntersection(boundarySP,spClip)})
        if(class(res)=="try-error") {
            cat("gIntersect::ERROR\n")
            bnd.s=gSimplify(boundarySP,tol=0.0000001)
            bnd.s=gBuffer(bnd.s,byid=T,width=0)
            clipBnd=gIntersection(bnd.s,spClip)
        }

    if(!is.null(clipBnd)) {
        if(class(clipBnd)=="SpatialCollections"){
            clipBnd=clipBnd@polyobj
        }
        return(clipBnd)
    } else {
        return(boundarySP)
    }

}

clipSPBoundaryToBBWithDraw<-function(boundarySP,bb,bColor=4) {

    clipBnd=clipBoundaryToBB(boundarySP,bb)

    if(!is.null(clipBnd)) {
        lines(clipBnd,col=bColor,lwd=3)
    }

    return(clipBnd)

}

#
# Note trim works on data frames as the result may be open
#
trimBoundaryToBB<-function(boundaryDF,bb) {

     xx=boundaryDF %>%
     filter(X>=bb$X0 & X<=bb$X1) %>%
     filter(Y>=bb$Y0 & Y<=bb$Y1)
     xx

}

trimDFToBB<-function(df,bb) {

     xx=df %>%
     filter(X>=bb$X0 & X<=bb$X1) %>%
     filter(Y>=bb$Y0 & Y<=bb$Y1)
     xx

}

subDivideLongSegments<-function(boundary,maxSegLength) {

    bline=boundary %>% dplyr::select(X,Y) %>% mutate(SegNo=seq(nrow(.)))
    newPts=list()
    for(is in seq(nrow(bline)-1)) {

        pt1=bline[is,]
        pt2=bline[is+1,]
        subDivide=ceiling(distance(pt1,pt2)/maxSegLength)
        if(subDivide>1) {
            ##cat(is,"Dist=",distance(pt1,pt2),"subDivide =",subDivide,"\n")
            dx=pt2$X-pt1$X
            dy=pt2$Y-pt1$Y

            newPts[[is]]=data.frame(
                    X=pt1$X+dx*seq(subDivide-1)/subDivide,
                    Y=pt1$Y+dy*seq(subDivide-1)/subDivide,
                    SegNo=is+seq(subDivide-1)/subDivide
                    )
        }
    }

    bind_rows(bline,bind_rows(newPts)) %>% arrange(SegNo)

}

findDistaceFromPointsToInterfacePoint<-function(pts,interface) {

    interfaceTree=createTree(interface)

    knn=knnLookup(interfaceTree,newdat=pts,k=1)

    as.vector(apply(cbind(pts,interface[knn,]),1,function(pair){distance(pair[c(1,2)],pair[c(3,4)])}))

}

getContours<-function(df,levels){

    contours=getContourLines(df,levels=levels)
    split(contours,contours$Group) %>% map(rename,X=x,Y=y,Z=z)

}

removeExcludedPoints<-function(ds,aFile,doPlot=F) {

    ds=ds %>%
        mutate(X=(XMax+XMin)/2,Y=-(YMax+YMin)/2) %>%
        dplyr::select(Sample,SPOT,UUID,X,Y,Marker,Value)

    boundaries=readHaloAnnotations(aFile)
    regionTable=boundaries %>% bind_rows %>% count(RegionNum,RegionCode)
    excB=boundaries[regionTable$RegionCode!="Tum"]

    spCells=ds %>% dplyr::select(UUID,X,Y)
    coordinates(spCells)=~X+Y

    if(len(excB)>0){
        spExcB=seq(excB) %>% map(function(b){boundaryToSpatialPolygon(excB[[b]],b)})
        excludedCellIdx=NULL
        for(jj in len(spExcB)) {
            excCellJJ=unlist(over(spExcB[[jj]],geometry(spCells),returnList=T))
            excludedCellIdx=union(excludedCellIdx,excCellJJ)
        }
        #points(spCells[excludedCellIdx,],col=3,cex=.7)
        if(len(excludedCellIdx)>0)
            ds=ds[-excludedCellIdx,]
    }

    return(ds)

}

#' Plot total density of each FOV for one sample
#' 
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
#' @return nothing
#' @export
plotTotalDensity <- function(dataFile, annotationsDir, cellTypeFile, cellTypeName,
                             fovBB, pad, plotBB, pdfFile, ymax=NULL,
                             logPlot=F, funcMarker="KI67", sortByMarker="CD3", sampleColor="#f16913", 
                             sampleColorDark="#80380A", exclude_sample_fov=NULL, exclude_sample_marker=NULL,
                             writeCSVfiles=TRUE){

    aFiles <- file.path(annotationsDir,dir(annotationsDir)[grep("\\.annotations$",dir(annotationsDir))])

    ## get density and area
    flog.debug("calculating interface stats")
    is <- calculateInterfaceStats(aFiles, dataFile, cellTypeFile, cellTypeName,      
                                 fovBB, pad, plotBB, writeCSVfiles=writeCSVfiles, funcMarker=funcMarker)
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
