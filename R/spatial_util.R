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
#' @param boundaryColors a list including hex colors for Tum, Exc, and Epi boundaries; 
#'                       default=list("65280"="Tum", "65535"="Exc", "255"="Epi")
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

    boundaryColors <- list("65280"="Tum", "65535"="Exc", "255"="Epi")

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
                          boundaryColors[[LineColor]], sep="\t")
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
#' @param annoteFile     XML file containing boundary annotations for a single FOV
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

#' Get min and max X and Y coordinates of a FOV
#' 
#' @param  spObj  spot object(?); data frame containing Xmin,Xmax and Ymin,Ymax 
#'                coordinates of all cells in a FOV
#' @return list where X0 = minimum X, X1 = maximum X, and same for Y 
getBoundingBox <- function(haloObj){
    list(X0=min(haloObj$XMin),
         X1=max(haloObj$XMax),
         Y1=-min(haloObj$YMin),
         Y0=-max(haloObj$YMax))
}

#' Replace Min/Max X and Y coordinates from Halo with cell midpoints
#'
#' Given a tibble containing columns XMax, Xmin, YMax, YMin, replace
#' those columns with 'X' and 'Y' columns which contain the midpoint
#' between those min and max values
#' 
#' @param ds  tibble containing columns XMax, Xmin, YMax, YMin, Sample,
#'            SPOT, UUID, Marker, Value
#' @return  tibble containing only those columns mentioned above
convertCellMinMaxToMidpoints <- function(ds){
    ds %>%
        mutate(X=(XMax+XMin)/2,Y=-(YMax+YMin)/2) %>%
        dplyr::select(Sample,SPOT,UUID,X,Y,Marker,Value)

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

findDistanceFromPointsToInterfacePoint<-function(pts,interface) {

    interfaceTree=createTree(interface)

    knn=knnLookup(interfaceTree,newdat=pts,k=1)

    as.vector(apply(cbind(pts,interface[knn,]),1,function(pair){distance(pair[c(1,2)],pair[c(3,4)])}))

}

getContours<-function(df,levels){

    contours=getContourLines(df,levels=levels)
    split(contours,contours$Group) %>% map(rename,X=x,Y=y,Z=z)

}

