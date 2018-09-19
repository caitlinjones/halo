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
                               NegativeROA=0, boundaryColors=NULL){
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

    if(is.null(boundaryColors)){
        ##BOUNDARY COLORS FOR COHORT 2
        boundaryColors <- list("65535"="Exc", "65280"="Epi", "255"="Gls", "16776960"="Tum") ### color for Tum is a guess for now!
    }

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
readHaloAnnotations <- function(annoteFile,boundaryColors=NULL,boundaryReassignmentFile=NULL){

    tfile <- tempfile(tmpdir=".")
    header = paste("File","AnnotationNum","LineColor","Name","Visible","RegionNum","Type",
                   "HasEndcaps","NegativeROA","X","Y","RegionCode",sep="\t")
    write(header, file=tfile, append=FALSE)

    rootNode <- xmlRoot(xmlParse(annoteFile))
    fovName <- gsub("\\.annotations","",basename(annoteFile))
    parseAnnotationXML(rootNode, fovName=fovName, outFile=tfile, boundaryColors=boundaryColors)

    aa <- read_tsv(tfile, col_types="ciccciciiiic")
    fixFile <- boundaryReassignmentFile 
    if(!is.null(fixFile) && file.exists(fixFile)){
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

#' Remove exclusion boundaries that are contained in another one
#' 
#' Remove data from boundaries table that represent exclusion boundaries
#' that are completely surrounded by another exclusion boundary
#' 
#' @param boundaries            table generated by readHaloAnnotations
#'                             
cleanBoundaries <- function(boundaries){

    regionTable <- boundaries %>% bind_rows %>% count(RegionNum,RegionCode)
    excB <- boundaries[regionTable$RegionCode=="Exc"]
    epiB <- boundaries[regionTable$RegionCode=="Epi"]
    glsB <- boundaries[regionTable$RegionCode=="Gls"]
    tumB <- boundaries[regionTable$RegionCode=="Tum"]

    excList <- list(excB=excB, epiB=epiB, glsB=glsB)
    for(x in names(excList)){
        ex <- excList[[x]]
        if(!is.null(ex) && len(ex) > 1){
            ## make list of spatial polygons, one element for each exclusion boundary
            spExcB <- seq(ex) %>% map(function(b){boundaryToSpatialPolygon(ex[[b]],b)})

            ## are any contained within another?
            containedBoundary <- rep(FALSE,len(ex))
            for(i in seq(len(ex))){
                containedBoundary[i] <- spExcB[-i] %>% map(gContains,spExcB[[i]]) %>% unlist %>% any
            }
            ## remove those
            excList[[x]] <- ex[!containedBoundary]
        }
    }
    excList[["tumB"]] <- tumB
    return(excList)
}


getAllHaloAnnotations <- function(samp, annotationsDirs, boundaryColors=NULL, boundaryReassignmentFile=NULL){

    allFiles <- NULL
    for(ad in annotationsDirs){
        aType <- basename(ad)
        fs <- file.path(ad, dir(ad)[grep(paste0("^",samp,"_"),dir(ad))])
        if(length(fs) > 0){
            for(f in fs){
                fov <- gsub(".*Spot","",f)
                fov <- gsub(".annotations","",fov)
                allFiles <- allFiles %>% bind_rows(list(File=f,Sample=samp,FOV=fov,BoundaryType=aType))
            }
        }
    }

    btCodes <- list('InterfaceCoordinates'='tumB', 'ExclusionCoordinates'='excB', 
                    'EpidermisCoordinates'='epiB', 'GlassCoordinates'='glsB')
    allAnnotations <- list()
    if(is.null(allFiles) || nrow(allFiles) == 0){
        return(allAnnotations)
    }
    for(fov in unique(allFiles$FOV)){
        fovBoundaries <- list()
        aFiles <- allFiles %>% filter(FOV == fov) %>% select(File,BoundaryType)

        for(x in 1:nrow(aFiles)){
            af <- aFiles$File[x]
            flog.debug(paste0("getting boundaries from file ",af))

            ##### TEMPORARY, FOR COHORT 1
            boundaryType <- aFiles$BoundaryType[x]
            bt <- btCodes[[boundaryType]]
            if(!is.null(boundaryColors)){
                boundaryColors <- list('65280'='Epi','65535'='Exc','255'='Gls')
                if(boundaryType == "InterfaceCoordinates"){
                    if(samp == "Untreated"){
                        boundaryColors <- list('65280'='Tum','65535'='IGNORE', '255'='IGNORE')
                    } else if(samp == "PR"){
                        if(as.numeric(fov) %in% c(3,13,14,21)){ 
                            boundaryColors <- list('65535'='Tum', '65280'='Exc', '255'='IGNORE')
                        } else {
                            boundaryColors <- list('65535'='IGNORE', '65280'='Tum')
                        }
                    }
                }
            } 

            ##########################           

            boundaries <- cleanBoundaries(readHaloAnnotations(af, boundaryColors=boundaryColors, 
                                                              boundaryReassignmentFile=boundaryReassignmentFile))
            fovBoundaries[[bt]] <- boundaries[[bt]] 
        }

        allAnnotations[[fov]] <- fovBoundaries
    }

    allAnnotations

}


#' Remove points from data set that fall inside exclusion boundaries
#' 
#' Read Halo boundaries annotations file and remove from data set any cells
#' that fall inside the exclusion regions
#' 
#' @param ds     tibble containing at minimum, Sample, SPOT, UUID, X, Y, Marker, Value
#' @param excB   a list of excluded region tables as returned by cleanBoundaries;
#'               if NULL, must provide the original Halo annotations file in XML format
#' @param aFile  Halo boundaries annotation file in XML format; if NULL, must provide
#'               the list of excluded region tables as returned by cleanBoundaries()
#' @return tibble filtered to remove points that fall inside exclusion boundaries
getHaloExclusions <- function(fov, dd, excType, haloAnnotations=NULL, aFile=NULL, exc=NULL, boundaryColors=NULL,
                             boundaryReassignmentFile=NULL){

    ds <- dd %>% filter(SPOT == fov)
    if(!all(c('X','Y') %in% colnames(ds))){
        ds <- convertCellMinMaxToMidpoints(ds)
    }

    annCodes <- list(Exc='excB', Epi='epiB', Gls='glsB', Tum='tumB')
    if(is.null(exc)){
        if(is.null(haloAnnotations)){
            boundaries <- readHaloAnnotations(aFile,boundaryColors=boundaryColors, boundaryReassignmentFile=boundaryReassignmentFile)
            allBoundaries <- cleanBoundaries(boundaries) #remove exc regions completely surrounded by another
            excB <- allBoundaries[[annCodes[[excType]]]]
        } else {
            excB <- haloAnnotations[[as.character(fov)]][[annCodes[[excType]]]]
            #excB <- cleanBoundaries(split(excB, excB$RegionNum))
        }
        if(is.null(excB)){
            return(c())
        }
#        #excB <- split(excB, excB$RegionNum)
    }

    spCells <- ds %>% dplyr::select(UUID,X,Y)
    coordinates(spCells) <- ~X+Y

    excludedCellIdx <- NULL
    if(len(excB)>0){
        spExcB <- seq(excB) %>% map(function(b){boundaryToSpatialPolygon(excB[[b]],b)})
        for(jj in 1:len(spExcB)) {
            excCellJJ <- unlist(sp::over(spExcB[[jj]],geometry(spCells),returnList=T))
            excludedCellIdx <- union(excludedCellIdx,excCellJJ)
        }
    }
    fovHaloIdx <- which(dd$SPOT == fov)
    fovHaloIdx[excludedCellIdx]

}

#' Get indices of cells to be excluded because of reasons determined manually by a lab member
#' 
#' Get indices of cells to be excluded because of reasons determined manually by a lab member
#'
#' @param  fovAnnotations  tibble of FOV annotations for one sample, read from xlsx file
#'                         included in study meta data
#' @param  dat             full data set to be marked for exclusions
#' @return vector of dat indices to be marked for exclusion
getStudyExclusions <- function(fovAnnotations, dat){
    excl <- fovAnnotations %>% select(FOV_number, matches("exclusion")) 
    
    ## get indices of FOVs to be excluded
    fovExcl <- excl %>% filter(FOV_exclusion_postHalo == "X") %>% pull(FOV_number)
    datExcl <- which(dat$SPOT %in% fovExcl)

    ## get indices of markers to be excluded
    markerExcl <- excl %>% 
                  filter(Marker_exclusion != "" & !is.na(Marker_exclusion)) %>% 
                  select(FOV_number,Marker_exclusion)
    if(nrow(markerExcl) > 0){
        for(x in 1:nrow(markerExcl)){
            excl <- markerExcl[x,]
            fov <- pull(excl[1])
            markers <- gsub("\\s+","",as.character(pull(excl[2])))
            markers <- unlist(strsplit(markers,","))
            for(m in markers){
                datExcl <- c(datExcl, which(dat$SPOT == fov & dat$Marker == m))
            }
        }
    }
    datExcl
}

getDriftExclusions <- function(drift, dat, threshold, samp=NULL){
    if(length(unique(gsub("_Spot.*.afi","",drift$image_location))) > 1){
        ## filter drift for current sample
        drift <- drift[grep(paste0("[\\|/]",samp,"_"),tmp,ignore.case=TRUE),]
    }

    excl <- c()
    drft <- drift %>% select(Image_Location=image_location, XMin=x_min, XMax=x_max,
                             YMin=y_min, YMax=y_max, drift_loss_pixel_pct)
    tmp <- full_join(dat, drft, by=c("Image_Location","XMin","XMax","YMin","YMax"))
    excl <- which(!is.na(tmp$drift_loss_pixel_pct) & tmp$drift_loss_pixel_pct > threshold) 
    excl
}

#' Get indices of cell to be excluded because they lie outside the border padding 
#' 
#' Get indices of cell to be excluded because they lie outside the border padding 
#' 
#' @param fov        fov to be marked for exclusions
#' @param samp       sample to be marked
#' @param dat        data tibble containing data to be marked
#' @param borderPad  number in microns to trim from data
getBorderPaddingExclusions <- function(fov, samp, dat, borderPad){

    fovDat <- dat %>% filter(Sample==samp, SPOT==fov)
    
    #bbFOV <- list(X0=1,Y0=-3375,X1=5363,Y1=1)
    bb <- list(X0 = bbFOV$X0 + borderPad,
               X1 = bbFOV$X1 - borderPad,
               Y0 = bbFOV$Y0 + borderPad,
               Y1 = bbFOV$Y1 - borderPad)

    which( dat$Sample == samp & dat$SPOT == fov &
         ( dat$X < bb$X0 | dat$X > bb$X1 | dat$Y < bb$Y0 | dat$Y > bb$Y1)) 
    
}

#' Add exclusion reasons to a vector of existing exclusion reasons
#' 
#' Given an exclusion reason and the indices of cells to be excluded, 
#' properly add the reason to vector of existing exclusions
#' 
#' @param exclCol   vector of existing exclusions (EXCLUDE column from object analysis data)
#' @param reason    string to be added to exclude column indicating why a cell is to be excluded
#' @param idxs     vector of indices at which exclusion reason should be added
#' @return vector of exclusion reasons with new reason added
writeExclusionReasons <- function(exclCol,idxs,reason){
    for(idx in idxs){
        exclCol[idx] <- ifelse(exclCol[idx] == "", reason, paste0(exclCol[idx],",",reason)) 
    }
    exclCol
}

#' Convert shared drive path to server path, and add subdirectory for exclusion coordinates
#' 
#' Do a simple string substitution to convert one path to another, and add exclusion 
#' coordinates subfolder to all paths
#' 
#' @param paths  a vector of paths to be converted
#' @param from   string to be replaced
#' @param to     string to replace {from}
#' @return a vector equal in length to input vector, with paths converted and ExclusionCoordinates
#'         subfolder added 
convertBoundaryFilePaths <- function(paths, subdir, from=".*bic.mskcc.org", to="/ifs/tcga/socci/Multiomyx"){
    notNA <- !is.na(paths)
    pths <- gsub(from,to,paths)
    pths[notNA] <- gsub("\\\\","/",pths[notNA])
    pths[notNA] <- file.path(pths[notNA],subdir)
    pths
}

#' Plot cell locations of exclusions
#' 
#' Plot cell locations of exclusions
#' @param dd               data for one sample with EXCLUDE column
#' @param haloAnnotations  halo annotations in list format, for one sample and one or multiple FOV
#' @param iFiles           if not providing haloAnnotations, a vector of interface XML files
#' @param aFiles           if not providing haloAnnotations, a vector of exclusion XML files 
#' @param epFiles          if not providing haloAnnotations, a vector of epidermis XML files
#' @param gFiles           if not providing haloAnnotations, a vector of glass XML files
#' @param outDir           directory where plots should be saved
#' @param boundaryColors   list of boundary colors as they appear in halo annotation XML files
#'                         where the values are the associated boundary code (e.g., 'epiB', 'glsB', etc.)
#' @param boundaryReassignmentFile  comma-separated file indicating exceptions to the color assignments
#'                                  in boundaryColors
#' @export
plotExclusions <- function(dd,haloAnnotations=NULL,iFiles=NULL,aFiles=NULL,epFiles=NULL,gFiles=NULL,
                           outDir=getwd(),boundaryColors=NULL,boundaryReassignmentFile=NULL){

    if(!all(c('X','Y') %in% colnames(dd))){
        dd <- dd %>% mutate(X=(XMin+XMax)/2,Y=-(YMin+YMax)/2)
    }

    allAnns <- list(Tum=iFiles,Exc=aFiles,Epi=epFiles,Gls=gFiles)
    allAnnCodes <- list(Tum='tumB', Exc='excB', Epi='epiB', Gls='glsB')
    info <- data.frame(boundaryType=names(allAnns),
                       plotTitle=c("Tumor","HALOExclusion","HALOEpidermis","HALOGlass"),
                       boundaryColor=c("black","lightblue","orange","magenta"),
                       cellColor=c("black","yellow","green","blue"))

    samp <- unique(dd$Sample)
    for(fov in unique(dd$SPOT)){
        fovHaloAnn <- NULL
        if(as.character(fov) %in% names(haloAnnotations)){
            fovHaloAnn <- haloAnnotations[[as.character(fov)]]
        } else {
            next
        }
        jpeg(file.path(outDir,paste0(samp,"_",fov,"_exclusion_debug.jpeg")),width=11,height=11,units="in",res=180)
        par(mfrow=c(3,3))
        flog.debug(paste0("plotting FOV ",fov))
        tmp <- dd %>% filter(SPOT == fov)

        for(bt in names(allAnns)){
            inf <- info[which(info$boundaryType == bt),]
            plot(NA, xlim=c(bb$X0,bb$X1), ylim=c(bb$Y0,bb$Y1), main=inf$plotTitle)
            rect(bb$X0,bb$Y0,bb$X1,bb$Y1, border="black",lty="dashed")
            points(tmp$X,tmp$Y,col="gray",cex=0.15,pch=19)
            if(!inf$plotTitle == "Tumor"){
                haloExc <- tmp %>% filter(grepl(inf$plotTitle,EXCLUDE))
                points(haloExc$X,haloExc$Y,col=inf$cellColor,pch=19,cex=0.15)
            }

            if(is.null(haloAnnotations)){
                cha <- NULL
                allBounds <- NULL
                aFile <- allAnns[[bt]][grep(paste0("Spot",fov,".annotations"),allAnns[[bt]])]
                if(length(aFile == 1)){
                    ha <- readHaloAnnotations(aFile, boundaryColors=boundaryColors, 
                                              boundaryReassignmentFile=boundaryReassignmentFile)
                    cha <- cleanBoundaries(ha)
                }
                if(!is.null(cha)){
                    for(x in names(cha)){
                        for(y in names(cha[[x]])){
                            allBounds <- allBounds %>% bind_rows(cha[[x]][[y]])
                        }
                    }
                }
            } else {
                allBounds <- fovHaloAnn[[allAnnCodes[[bt]]]]
            }
            if(!is.null(allBounds) && length(allBounds) > 0){
                for(region in names(allBounds)){
                    bounds <- allBounds[[region]]
                    if(!is.null(bounds) && nrow(bounds) > 0){
                        color <- info$boundaryColor[which(info$boundaryType == bt)]
                        points(bounds$X,bounds$Y,col=color,cex=0.15,pch=19)
                    }    
                }
            }
        }

        plot(NA, xlim=c(bb$X0,bb$X1), ylim=c(bb$Y0,bb$Y1),main="PADDED_20")
        rect(bb$X0,bb$Y0,bb$X1,bb$Y1, border="black",lty="dashed")
        points(tmp$X,tmp$Y,col="gray",cex=0.15,pch=19)
        borderExc <- tmp %>% filter(grepl("PAD",EXCLUDE))
        points(borderExc$X,borderExc$Y,col="purple",cex=0.15,pch=19)

        plot(NA, xlim=c(bb$X0,bb$X1), ylim=c(bb$Y0,bb$Y1),main="DRIFT_10")
        rect(bb$X0,bb$Y0,bb$X1,bb$Y1, border="black",lty="dashed")
        points(tmp$X,tmp$Y,col="gray",cex=0.155555,pch=19)
        driftExc <- tmp %>% filter(grepl("DRIFT",EXCLUDE))
        points(driftExc$X,driftExc$Y,col="blue",cex=0.15,pch=19)

        plot(NA, xlim=c(bb$X0,bb$X1), ylim=c(bb$Y0,bb$Y1),main="LabExclusion")
        rect(bb$X0,bb$Y0,bb$X1,bb$Y1, border="black",lty="dashed")
        points(tmp$X,tmp$Y,col="gray",cex=0.15,pch=19)
        studyExc <- tmp %>% filter(grepl("Lab",EXCLUDE))
        points(studyExc$X,studyExc$Y,col="red",cex=0.15,pch=19)

        dev.off()
    }
}

#' Given a tibble of object analysis data, add a column EXCLUDE to indicate
#' which cells should be excluded from further analysis
#'
#' Given a tibble of object analysis data, add a column EXCLUDE to indicate
#' which cells should be excluded from further analysis, according to various upstream
#' factors including drift/loss, border padding, and other technical reasons
#' determined by investigators
#'
#' @param samp            sample name as it appears in both SampleAnnotations.xlsx and data
#' @param dat             tibble containing Halo object analysis data, to which EXCLUDE column
#'                        will be added
#' @param drift           tibble of drift/loss summary 
#' @param fovAnn          FOV annotations for one sample
#' @param cellDiveId      CELL_DIVE_ID mapping sample name to fov annotations
#' @param borderPad       number in pixels indicating the minimum distance between a cell and the FOV
#'                        border in order for that cell NOT to be excluded
#' @param driftThreshold  maximum pixel percent determined to have drifted in order to NOT be excluded
#' @return tibble matching {dat} exactly, with an extra column, EXCLUDE, added to indicate which
#'         cells should not be analyzed
#' @export
markExclusions <- function(samp, dat, drift, fovAnn, cellDiveId, haloAnn=NULL, borderPad=20, driftThreshold=0.1, 
                           printPlots=FALSE, debugDir=getwd(), boundaryColors=NULL, boundaryReassignmentFile=NULL){

    borderPad_um <- borderPad/pixel2um
    aFiles <- epFiles <- gFiles <- iFiles <- NULL
    fovAnn <- fovAnnotations %>% filter(CELL_DIVE_ID == cellDiveId)

    if(is.null(haloAnn)){
        ## filter fov annotations for this sample
        fovAnn$Path_to_exclusion_boundaries <- convertBoundaryFilePaths(fovAnn$Path_to_boundary_files,"ExclusionCoordinates")
        fovAnn$Path_to_epidermis_boundaries <- convertBoundaryFilePaths(fovAnn$Path_to_boundary_files,"EpidermisCoordinates")
        fovAnn$Path_to_glass_boundaries <- convertBoundaryFilePaths(fovAnn$Path_to_boundary_files,"GlassCoordinates")
        fovAnn$Path_to_interface_boundaries <- convertBoundaryFilePaths(fovAnn$Path_to_boundary_files,"InterfaceCoordinates")

        ## get halo exclusion files for this sample
        pths <- unique(fovAnn$Path_to_exclusion_boundaries)
        for(pth in pths[!is.na(pths)]){
            pat <- paste0("^",samp,"_")
            aFiles <- c(aFiles, file.path(pth,dir(pth)[grep(pat,dir(pth))]))
        }
        ## get epidermis annotation files for this sample
        pths <- unique(fovAnn$Path_to_epidermis_boundaries)
        for(pth in pths[!is.na(pths)]){
            pat <- paste0("^",samp,"_")
            epFiles <- c(epFiles, file.path(pth,dir(pth)[grep(pat,dir(pth))]))
        }
        ## get glass annotation files for this sample
        pths <- unique(fovAnn$Path_to_glass_boundaries)
        for(pth in pths[!is.na(pths)]){
            pat <- paste0("^",samp,"_")
            gFiles <- c(gFiles, file.path(pth,dir(pth)[grep(pat,dir(pth))]))
        }
    }

    ## get interface annotation files for this sample (for debugging only at this point)
    fovAnn$Path_to_interface_boundaries <- convertBoundaryFilePaths(fovAnn$Path_to_boundary_files,"InterfaceCoordinates")

    pths <- unique(fovAnn$Path_to_interface_boundaries)
    for(pth in pths[!is.na(pths)]){
        pat <- paste0("^",samp,"_")
        iFiles <- c(iFiles, file.path(pth,dir(pth)[grep(pat,dir(pth))]))
    }

    flog.debug("Converting min/max coordinates to midpoint coordinates")
    dd <- dat %>% mutate(X=(XMax+XMin)/2,Y=-(YMax+YMin)/2) ## need to keep X|Y min|max for drift loss 
    dd$EXCLUDE <- ""

    ## get exclusions determined by lab personel, found in FOV annotations meta data spreadsheet   
    studyExcl <- getStudyExclusions(fovAnn, dd)
    if(length(studyExcl) > 0){
        dd$EXCLUDE <- writeExclusionReasons(dd$EXCLUDE, studyExcl, "LabExclusion")
    } else {
        flog.info("  No study exclusions found.")
    }

    ## get drift/loss exclusions
    if(!is.null(drift)){
        driftExcl <- getDriftExclusions(drift, dd, driftThreshold)
        if(length(driftExcl) > 0){
            dd$EXCLUDE <- writeExclusionReasons(dd$EXCLUDE, driftExcl, paste0("DRIFT_",driftThreshold*100))
        } else {
            flog.info("  No drift exclusions found.")
        }
    } else {
        flog.info("  No drift exclusions marked.")
    }

    ## mark cells that fall inside Halo exclusion boundaries (exclusions, epidermis, glass) and also
    ## those that fall outside limits set by border pad 
    for(fov in unique(dd$SPOT)){ 
        aFile <- epFile <- gFile <- iFile <- NULL

        ## exclude points that fall within outside border
        flog.info(paste0("  getting exclusions for FOV ",fov))
        flog.info("    excluding points that fall within padding")
        borderExcl <- getBorderPaddingExclusions(fov, samp, dd, borderPad_um)
        if(length(borderExcl) > 0){
            dd$EXCLUDE <- writeExclusionReasons(dd$EXCLUDE, borderExcl, paste0("PADDED_",borderPad))
        } else {
            flog.info("      No border exclusions found.")
        }

        ## if parsed halo annotation is not given, get all halo annotations files for this FOV
        if(is.null(haloAnn)){
            sampAnns <- aFiles[grep(paste0(samp,"_"),aFiles)]
            aFile <- sampAnns[grep(paste0("Spot",fov,".annotations"),sampAnns)]
            sampAnns <- epFiles[grep(paste0(samp,"_"),epFiles)]
            epFile <- sampAnns[grep(paste0("Spot",fov,".annotations"),sampAnns)]
            sampAnns <- gFiles[grep(paste0(samp,"_"),gFiles)]
            gFile <- sampAnns[grep(paste0("Spot",fov,".annotations"),sampAnns)]
            sampAnns <- iFiles[grep(paste0(samp,"-"),iFiles)]
            iFile <- sampAnns[grep(paste0("Spot",fov,".annotations"),sampAnns)]
        }

        ## get halo exclusions 
        flog.info("    excluding points in Halo exclusion annotation files")
        haloExcl <- getHaloExclusions(fov, dd, "Exc", haloAnnotations=haloAnn, aFile=aFile, 
                                      boundaryColors=boundaryColors, boundaryReassignmentFile=boundaryReassignmentFile)
        if(length(haloExcl) > 0){
            dd$EXCLUDE <- writeExclusionReasons(dd$EXCLUDE, haloExcl, "HALOExclusion")
        } else {
            flog.info("      No Halo exclusion boundaries found.")
        }

        ## get halo epidermis exclusions 
        flog.info("    excluding points in Halo epidermis annotation files")
        haloEpi <- getHaloExclusions(fov, dd, "Epi", haloAnnotations=haloAnn, aFile=epFile, 
                                     boundaryColors=boundaryColors, boundaryReassignmentFile=boundaryReassignmentFile)
        if(length(haloEpi) > 0){
            dd$EXCLUDE <- writeExclusionReasons(dd$EXCLUDE, haloEpi, "HALOEpidermis")
        } else {
            flog.info("      No Halo epidermis boundaries found.")
        }

        ## get halo glass exclusions 
        flog.info("    excluding points in Halo glass annotation files")
        haloGls <- getHaloExclusions(fov, dd, "Gls", haloAnnotations=haloAnn, aFile=gFile, 
                                     boundaryColors=boundaryColors, boundaryReassignmentFile=boundaryReassignmentFile)
        if(length(haloGls) > 0){
            dd$EXCLUDE <- writeExclusionReasons(dd$EXCLUDE, haloGls, "HALOGlass")
        } else {
            flog.info("      No Halo glass boundaries found.")
        }
 
        ## debug with plots
        if(printPlots){
            flog.info("Printing plots for debugging")
            plotExclusions(dd[which(dd$SPOT==fov),],haloAnnotations=haloAnn,iFiles,aFiles,epFiles,gFiles,outDir=debugDir,
                           boundaryReassignmentFile=boundaryReassignmentFile)
        }
    }
    dd %>% select(-(X:Y))
}


