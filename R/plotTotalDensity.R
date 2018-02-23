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
