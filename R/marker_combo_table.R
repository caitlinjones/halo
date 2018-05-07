library(tidyverse)
library(xlsx)


#' Find all marker combinations existing in data and count cells that
#' express each combination and sort combinations by the most to least frequent
#' 
#' Each row contains cell count, percentage and cumulative percentage of all cells,
#' and a column for each marker, with an "X" indicating that marker is part of the
#' combination in that row
#' 
#' @param dataFiles   vector of *.rda files containing Halo data
#' @param markerFile  file containing a single identity marker on each line     
#' @param outFile     name of *.xlsx file to write final tables; default=NULL
#' @return list of tables, one for all samples and one for each sample
markerComboCounts <- function(dataFiles, markerFile, outFile=NULL){

    samps <- c()
    allSamps <- tibble()

    cellTypeMarkers <- scan(markerFile,"")

    ## get table of "Positive" values for all cell type markers
    ## rows are cells, columns are Sample,UUID,cellTypeMarkers
    for(df in dataFiles){

        dat <- readRDS(df)
        dat$Sample <- gsub("_ObjectAnalysisData","",dat$Sample)
        samps <- unique(c(samps, unique(dat$Sample)))

        pTbl <- dat %>%
                filter(Marker %in% cellTypeMarkers, ValueType=="Positive") %>%
                mutate(Marker=factor(Marker,levels=cellTypeMarkers)) %>%
                arrange(Marker) %>%
                dplyr::select(Sample,UUID, Marker, Value) %>%
                spread(Marker,Value)

        allSamps <- allSamps %>% bind_rows(pTbl)
    }

    allTbls <- list()

    ## get counts of each combination of 0's and 1's (representing
    ## a marker combination) & calculate total percentage of all cells
    ## and cumulative percentage
    for(s in c("AllSamples",samps)){

        if(s!="AllSamples"){
            pTbl <- allSamps %>% filter(Sample==s)
            append <- TRUE
        } else {
            pTbl <- allSamps
            append <- FALSE
        }
        cTbl <- pTbl %>% dplyr::select(-(Sample:UUID)) %>%
                group_by_at(vars(one_of(cellTypeMarkers))) %>%
                summarise(Count=n()) %>%
                arrange(desc(Count))
    
        ### these things can't be piped for some reason
        cTbl$PCT <- cTbl$Count/sum(cTbl$Count)
        cTbl$CUM.PCT <- cumsum(cTbl$PCT)

        ### reorder columns 
        cTbl <- cTbl %>% dplyr::select(Count, PCT, CUM.PCT, everything())

        cTbl[,4:ncol(cTbl)][cTbl[,4:ncol(cTbl)] == 0] <- ""
        cTbl[,4:ncol(cTbl)][cTbl[,4:ncol(cTbl)] == 1] <- "X"
 
        allTbls[[s]] <- cTbl
        if(!is.null(outFile)){
            write.xlsx(as.data.frame(cTbl), outFile, sheet=s, append=append, row.names=FALSE)
        }
    }

    return(allTbls)
}

#' Build a list of all marker combinations described in "complex" combination
#' format in cell types spreadsheet (described in docs)
#'
#' Based on combination type ("ALL", "ANY", "<2", ">3", etc.), generate all
#' possible combinations of markers given
#' 
#' @param comboType  tells what kind of combinations to build ("ALL" indicates
#'                   all markers given must stay together, "ANY" means any combination
#'                   of given markers, including all different lengths, and 
#'                   types like ">2" or "<=3" tells how many markers each combination must 
#'                   contain
#' @param markers    vector of markers to combine
#' @return a list of all possible combinations matching criteria in comboType
getMarkerCombos <- function(comboType, markers){
    #if(markers == "TOTAL"){
        ## handle differently
    #}

    markers <- unlist(strsplit(markers,","))
    if(comboType == "ALL"){
        nums <- length(markers)
    } else if(comboType == "ANY") {
        nums <- c(1:length(markers))
    } else if( length(grep(">=",comboType)) > 0 ){
        num <- as.integer(gsub(">=","",comboType))
        nums <- c(num:length(markers))
    } else if( length(grep("<=",comboType)) > 0 ){
        num <- as.integer(gsub("<=","",comboType))
        nums <- c(1:num) 
    } else if( length(grep(">",comboType)) > 0 ){
        num <- as.integer(gsub(">","",comboType))
        nums <- c((num+1):length(markers))
    } else if( length(grep("<",comboType)) > 0 ){
        num <- as.integer(gsub("<","",comboType))
        nums <- c(1:(num-1)) 
    }
    mSets <- list()
    for(i in nums){
        cmb <- t(combn(sort(markers),m=i))
        for(c in 1:nrow(cmb)){
            mSets[[length(mSets)+1]] <- cmb[c,]
        }
    }
    return(mSets)
}



getCellTypes <- function(cellTypesFile, unreasonableCombinationsFile){

    simpleTypes <- read.xlsx(cellTypesFile,sheetName="Simple",stringsAsFactors=FALSE)
    complexTypes <- read.xlsx(cellTypesFile,sheetName="Complex",stringsAsFactors=FALSE)
    unCombos <- read.xlsx(unreasonableCombinationsFile,sheetIndex=1)
    conditionalTypes <- NULL

    cellTypes <- simpleTypes

    ## sort combos from simple sheet file
    for(x in 1:nrow(cellTypes)){
        combo <- unlist(strsplit(cellTypes$Marker_combination[x],","))
        cellTypes$Marker_combination[x] <- paste(sort(combo),collapse=",")
    }

    ## parse complex types into "simple" types by generating all possible
    ## combinations for the type given
    for(x in 1:nrow(complexTypes)){
        ct <- complexTypes[x,]
        ## pull out conditional types to be parsed later
        if("TOTAL" %in% c(ct$OF,ct$OF.1)){
            conditionalTypes <- rbind(conditionalTypes,ct)
            next
        } 
        set1 <- getMarkerCombos(ct$IF,ct$OF)
        set2 <- getMarkerCombos(ct$MATCHES,ct$OF.1)    
        allCombos <- expand.grid(set1,set2)
        finalCombos <- c()
        for(i in 1:length(allCombos[[1]])){
            combo <- unique(sort(unlist(c(allCombos[[1]][i],allCombos[[2]][i]))))
            finalCombos <- c(finalCombos,paste(combo,collapse=","))
        }
        finalCombos <- unique(finalCombos)
        overlaps <- finalCombos[which(finalCombos %in% cellTypes$Marker_combination)]
        if(length(overlaps) > 0){
            stop(paste0("\n\nThe following 'Complex' marker combinations also occur in the list of 'Simple' combinations: ",paste(overlaps,collapse="; "),". Please correct and rerun.\n\n"))
        }
        cts <- data.frame(Marker_combination = finalCombos, 
                          Cell_type = rep(ct$CELL_TYPE,length(finalCombos)),
                          Subtype = rep(ct$CELL_TYPE,length(finalCombos)))
        cellTypes <- rbind(cellTypes,cts) 
    }
 
    ## parse unreasonable combos matrix to get the "Unknown" combinations
    rownames(unCombos) <- unCombos[,1]
    unCombos <- unCombos[,-1]    
    allUnknowns <- c()
    for(x in rownames(unCombos)){
        unknowns <- colnames(unCombos)[which(unCombos[x,] == "Unknown")]
        if(!is.null(unknowns) & length(unknowns) > 0){
            combos <- unlist(lapply(unknowns,function(y){ paste(sort(c(y,x)),collapse=",") })) 
            overlaps <- combos[which(combos %in% cellTypes$Marker_combination)]
            if(length(overlaps) > 0){
                stop(paste0("\n\nThe following 'Unknown' marker combinations also occur in the CellTypes xlsx file: ",paste(overlaps,collapse="; "),". Please correct and rerun.\n\n"))
            }
            allUnknowns <- c(allUnknowns,combos)
        }
    }
    uc <- data.frame(Marker_combination = sort(allUnknowns),
                     Cell_type = rep("Unknown", length(allUnknowns)),
                     Subtype = rep("Unknown", length(allUnknowns)))
    cellTypes <- rbind(cellTypes,uc)   

    return(list(CellTypes=cellTypes,Conditional=conditionalTypes))
}

assignCellType <- function(combinationRow,cellTypes,conditionalTypes){
    markers <- names(combinationRow)[which(combinationRow == "X")]
    combo <- paste(sort(markers),collapse=",")
    if(combo %in% cellTypes$Marker_combination){
        cellType <- cellTypes$Cell_type[which(cellTypes$Marker_combination==combo)]
    } else if(is.null(markers) | length(markers) == 0){
        cellType <- "NEGATIVE"
    } else if(combo %in% conditionalTypes$Marker_combination){
        cellType <- "NOT SURE YET"   
    } else {
        ## all combos not in cell types file are considered unreasonable for now
        cellType <- "Unreasonable"
    }
    return(cellType)
}

#' Interpret marker combinations
#' 
#' Given a file of cell types defined by certain marker combinations (format
#' in docs), assign cell types to each combination in a marker combination counts
#' file
#' 
#' @param markerComboCounts       table generated by markerComboCounts()
#' @param cellTypesFile           *.xlsx file defining cell types based on known marker
#'                                combination (format in docs)
#' @param unreasonableCombosFile  *.xlsx file containing a matrix indicating both known 
#'                                unreasonable pairwise marker combinations and also
#'                                "unknown" pairwise marker combinations
#' @param outFile                 if given, markerComboCounts table with interpretations
#'                                added will be printed to this file
interpretMarkerCombos <- function(markerComboCts, cellTypesFile,
                                  unreasonableCombosFile, outFile=NULL){

    allCellTypes <- getCellTypes(cellTypesFile, unreasonableCombosFile)
    
    cellTypes <- allCellTypes$CellTypes
    conditionalTypes <- allCellTypes$Conditional

    ctCounts <- list()
    cellTypeList <- c()
    cellTypeNums <- c()
    comboStrings <- c()
    for(rowNum in 1:nrow(markerComboCts)){
        row <- markerComboCts[rowNum,]
        ct <- assignCellType(row, cellTypes, conditionalTypes)
        if(!ct %in% names(ctCounts)){ 
            ctCounts[[ct]] <- 0
        }
        ctCounts[[ct]] <- ctCounts[[ct]] + 1
        markers <- names(row)[which(row == "X")]
        if(length(markers) == 0){
            comboStr <- ""
        } else if(length(markers) == 1){ 
            comboStr <- paste0(markers, " only") 
        } else {
            comboStr <- paste(markers, collapse = "/")
        }
        cellTypeList <- c(cellTypeList,ct)
        cellTypeNums <- c(cellTypeNums, ifelse(ct == "NEGATIVE","",paste0("#",ctCounts[[ct]])))
        comboStrings <- c(comboStrings, comboStr)
    }
    
    ## add a separate column for each piece, to be pasted together later
    markerComboCts$cellType <- cellTypeList
    markerComboCts$cellTypeNum <- cellTypeNums
    markerComboCts$comboString <- comboStrings
    
    return(markerComboCts)
}

#' Create XLSX workbook and style all sheets for marker combo tables
#'
#' For each marker combo table, create a XLSX sheet and style
#' according to original specifications (TO DO: add options?)
#'
#' @param allTbls  list of tables, one for each sheet, each containing
#'                 the following columns: Count, PCT, CUM.PCT, Interpretation, 
#'                 and a column for each marker being counted
#' @return a completely styled XLSX workbook
markerComboWorkbook <- function(allTbls, cellTypes){

    ## tmp
    cellTypes <- sort(c("Tumor","Macrophage","NEGATIVE","UNKNOWN",
                      "Unreasonable","T cell","B cell", "Natural killer cell"))
    ## alternating background colors for different cell types
    cellTypeColors <- rep(9,length(cellTypes)) # white
    cellTypeColors[which(1:length(cellTypes) %% 2 == 0)] <- 55 #grey
    
    wb <- createWorkbook(type="xlsx")
    MAIN_STYLE <- CellStyle(wb, border=NULL) + 
                  Alignment(horizontal="ALIGN_CENTER")
    HEADER_STYLE <- MAIN_STYLE + 
                    Font(wb, isBold=TRUE) + 
                    Fill(backgroundColor="#A9A9A9") + 
                    Border(color="black", position="BOTTOM", pen="BORDER_THICK")
    INTERP_STYLE <- MAIN_STYLE + 
                          Font(wb, color="#FF0000", isBold=TRUE)
    INTERP_STYLE_5PCT <- MAIN_STYLE + 
                          Font(wb, color="#6bccf9", isBold=TRUE)
    COUNT_STYLE <- MAIN_STYLE + 
                   Font(wb, color="#FF0000", isBold=TRUE)
    TOTALS_BORDER <- Border(color="black", position=c("BOTTOM"))

    for(t in 1:length(allTbls)){
        s <- names(allTbls)[t]
print(s)
        tbl <- allTbls[[s]]
        sheet <- createSheet(wb, sheetName = s)
        interpCol <- which(names(tbl)=="Interpretation") 
        countCol <- which(names(tbl)=="Count")
        totalsRows <- which(is.null(tbl$CUM.PCT)|nchar(tbl$CUM.PCT)==0)
        totalsRows <- totalsRows[-length(totalsRows)]
        columnStyles <- rep(list(MAIN_STYLE), ncol(tbl))
        names(columnStyles) <- 1:ncol(tbl)
        columnStyles[[interpCol]] <- INTERP_STYLE
        columnStyles[[countCol]] <- COUNT_STYLE
        ## add data
        addDataFrame(as.data.frame(tbl), sheet, row.names=FALSE, startRow=1, startColumn=1,
                     colnamesStyle=HEADER_STYLE, colStyle=columnStyles)
        ## fit interpretation column
        autoSizeColumn(sheet, interpCol)

        ## put border around totals columns (really just for sheets that are by cell type)
        for(tr in totalsRows){
            ## add one to row number bc it looks like once in a sheet, the header
            ## is no longer row 1? (or rows are indexed by 0?)
            cb <- CellBlock(sheet, tr+1, 1, 1, ncol(tbl), create=FALSE)
            CB.setBorder(cb, TOTALS_BORDER, 1, 1:ncol(tbl)) 
        }   

        if(grepl("by_cell_type",s)){
            for(x in seq(cellTypes)){
                ctRows <- grep(cellTypes[x],tbl$Interpretation)
                if(length(ctRows) == 0){ next }

                rowIdxs <- c()
                for(r in ctRows){ rowIdxs <- c(rowIdxs,rep(r+1,ncol(tbl))) }

                colIdxs <- rep(1:ncol(tbl), length(ctRows))

                fill <- Fill(backgroundColor=cellTypeColors[x])
                cb <- CellBlock(sheet, ctRows[1]+1, 1, length(ctRows), ncol(tbl), create=FALSE)
                CB.setFill(cb, fill, rowIdxs, colIdxs)
            }
        }
    }

    return(wb)
}

#' Generate XLSX file of all marker combinations represented in
#' given halo data files
#' 
#' XLSX file generated will contain a detailed count sheet and a summary sheet for ALL data, 
#' incluing all samples, and for each sample. Each sheet includes marker count, percentage 
#' and cumulative percentage of all cells that contain each marker combination, and interpretation
#' or cell type of each marker combination
#' 
#' @param dataFiles               vector of *.rda files containing Halo data
#' @param markerFile              file containing a single identity marker on each line     
#' @param cellTypesFile           a *.xlsx file containing sheets for Simple and Complex
#'                                marker combinations that identify certain cell types (format
#'                                described in docs)
#' @param unreasonableCombosFile  a *.xlsx file containing a matrix that indicated invalid or
#'                                "unknown" pairwise combinations of markers in markerFile (format
#'                                described in docs)
#' @param outFile                 name of *.xlsx file to write final tables
#' @export
markerComboXLSX <- function(dataFiles, markerFile, cellTypesFile, unreasonableCombosFile,outFile){
    allCountTbls <- markerComboCounts(dataFiles, markerFile)

    for(x in 1:length(allCountTbls)){
        s <- names(allCountTbls)[x]
        tbl <- allCountTbls[[s]]
        iTbl <- interpretMarkerCombos(tbl, cellTypesFile, unreasonableCombosFile)

        allCountTbls[[s]] <- iTbl %>% 
                unite("Interpretation", cellType, cellTypeNum, comboString, sep = " ") %>%
                dplyr::select(Count, PCT, CUM.PCT, Interpretation, everything())

        totalsByType <- as.tibble(iTbl) %>%
                        filter(CUM.PCT <= 0.95) %>% 
                        group_by(cellType) %>%
                        summarize(Count=sum(Count),PCT=sum(PCT))

        rep_na <- as.list(rep("",ncol(iTbl)))
        names(rep_na) <- names(iTbl)
        allCountTbls[[paste0(s,"_by_cell_type")]] <- iTbl %>% 
                        bind_rows(totalsByType) %>% 
                        replace_na(rep_na) %>%
                        arrange(cellType) %>%
                        unite("Interpretation", cellType, cellTypeNum, comboString, sep = " ") %>%
                        dplyr::select(Count, PCT, CUM.PCT, Interpretation, everything())

    }

    mcWorkbook <- markerComboWorkbook(allCountTbls)

    saveWorkbook(mcWorkbook, "testing_new.xlsx")

    return(allCountTbls)
}
