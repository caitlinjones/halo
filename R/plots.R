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
#' @param   v                   verbose - when set to TRUE, messages
#'                              will be printed to screen as well as to file;
#'                              DEFAULT = TRUE
#' @param   logFile             file to write log messages to; DEFAULT=NULL
#' @param   debug               print debug messages; DEFAULT=FALSE
#' @export
plotMarkerPercentages <- function(countsTable, markerNames, type="bar", other_threshold=0.05, 
                                    exclude_sample_fov=NULL, pdfFile=NULL, v=TRUE, 
                                    logFile=NULL, debug=FALSE,custom_colors=FALSE){

    ## set up colors
    markerColors <- c()

    if(!is.null(pdfFile)){
        pdf(pdfFile,width=10,height=8)#paper="a4r")
    }

    ## remove median rows from counts table
    countsTable <- filter(countsTable, !is.na(FOV) & !is.na(SLICE))

    ## remove any exclusions
    if(!is.null(exclude_sample_fov)){
        countsTable <- removeExclusions(countsTable, exclude_sample_fov, v=v, debug=debug)
    }

    allPlots <- c()
    for(samp in unique(countsTable$Sample)){
        if(debug){ logMsg(samp,v,logFile,"DEBUG") }

        pieTbl <- tibble()
        for(m in markerNames){
            if(debug){ logMsg(m,v,logFile,"DEBUG") }
            if(length(grep(m,names(countsTable))) == 0){
                if(debug){ logMsg(paste0("skipping ",m),v,logFile,"DEBUG") }
                next
            }

            tmp <- countsTable
            ## remove marker combos in counts table that do not have this POSITIVE marker
            posMarkerCombos <- names(tmp)[!grepl(paste0(m,"-"),names(tmp))]
            if(length(posMarkerCombos) == 3){
                logMsg(paste0("No ",m,"+ cells. Skipping."),v,logFile)
                next
            }
            ## get total marker counts and percentage of each marker for this sample
            tmp <- select(tmp, posMarkerCombos) %>%
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
                     select(-PercentAllSampleCells)  %>%
                     spread(Marker,TotalMarkerCount)
            if(ncol(all) > 1){
                all$Marker <- paste0(m," : ",rowSums(all[1,2:ncol(all)]))
            } else {
                all$Marker <- paste0(m," : ",all$Other)
            } 
            all <- select(all, Marker, everything())
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

        logMsg(paste0("Sample ", samp, " has ",length(unique(pie_data$Combo))," marker combinations"),v,logFile) 

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
        ## and readable; then print again with legend

        if(length(unique(pie_data$Combo)) > 70){
            ## plot WITHOUT legend and print
            p <- ggplot(pie_data, aes(x = factor(1), y = Count, fill=Combo)) +
                   facet_wrap(~Marker,nrow = numRows) +
                   geom_bar(width=1,stat="identity",position="fill") +
                   labs(x = "", y = "", title = paste0(samp,"\n")) +
                   scale_fill_manual(values = sampleColors) + 
                   theme(axis.text.x = element_blank(), 
         #            legend.text = element_text(size=6), 
         #            legend.key.size = unit(0.3,"cm"),
                     legend.position = "none",
                     strip.text.x = element_text(size=6)) +
                   guides(fill=guide_legend(title="",nrow=length(unique(pie_data$Combo))/4))
            if(type == "pie"){
                p <- p + coord_polar("y")
            }
        print(p)
        }

        ## add legend and print
    #    p <- ggplot(pie_data, aes(x = factor(1), y = Count, fill=Combo)) +
    #           facet_wrap(~Marker,nrow = numRows) +
    #           geom_bar(width=1,stat="identity",position="fill") +
    #           labs(x = "", y = "", title = paste0(samp,"\n")) +
    #           scale_fill_manual(values = sampleColors) +
    #           theme(axis.text.x = element_blank(),
    #                 legend.text = element_text(size=6),
    #                 legend.key.size = unit(0.3,"cm"),
    #                 legend.position = "right",
    #                 strip.text.x = element_text(size=6)) +
    #           guides(fill=guide_legend(title="",nrow=length(unique(pie_data$Combo))/4))
    #    if(type == "pie"){
    #        p <- p + coord_polar("y")
    #    }
        p <- p + theme(legend.text = element_text(size=6),
                       legend.key.size = unit(0.3,"cm"),
                       legend.position = "right")
        print(p)

    }

    if(!is.null(pdfFile)){
        dev.off()
    }
}


