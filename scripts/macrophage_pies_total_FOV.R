 plotTheme <- theme(legend.title = element_blank(),
                  legend.position = "right",
                  legend.text = element_text(size=12),
                  axis.text = element_text(size=12),
                  axis.text.x = element_text(size=8),
                  axis.title.y = element_text(size=14),
                  axis.title.x = element_text(size=14),
                  axis.ticks = element_blank(),
                  axis.line = element_line(color = "black", size=0.75),
                  strip.text.x = element_text(size=12, margin=margin(.1, .1, .1, .1, "cm")),
                  strip.placement = "outside",
                  plot.background = element_blank(),
                  plot.title = element_text(size=20, margin=margin(b=0.1, t=0.2, r=0, l=0, "cm")),
                  plot.margin = unit(c(0,0,0,0), "cm"), 
                  panel.background = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()
                  )

outDir <- file.path(pp$fov_density_dir,"macrophage_pies")
sampleOrder <- names(sort(sapply(unique(fovDensity$Sample),function(x){ as.numeric(unlist(strsplit(x,"_"))[2]) })))
overviewTitle <- "Macrophage population cell types"
samplePlotTitle <- overviewTitle

macrophage_pies <- function(fovDensity, markerConfig, sampleOrder=NULL, outDir=NULL, overviewTitle=NULL, samplePlotTitle=NULL, functionalPlotTitle=NULL){

    allPlots <- list()

    ## get config for cell types that do NOT included a functional marker
    oviewCfg <- markerConfig %>% filter(CellType == Population)
    ## set up labels and colors based on config
    ctLabels <- oviewCfg$Label
    names(ctLabels) <- oviewCfg$CellType
    clrLbls <- unique(oviewCfg[,c("Color","Label")])
    ctClrs <- pull(clrLbls[,"Color"])
    names(ctClrs) <- pull(clrLbls[,"Label"])
    cellTypeLabels <- ctLabels
    clrs <- ctClrs



    allCounts <- fovDensity %>%
                 filter(CellType %in% unique(oviewCfg$Population)) %>%
                 select(Sample,SPOT,CellType,Counts)

    totalFOVcounts <- allCounts %>%
                      group_by(Sample,SPOT) %>%
                      summarise(TotalFOVCounts=sum(Counts))
    allCounts <- allCounts %>%
                 left_join(totalFOVcounts, by=c("Sample","SPOT")) %>%
                 mutate(PctFOV=Counts/TotalFOVCounts)

    ctCountsBySample <- allCounts %>%
                        group_by(Sample,CellType) %>%
                        summarise(TotalSampleCellTypeCounts=sum(Counts))

    allCounts <- allCounts %>%
                 left_join(ctCountsBySample, by=c("Sample","CellType"))

    totalSampleCounts <- allCounts %>%
                         group_by(Sample) %>%
                         summarise(TotalSampleCounts=sum(Counts))

    allCounts <- allCounts %>%
                 left_join(totalSampleCounts, by=c("Sample")) %>%
                 mutate(PctSample=TotalSampleCellTypeCounts/TotalSampleCounts)

    flog.debug("setting up labels")
    allCounts$CellTypeLabels <- unlist(allCounts$CellType)
    if(!is.null(cellTypeLabels)){
        for(x in names(cellTypeLabels)){
            allCounts$CellTypeLabels[which(allCounts$CellType == x)] <- cellTypeLabels[[x]]
        }
    }
    allCounts$CellTypeLabels <- factor(allCounts$CellTypeLabels, levels=cellTypeLabels)
    allCounts[is.na(allCounts)] <- 0

    if(!is.null(sampleOrder)){
        allCounts$Sample <- factor(allCounts$Sample, levels=sampleOrder)
    }


    pdf(file.path(outDir, "sample_population_overview.pdf"), height=8.5, width=11, onefile=TRUE)
    popOverview <- allCounts %>% select(Sample,PctSample,CellTypeLabels) %>% unique()
    p1 <- ggplot(popOverview, aes(x="", y=PctSample, fill=CellTypeLabels)) +
                  geom_bar(stat="identity") +
                  scale_fill_manual(values=clrs, labels=names(clrs)) +
                  coord_polar("y") +
                  facet_wrap(~Sample) +
                  #getPlotTheme("infiltration density") + 
                  plotTheme +
                  theme(axis.text.x = element_blank(),
                        axis.line = element_blank(),
                        axis.title.x = element_blank(),
                        axis.title.y = element_blank()) +
                  labs(title=paste0(overviewTitle, " by sample"))
    print(p1)
    dev.off()

    allPlots[["overview"]] <- p1



    funcConfig <- parsedMarkerConfig %>% filter(CellType != Population)
    funcDat <- fovDensity %>%
               filter(CellType %in% funcConfig$CellType) %>%
               select(Sample,SPOT,CellType,Counts) %>%
               group_by(Sample,CellType) %>% 
               summarize(SampleCellTypeCounts=sum(Counts))

    funcDat$Functional <- unlist(lapply(strsplit(funcDat$CellType,","),
                          function(x){ gsub("-","",x[length(x)]) }))
    funcDat$Population <- unlist(lapply(strsplit(funcDat$CellType,","),
                          function(x){ paste0(x[-length(x)],collapse=",") }))
    funcDat$PopulationLabel <- unlist(lapply(funcDat$Population, function(x){ removeNegativeMarkers(x) }))


    ## get totals of each functional marker (pos and neg)
    funcSum <- funcDat %>%
               group_by(Sample,Population,Functional) %>%
               summarize(SamplePopFunctionalCount=sum(SampleCellTypeCounts))

    funcDat <- funcDat %>%
               left_join(funcSum, by=c("Sample","Population","Functional")) %>%
               mutate(Pct=SampleCellTypeCounts/SamplePopFunctionalCount)
               funcDat[is.na(funcDat)] <- 0

    cellTypeLabels <- markerConfig$Label
    names(cellTypeLabels) <- markerConfig$CellType
    cellTypeLabels <- lapply(cellTypeLabels, function(x){ removeNegativeMarkers(x) })
    cellTypeLabels[which(cellTypeLabels == "Negative")] <- markerConfig$Population[which(cellTypeLabels == "Negative")]
    cellTypeLabels <- lapply(cellTypeLabels, function(x){ removeNegativeMarkers(x) })

    clrLbls <- unique(markerConfig[,c("Color","Label")])
    clrs <- pull(clrLbls[,"Color"])
    names(clrs) <- pull(clrLbls[,"Label"])
    clrs <- clrs[!grepl("-",names(clrs))]
    clrs[which(clrs == "gray")] <- "#914bf6"
    clrs[["TIM3"]] <- "#914bf6"
    clrs[["PD1"]] <- "#914bf6"
    clrs[["LAG3"]] <- "#914bf6"


    funcDat$CellTypeLabels <- unlist(funcDat$CellType)
    if(!is.null(cellTypeLabels)){
        for(x in names(cellTypeLabels)){
            funcDat$CellTypeLabels[which(funcDat$CellType == x)] <- cellTypeLabels[[x]]
        }
    }
    funcDat$CellTypeLabels <- factor(funcDat$CellTypeLabels, levels=unique(cellTypeLabels))
    funcDat[is.na(funcDat)] <- 0

    funcSummarykeyLabels <- clrs[-which(clrs == "#914bf6")]
    funcDat$PopulationLabel <- factor(funcDat$PopulationLabel, levels=names(clrs))

    if(!is.null(sampleOrder)){
        funcDat$Sample <- factor(funcDat$Sample, levels=sampleOrder)
    }


    for(s in unique(sampleOrder)){

        print(s)

        ## get config for cell types that do NOT included a functional marker
        oviewCfg <- markerConfig %>% filter(CellType == Population)
        ## set up labels and colors based on config
        cellTypeLabels <- oviewCfg$Label
        names(cellTypeLabels) <- oviewCfg$CellType
        clrLbls <- unique(oviewCfg[,c("Color","Label")])
        clrs <- pull(clrLbls[,"Color"])
        names(ctClrs) <- pull(clrLbls[,"Label"])

        allPlots[[s]] <- list()

        pdf(file.path(outDir, paste0(s,"_macrophages.pdf")), height=8.5, width=11, onefile=TRUE)

        sDat <- allCounts %>% filter(Sample == s)

        ### sample summary (all FOV combined)
        sSum <- sDat %>% select(Sample,CellTypeLabels,PctSample) %>% unique()
        keyLabels <- paste0(sSum$CellTypeLabels, " (",round(sSum$PctSample*100),"%)")
        names(keyLabels)=sSum$CellTypeLabels

        p1 <- ggplot(sSum, aes(x="", y=PctSample, fill=CellTypeLabels)) +
                      geom_bar(stat="identity") +
                      scale_fill_manual(values=clrs, labels=keyLabels) +
                      coord_polar("y") +
                      plotTheme +
                      theme(axis.text.x = element_blank(),
                            axis.line = element_blank(),
                            axis.title.x = element_blank(),
                            axis.title.y = element_blank()) +
                      labs(title=paste0(overviewTitle, "\n", s))
        print(p1)
        allPlots[[s]][["summary"]] <- p1

        #### by FOV
        p1 <- ggplot(sDat, aes(x="", y=PctFOV, fill=CellTypeLabels)) +
                      geom_bar(stat="identity") +
                      scale_fill_manual(values=clrs, labels=names(clrs)) +
                      coord_polar("y") +
                      facet_wrap(~SPOT) +
                      plotTheme +
                      theme(axis.text.x = element_blank(),
                            axis.line = element_blank(),
                            axis.title.x = element_blank(),
                            axis.title.y = element_blank()) +
                      labs(title=paste0(overviewTitle, "\n", s, " by FOV"))
        print(p1)
        allPlots[[s]][["by_fov"]] <- p1
  
  
        for(pop in unique(sDat$CellType)){

            print(pop)

            popCfg <- funcConfig %>% filter(Population == pop)
            ctLabels <- popCfg$Label
            names(ctLabels) <- popCfg$CellType
            cellTypeLabels <- ctLabels

            popCfg$Color <- "gray"
            negIdx <- which(grepl("-",popCfg$Label))
            popClr <- unique(oviewCfg$Color[which(oviewCfg$CellType == pop)])
            clrs <- c(popClr, "#914bf6")
            names(clrs) <- c("-","+")
            popCfg$Color[negIdx] <- popClr
 
            popCfg <- popCfg %>% filter(CellType != pop)
        
            popDat <- funcDat %>%
                      filter(Sample == s, CellType %in% popCfg$CellType) %>%
                      mutate(Functional=gsub("-","",gsub(paste0(pop,","),"",CellType)))
     
            #popDat$Functional <- gsub("-","",gsub(paste0(pop,","),"",popDat$CellType))

            popDatSum <- popDat %>%
                         group_by(Sample,Functional) %>%
                         summarize(TotalSamplePopulationCount=sum(SampleCellTypeCounts))

            popDat <- popDat %>% 
                     left_join(popDatSum, by=c("Sample","Functional")) %>%
                     mutate(Pct=SampleCellTypeCounts/TotalSamplePopulationCount)

            popDat[is.na(popDat)] <- 0

            popDat$CellTypeLabels <- unlist(popDat$CellType)
            if(!is.null(cellTypeLabels)){
                for(x in names(cellTypeLabels)){
                    popDat$CellTypeLabels[which(popDat$CellType == x)] <- cellTypeLabels[[x]]
                }
            }

            popDat$CellTypeLabels[which(grepl("-",popDat$CellTypeLabels))] <- "-"
            popDat$CellTypeLabels[-which(grepl("-",popDat$CellTypeLabels))] <- "+"
            popDat$CellTypeLabels <- factor(popDat$CellTypeLabels, levels=names(clrs))

            popDat$Functional <- gsub("-","",gsub(paste0(pop,","),"",popDat$CellType))
            popDat$Label <- ""
            pos <- which(popDat$CellTypeLabels == "+")
            popDat$Label[pos] <- paste0(round(popDat$Pct[pos]*100),"%")

            titles <- popDat %>%
                     filter(Label != "") %>%
                     select(Sample,Functional, Label) %>%
                     mutate(Title=paste0(Functional," (",Label,")")) %>%
                     select(-(Label))

            popDat <- popDat %>%
                     left_join(titles, by=c("Sample","Functional"))

            p1 <- ggplot(popDat, aes(x="", y=Pct, fill=CellTypeLabels)) +
                  geom_bar(stat="identity") +
                  facet_wrap(~Title) +
                  scale_fill_manual(values=clrs, labels=names(clrs)) +
                  coord_polar("y") +
                  plotTheme +
                  theme(axis.text.x = element_blank(),
                        axis.line = element_blank(),
                        axis.title.x = element_blank(),
                        axis.title.y = element_blank()) +
                  labs(title=s, subtitle=pop)
            print(p1)

            allPlots[[s]][[pop]] <- p1
        }


        #### full sample summary of functional markers in populations
        ## get config for cell types that do NOT included a functional marker
        ## set up labels and colors based on config
        sDat <- funcDat %>% filter(Sample == s)
        clrLbls <- unique(markerConfig[,c("Color","Label")])
        clrs <- pull(clrLbls[,"Color"])
        names(clrs) <- pull(clrLbls[,"Label"])
        clrs <- clrs[!grepl("-",names(clrs))]
        clrs[which(clrs == "gray")] <- "#914bf6"
        clrs[["TIM3"]] <- "#914bf6"
        clrs[["PD1"]] <- "#914bf6"
        clrs[["LAG3"]] <- "#914bf6"


        p1 <- ggplot(sDat, aes(x="", y=Pct, fill=CellTypeLabels)) +
             geom_bar(stat="identity") +
             facet_grid(Functional ~ PopulationLabel, switch="both") +
             scale_fill_manual(values=clrs, labels=names(clrs)) +
             coord_polar("y") +
             plotTheme +
             theme(axis.text.x = element_blank(),
                   axis.line = element_blank(),
                   axis.title.x = element_blank(),
                   axis.title.y = element_blank(),
                   strip.text.x = element_text(size = 10, angle=90, hjust=1),
                   strip.text.y = element_text(size = 10, angle=180),
                   strip.background = element_blank(),
                   plot.margin = unit(c(0,0,0,0), "cm"),
                   legend.position="none",
                   panel.spacing = unit(0.1, "cm")) +
             labs(title=s)
        print(p1)

        dev.off()

    }

    mdir <- file.path(outDir,"by_marker")

    allPlots[["by_marker"]] <- list()
    ####### summarize by marker
    for(m in unique(funcDat$Functional)){
        print(m)
        pdf(file.path(mdir,paste0(m,".pdf")),height=8.5,width=11,onefile=TRUE)
        tst <- funcDat %>% filter(Functional == m)
        p1 <- ggplot(tst, aes(x="", y=Pct, fill=CellTypeLabels)) +
             geom_bar(stat="identity") +
             #facet_wrap(~Title) +
             facet_grid(Sample ~ PopulationLabel, switch="both") +
             scale_fill_manual(values=clrs, labels=names(clrs)) +
             coord_polar("y") +
             plotTheme +
             #getPlotTheme("infiltration density") + 
             theme(axis.text.x = element_blank(),
                   axis.line = element_blank(),
                   axis.title.x = element_blank(),
                   axis.title.y = element_blank(),
                   strip.text.x = element_text(size = 10, angle=90, hjust=1),
                   strip.text.y = element_text(size = 10, angle=180),
                   strip.background = element_blank(),
                   plot.margin = unit(c(0,0,0,0), "cm"),
                   legend.position="none",
                   panel.spacing = unit(0.1, "cm")) +
             labs(title=m)
        print(p1)
        dev.off() 
        allPlots[["by_marker"]][[m]] <- p1
   
    }

    return(allPlots)
} # end function



macrophage_pies(fovDensity, markerConfig, sampleOrder=sampleOrder, outDir=outDir, overviewTitle=overviewTitle, samplePlotTitle=samplePlotTitle, functionalPlotTitle=NULL)

