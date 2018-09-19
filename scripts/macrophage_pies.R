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

outDir <- file.path(pp$infiltration_density_dir,"macrophage_pies")

macrophage_pies <- function(infiltrationDensity, infiltrationAreas, markerConfig, sampleOrder=NULL, outDir=NULL, overviewTitle=NULL, samplePlotTitle=NULL, functionalPlotTitle=NULL){

    allPlots <- list()

    ## get config for cell types that do NOT included a functional marker
    oviewCfg <- markerConfig %>% filter(CellType == Population)
    ## set up labels and colors based on config
    ctLabels <- oviewCfg$Label
    names(ctLabels) <- oviewCfg$CellType
    clrLbls <- unique(cfg[,c("Color","Label")])
    ctClrs <- pull(clrLbls[,"Color"])
    names(ctClrs) <- pull(clrLbls[,"Label"])
    cellTypeLabels <- ctLabels
    clrs <- ctClrs

    ###
    ### plot overview: percentages of all macrophage types in each sample
    ###

    #sia <- infiltrationAreas %>% 
           gather(3:ncol(infiltrationAreas), key="Band", value="Area") %>% 
           group_by(Sample) %>%
           summarise(TotalSampleArea=sum(Area)) %>% 
           ungroup()

    den <- infiltrationDensity %>% 
           filter(CellType %in% unique(oviewCfg$Population)) %>%
           select(Sample,SPOT,Band,CellType,Counts) %>%
           group_by(Sample,CellType) %>%
           summarise(TotalCounts=sum(Counts))

    totCts <- den %>% group_by(Sample) %>% summarise(TotalSampleCounts=sum(TotalCounts))

    den <- den %>% 
           left_join(sia, by=c("Sample")) %>% 
           mutate(Density=TotalCounts/TotalSampleArea) %>% 
           left_join(totCts, by=c("Sample")) %>% 
           mutate(Pct=TotalCounts/TotalSampleCounts)

    flog.debug("setting up labels")
    den$CellTypeLabels <- unlist(den$CellType)
    if(!is.null(cellTypeLabels)){
        for(x in names(cellTypeLabels)){
            den$CellTypeLabels[which(den$CellType == x)] <- cellTypeLabels[[x]]
        }
    }
    den$CellTypeLabels <- factor(den$CellTypeLabels, levels=cellTypeLabels)
    den[is.na(den)] <- 0

    if(!is.null(sampleOrder)){
        den$Sample <- factor(den$Sample, levels=sampleOrder)
    }


    pdf(file.path(outDir, "sample_population_overview.pdf"), height=8.5, width=11)
    p1 <- ggplot(den, aes(x="", y=Pct, fill=CellTypeLabels)) + 
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
                  labs(title=overviewTitle)
    print(p1)
    dev.off()

    allPlots[["overview"]] <- p1
    #########################
    ####### print by sample
    #########################
    for(s in sampleOrder){
        print(s)
        allPlots[[s]] <- list()

        ctLabels <- oviewCfg$Label
        names(ctLabels) <- oviewCfg$CellType
        clrLbls <- unique(cfg[,c("Color","Label")])
        clrs <- pull(clrLbls[,"Color"])
        names(clrs) <- pull(clrLbls[,"Label"])
        cellTypeLabels <- ctLabels

        sDir <- file.path(outDir,s)
        dir.create(sDir, showWarnings=FALSE, recursive=T)
        sDen <- infiltrationDensity %>% 
                  filter(Sample == s, CellType %in% oviewCfg$CellType) %>%
                  select(Sample,SPOT,Band,CellType,Counts) %>%
                  group_by(Sample,CellType) %>%
                  summarise(TotalCounts=sum(Counts))
  
        totCts <- sDen %>% group_by(Sample) %>% summarise(TotalSampleCounts=sum(TotalCounts))

        sDen <- sDen %>% 
               left_join(sia, by=c("Sample")) %>% 
               mutate(Density=TotalCounts/TotalSampleArea) %>% 
               left_join(totCts, by=c("Sample")) %>% 
               mutate(Pct=TotalCounts/TotalSampleCounts)


        sDen$CellTypeLabels <- unlist(sDen$CellType)
        if(!is.null(cellTypeLabels)){
            for(x in names(cellTypeLabels)){
                sDen$CellTypeLabels[which(sDen$CellType == x)] <- cellTypeLabels[[x]]
            }
        }
        sDen$CellTypeLabels <- factor(sDen$CellTypeLabels, levels=unique(cellTypeLabels))
        sDen[is.na(sDen)] <- 0

        keyLabels <- paste0(sDen$CellTypeLabels, " (",round(sDen$Pct*100),"%)")
        names(keyLabels)=sDen$CellTypeLabels

        pdf(file.path(sDir, paste0(s,"_summary.pdf")), height=8.5, width=11)
        p1 <- ggplot(sDen, aes(x="", y=Pct, fill=CellTypeLabels)) + 
                  geom_bar(stat="identity") + 
                  scale_fill_manual(values=clrs, labels=keyLabels) + 
                  coord_polar("y") + 
                  #getPlotTheme("infiltration density") + 
                  plotTheme +
                  theme(axis.text.x = element_blank(), 
                        axis.line = element_blank(), 
                        axis.title.x = element_blank(), 
                        axis.title.y = element_blank()) + 
                  labs(title=paste0(samplePlotTitle,"\n",s))
        print(p1)
        dev.off()
        allPlots[[s]][["summary"]] <- p1

        for(pop in unique(config$Population)){
            allPlots[[s]][[pop]] <- NULL
            print(paste0("  ",pop))

            popCfg <- config %>% filter(Population==pop)

            ctLabels <- popCfg$Label
            names(ctLabels) <- popCfg$CellType
            cellTypeLabels <- ctLabels

            popCfg$Color <- "gray"
            negIdx <- which(grepl("-",popCfg$Label))
            popClr <- unique(oviewCfg$Color[which(oviewCfg$CellType == pop)])
            clrs <- c(popClr, "darkgray")
            names(clrs) <- c("-","+")
            popCfg$Color[negIdx] <- popClr

            popCfg <- popCfg %>% filter(CellType != pop)            

            funcDen <- infiltrationDensity %>% filter(Sample==s, CellType %in% unique(popCfg$CellType))
            funcDen$Functional <- gsub("-","",gsub(paste0(pop,","),"",funcDen$CellType))
            funcDenSum <- funcDen %>% group_by(Sample,Functional,CellType) %>% summarise(Cts=sum(Counts))
 
            funcSum <- funcDen %>% 
                       group_by(Sample,Functional) %>% 
                       summarise(TotalCounts=sum(Counts))

            funcDen <- funcDenSum %>% 
                       left_join(funcSum, by=c("Sample","Functional")) %>%
                       left_join(sia, by=c("Sample")) %>%
                       mutate(Density=Cts/TotalSampleArea, Pct=Cts/TotalCounts)

            funcDen[is.na(funcDen)] <- 0 

            funcDen$CellTypeLabels <- unlist(funcDen$CellType)
            if(!is.null(cellTypeLabels)){
                for(x in names(cellTypeLabels)){
                    funcDen$CellTypeLabels[which(funcDen$CellType == x)] <- cellTypeLabels[[x]]
                }
            }

            funcDen$CellTypeLabels[which(grepl("-",funcDen$CellTypeLabels))] <- "-"
            funcDen$CellTypeLabels[-which(grepl("-",funcDen$CellTypeLabels))] <- "+"
            funcDen$CellTypeLabels <- factor(funcDen$CellTypeLabels, levels=names(clrs))

            funcDen$Functional <- gsub("-","",gsub(paste0(pop,","),"",funcDen$CellType))
            funcDen$Label <- ""
            pos <- which(funcDen$CellTypeLabels == "+")
            funcDen$Label[pos] <- paste0(round(funcDen$Pct[pos]*100),"%")
            

            titles <- funcDen %>%
                      filter(Label != "") %>% 
                      select(Sample,Functional, Label) %>% 
                      mutate(Title=paste0(Functional," (",Label,")")) %>% 
                      select(-(Label))

            funcDen <- funcDen %>% 
                       left_join(titles, by=c("Sample","Functional"))
            funcDen[is.na(funcDen)] <- 0

            pdf(file.path(sDir,paste0(s,"_",pop,".pdf")),height=8.5,width=11)
            p1 <- ggplot(funcDen, aes(x="", y=Pct, fill=CellTypeLabels)) + 
                  geom_bar(stat="identity") + 
                  facet_wrap(~Title) + 
                  scale_fill_manual(values=clrs, labels=names(clrs)) + 
                  coord_polar("y") + 
                  plotTheme + 
                  #getPlotTheme("infiltration density") + 
                  theme(axis.text.x = element_blank(), 
                        axis.line = element_blank(), 
                        axis.title.x = element_blank(), 
                        axis.title.y = element_blank()) +
                  labs(title=s, subtitle=pop)
            print(p1)
            dev.off()
            allPlots[[s]][[pop]] <- p1
        }
    }
    return(allPlots)
} # end function
