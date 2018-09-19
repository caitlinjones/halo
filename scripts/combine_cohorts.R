
options(stringsAsFactors = FALSE)
options( java.parameters = c("-Xss2560k", "-Xmx8g") )
options('python_cmd'='/opt/common/CentOS_6-dev/python/python-3.5.1/bin/python')

write("Loading libraries...",stdout())
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("halodev"))

#### for COMBINED cohorts
args <- list(manifest="config/study_config.yaml",
             fovStats=TRUE,
             noFOVplots=FALSE,
             forceAllPlotting=FALSE,
             infiltrationStats=FALSE)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

################################################
####           INITIALIZE PROJECT           ####
################################################

## load all existing project data to be analyzed
proj    <- cellDive.initCombinedProject(args)
pp      <- proj$pp
updated <- proj$updated
## pull values in dat list out into their own variables
for(var in names(proj$dat)){
    assign(var, proj$dat[[var]])
}


################################################
#########        START  ANALYSES       #########
################################################

annotationCols <- c("IL2 treated", "Lesion Response")
#annotationCols <- NULL

###
### plot total FOV densities
###
if(!args$noFOVplots && (updated$FOVdensity || args$forceAllPlotting)){
    flog.info("Printing total FOV density plots...")
    for(ms in unique(parsedMarkerConfig$MarkerSet)){
        msConfig <- parsedMarkerConfig %>% filter(MarkerSet == ms)
        den <- fovDensity %>% filter(CellType %in% unique(msConfig$CellType), !is.na(Counts))        
        den[["IL2 treated"]] <- den$Treatment
        den[["IL2 treated"]][!grepl("IL2",den$Treatment)] <- "-"
        den[["IL2 treated"]][grepl("IL2",den$Treatment)] <- "+"
        printTotalFOVDensityPlots(den, fovAreas, msConfig, yScaleConsistency="population", 
                           absoluteDensity=TRUE, densityPercentage=TRUE, summarize=TRUE, 
                           stacked=TRUE, sampleOrder=sampleOrder, outDir=pp$fov_density_dir, 
                           forceStack=FALSE, annotationCols=annotationCols)
    }

    ### plot cell identity densities
    if("cell_type_config_file" %in% names(pp)){
        parsedCellTypeConfig <- getMarkerConfig(pp$cell_type_config_file, pp$plot_config_file) %>%
                                filter(CellType != "Tumor")
        cellIdDen <- summarizeFOVDataByCellTypeDefinition(fovDensity, parsedCellTypeConfig, 
                                                          fovAreas, pp$plot_config_file)
        den <- cellIdDen$den %>% filter(CellType %in% cellIdDen$config$CellType, !is.na(Counts))
        den$Sample <- factor(den$Sample, levels=sampleOrder)
        den[["IL2 treated"]] <- den$Treatment
        den[["IL2 treated"]][!grepl("IL2",den$Treatment)] <- "-"
        den[["IL2 treated"]][grepl("IL2",den$Treatment)] <- "+"        

        printTotalFOVDensityPlots(den, fovAreas, cellIdDen$config, yScaleConsistency="population",
                               absoluteDensity=TRUE, densityPercentage=TRUE, summarize=TRUE,
                               stacked=TRUE, sampleOrder=sampleOrder, outDir=pp$fov_density_dir, 
                               forceStack=TRUE, annotationCols=annotationCols) 
    }

    ### cell type density heatmap

    ##### TO DO: PULL ALL OF THIS INFO OUT INTO A CONFIG FILE
    ## Prior Rx (systemic)
    prx <- c("#0c2c84","#df65b0","gray")
    names(prx) <- c("ICI", "non-ICI", "None")
    ## Num. IL2 injections
    ninj <- c("#9ed891", "#64875d", "#344930", "gray")
    names(ninj) <- c("<4", "4-10", ">10", "None") 
    ## Injection interval
    inj <- c("#084594","#4292c6","#9ecae1","gray")
    names(inj) <- c("<15", "15-60", ">60", "None")
    ## Lesion Response
    resp <- c("gray","orange","#32cd32")
    names(resp) <- sort(unique(allStudyAnn$`Lesion Response`), decreasing=TRUE) 
    ## Treatment
    tmt <- c("#0570b0","#d73027","#01665e","gray")
    names(tmt) <- c("IL2","IL2+","TVEC+","UT")
    ## Patient
    pt <- c("red","orange","yellow","#32cd32","purple","blue","lightblue")
    names(pt) <- sort(unique(allStudyAnn$Patient))

    annColors <- list(`Prior Rx (systemic)` = prx,
                      `Num. IL2 injections` = ninj,
                      `IL2 interval (days)` = inj,
                      `Lesion Response` = resp,
                      `Treatment` = tmt,
                      `Patient` = pt)                          
    htmpCfg <- completeMarkerConfig
    htmpCfg$CellType <- gsub(",PCK26-","",htmpCfg$CellType)
    htmpCfg$CellType <- gsub("CD45-,","",htmpCfg$CellType)
    fovDensity$CellType <- gsub("CD45-,","",fovDensity$CellType)

    htmpCfg <- htmpCfg %>% unique()
    htmpCfg <- htmpCfg %>% filter(!MarkerSet %in% c("exhaustion", 
                                                    "tim3_macrophages",
                                                    "tim3_nk_and_unknown_cells",
                                                    "T_cell"))
    htmpCfg$MarkerSet[which(htmpCfg$MarkerSet == "Macrophage_or_myeloid")] <- "Macrophages" 
    den <- fovDensity %>% filter(CellType %in% htmpCfg$CellType)

    pdf(file.path(pp$fov_density_dir, "cell_type_density_heatmaps.pdf"), height=14, width=24)
    makeDensityHeatmap(fovDensity, annot=allStudyAnn, msToSummarize=c("Macrophages"), annotColors=annColors, 
                       ctConfig=htmpCfg, areas=fovAreas, separateLegend=FALSE,byFOV=TRUE,bySample=TRUE)
    dev.off()

} else {
    flog.info("No changes to total FOV densities. No plotting necessary.")
}


if(args$infiltrationStats){

    ###
    ### plot infiltration densities
    ###
    if(!args$noInfiltrationPlots && (updated$InfiltrationDensity || args$forceAllPlotting)){

        ###
        ### plot infiltration densities
        ###

        byFOV <- TRUE
        indivPopulations <- TRUE
        ptsPerPage <- 7
        #sumCols <- 3
        #sumRows <- 3

        ### NOTE: it is not necessary to do this for each marker set (i.e., you can pass the entire
        ###       parsedMarkerConfig to the printInfiltrationDensityPlots() function, but this way, 
        ###       plots will be printed for each marker set before generating plots for the next set
        for(ms in unique(parsedMarkerConfig$MarkerSet)){
            msConfig <- parsedMarkerConfig %>% filter(MarkerSet == ms)
            den <- infiltrationDensity %>% filter(CellType %in% unique(msConfig$CellType))
            printInfiltrationDensityPlots(den, bandWidth=pp$band_width, msConfig, 
                   yScaleConsistency="population", absoluteDensity=TRUE, densityPercentage=TRUE, byFOV=byFOV, 
                   summarize=TRUE, stacked=TRUE, sampleOrder=sampleOrder, infiltrationAreas=infiltrationAreas, 
                   outDir=pp$infiltration_density_dir, indivPopulations=indivPopulations, sampleSummaryPtsPerPage=ptsPerPage,
                   sampleSummaryCols=sumCols,sampleSummaryRows=sumRows)
        }

        ### plot cell identity densities
        if("cell_type_config_file" %in% names(pp)){

            cellIdDen <- summarizeByCellTypeDefinition(bandAssignments, infiltrationAreas, parsedCellTypeConfig, pp$plot_config_file)
            den <- cellIdDen$den %>% filter(CellType != "Tumor")
            den$Sample <- factor(den$Sample, levels=sampleOrder)
            den <- den %>% left_join(fovAreas, by=colnames(den)[which(colnames(den) %in% colnames(fovAreas))])

            allStudyAnnot$Sample <- allStudyAnnot$SampleLabel
            den <- den %>% left_join(allStudyAnn, by=c("Sample","SPOT"))

            config <- cellIdDen$config %>% filter(CellType != "Tumor")
            printInfiltrationDensityPlots(den, bandWidth=pp$band_width, config, 
                       yScaleConsistency="population", absoluteDensity=TRUE, densityPercentage=TRUE, byFOV=byFOV, 
                       summarize=TRUE, stacked=TRUE, sampleOrder=sampleOrder, infiltrationAreas=infiltrationAreas, 
                       outDir=pp$infiltration_density_dir, forceStack=TRUE, indivPopulations=indivPopulations, 
                       sampleSummaryPtsPerPage=ptsPerPage,sampleSummaryCols=sumCols,sampleSummaryRows=sumRows)    
        }
    }
}
