validateConfig <- function(config){

    ## DEFAULTS
    pad_default <- 20
    drift_default <- 0.1
    write_csv_files <- TRUE
    max_g <- 5
    band_width <- 10
    maxDistanceFromInterface <- 360


    ## make sure all directories and files exist and are readable/writeable
    dir_settings <- names(config)[grep("_dir",names(config))]
    for(dir in dir_settings){
        if(is.null(config[[dir]])){ 
            warning(paste0("NULL parameter: ",dir))
            next 
        }
        assert_that(is.dir(config[[dir]]))
        assert_that(is.readable(config[[dir]]))
        if(!grepl("drift|raw",dir)){ ## these dirs should NOT be writeable
            assert_that(is.writeable(config[[dir]]))
        }
    }

    ## make sure all files exist, and if any file param is NULL, 
    ## make sure the equivalent _dir param is NOT NULL
    file_settings <- names(config)[grep("_file",names(config))]
    for(fs in file_settings){
        if(is.null(config[[fs]])){ 
            warning(paste0("NULL parameter: ",fs))
            assert_that(!is.null(config[[gsub("_files","_dir",fs)]]))
            next 
        }
        for(f in config[[fs]]){
            assert_that(file.exists(f))
            assert_that(is.readable(f))
        }
    }

    ## REQUIRED: 
    #study config file

    if(is.null(config$meta_dir) && is.null(config$meta_files)){
        stop("Please provide either meta_dir or meta_files in config")
    }
    if(config$markExclusions){ 
        if(is.null(config$raw_data_dir) && is.null(config$raw_data_files)){
            stop("markExclusions is set to TRUE, but no raw_data_[dir|files] are set in config.")
        }
        if(is.null(config$drift_dir) && is.null(config$drift_files)){
            warning("No drift_[dir|files] set in config. No drift exclusions will be considered.")
        }
        if(is.null(config$drift_threshold)){
            config$drift_threshold <- 0.1
            warning(paste0("No drift threshold set in study config. Setting to default: ",drift_default))
        }
        if(is.null(config$pad)){
            config$pad <- 20
            warning(paste0("No padding set in study config. Setting to default: ",pad_default))
        }
    } else {
        if(is.null(config$data_dir) && is.null(data_files)){
            stop("Please provide at least one data_dir or data_files in config")
        }
    }
    if(is.null(config$annotations_dirs)){
        if(is.null(config$annotations_dirs)){
            stop("Please provide at least one annotations_dirs in config")
        }
    }

    return(config)
}


validateExclusions <- function(dat){
    excTypes <- c("LabExclusion","HALOExclusion","HALOGlass","HALOEpidermis","DRIFT","PAD")

    tmp <- dat %>% 
           group_by(Sample) %>% 
           select(UUID,EXCLUDE) %>% 
           unique() %>% 
           group_by(EXCLUDE) %>% 
           summarise(Count=n()) %>% 
           mutate(FracTotal=Count/sum(Count))

    for(t in excTypes){
        tmp %>% 
        filter(grepl(t, EXCLUDE)) %>% 
        summarise(TotalCoun=sum(Count),TotFrac=sum(FracTotal))        
    }
}

#' Make sure list of unique FOVs contains only one element
#'
#' If the list is longer than one element, print error message and EXIT 
#'
#' @param uniqFOVs    vector containing list of unique FOVs in dataset
#' @export
validateSingleFOVFile <- function(uniqFOVs){
    ## ensure that file contains only ONE FOV
    if(length(uniqFOVs) > 1){
        flog.fatal("Multiple FOVs found in this file.")
        stop("MULTIPLE FOVs found in file. Please split file by FOV and rerun")
    }
}

