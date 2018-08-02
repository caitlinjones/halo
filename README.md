## Analyze Halo data
This package generates various statistics and plots for Halo data. Currently in its 
beginning stages, it has been written specifically for melanoma data but as it evolves
it ultimately will support all data sets.

### Step 1: Install
On MSK server, with R 
```{r eval=FALSE}
/home/byrne/R/R-3.4.3/bin/R
```
Run
```{r eval=FALSE}
library(devtools)
install_github("caitlinjones/halo")
```
### Step 2: Create manifest with all project parameters including input files
Example [here](example/manifest.txt) 

Run
```{r eval=FALSE}
cd [WORKING_DIRECTORY]
Rscript scripts/configure_halo_pipeline.R \
     --rawDataDir /home/byrne/halo/dev/halodev/tests/data \
     --studyName halodevTest \
     --dataDir $PWD/objectAnalysisData \
     --metaDir /ifs/tcga/socci/Multiomyx/HaloData/Melanoma_IL2__Final/Cohort2/MetaData \
     --driftDir /ifs/tcga/socci/Multiomyx/Cell_drift_loss_masks/melanoma_drift_result/drift_summary \
     --markerConfigFile /home/byrne/halo/data/template_configs/template_marker_config.yaml \
     --cellTypeConfigFile /home/byrne/halo/data/template_configs/template_marker_config.yaml \
     --plotConfigFile /home/byrne/halo/data/template_configs/plot_config.yaml \
     --annotationsDirs /ifs/tcga/socci/Multiomyx/HaloData/Melanoma_IL2__Final/Cohort2/HaloCoordinates \
     --setDefaultDirectoryStructure
```

Above are the minimum required arguments to start a new pipeline run FROM SCRATCH. For a full list of options, run
```{r eval=FALSE}
Rscript scripts/configure_halo_pipeline.R -h
```

Use '--setDefaultDirectoryStructure' to do just that. This will create default folders and subfolders for all possible analyses. 

The result of this script is a YAML file with all parameters needed to run entire pipeline. Default location for this file is studyName/config/study_config.yaml. This file can then be manually edited as needed or used as a template for future pipeline runs.

### Step 3: Run pipeline
Move to study directory, e.g., 
```{r eval=FALSE}
cd {studyName}
```

Run from scratch, including marking exclusions
```{r eval=FALSE}
Rscript scripts/final_pipeline.R -m config/study_config.yaml --markExclusions
```

NOTES: 
* Currently the pipeline runs the following steps:
   - parse and store all halo boundaries
   - mark exclusions & generate debug plots
   - generate cell type marker combination counts spreadsheet with cell type interpretations
   - calculate total FOV area and marker densities
   - plot total FOV marker densities
   - calculate infiltration area and marker densities (by distance intervals from tumor interfaces)
   - plot infiltration marker densities

* Exclusions need to be marked only once. Once run, study_config.yaml will be updated to include data_dir, which will point to the directory containing \*.rda files that include exclusions. For any subsequent pipeline runs, do NOT include markExclusions unless meta data, drift data or Halo boundary data has changed. 

* study_config.yaml will be updated as the pipeline runs to point to any new files generated during the run. This will allow steps to be skipped in future runs if their dependencies are unchanged.
 

 
