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
### Step 2: Create manifest with all project parameters
Example [here](example/manifest.txt) 

### Step 3: Run pipeline or one step at a time
To generate counts and/or plots of Object Analysis data from Halo, run 
```{r eval=FALSE}
halo_pipeline.R -m my_project_manifest.txt 
```
Note: The pipeline currently only runs counts and pie charts.


Alternatively, run each step separately. For example,
```{r eval=FALSE}
Rscript scripts/counts.R -m example/counts/counts_manifest.txt
```
then
```{r eval=FALSE}
Rscript scripts/pie_charts.R -m example/pie_charts/pie_charts_manifest.txt
```

### Spatial Plots
Spatial plots are generated per sample, so for a single sample, "CR", for example, run
```{r eval=FALSE}
Rscript scripts/spatial_plots.R -m example/spatial_plots/CR_manifest.txt
```
Example shell script to run all samples: example/spatial_plots/run_spatial_plots.R
