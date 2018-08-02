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
NOTE: Manifests will soon be built automatically using meta data files, but for now,
use template manifests in example/ folder for each step below and modify manually
as needed

### Step 3 (if applicable): Run ONLY ONE TIME - Mark exclusions
Given raw object analysis \*.rda files and exclusion files (including Halo XML, drift summaries
and FOV annotations), add an EXCLUDE column to data indicating which cells should not be analyzed. 
If data has already been marked, skip to Step 4.
```{r eval=FALSE}
Rscript scripts/mark_exclusions.R -m example/counts/exclusion_manifest.txt
```

### Step 3: Run one or more pre-written scripts 

#### * Generate counts of given markers
WARNING: This script has not been tested recently
```{r eval=FALSE}
Rscript scripts/counts.R -m example/counts/counts_manifest.txt
```

#### * Spatial Plots
Plot both total density and density by band
```{r eval=FALSE}
Rscript scripts/spatial_plots_2.R -m example/spatial_plots/density_manifest.txt --plotDensity --plotDensityByBand
```

NOTES: 
* Plotting total density has not been tested since modifying and adding code for plotting 
density by band
* Other scripts in scripts/ folder also have not been tested recently
* Manifests are a bit painful to create right now - they will soon be built automatically
  using meta data files and can be modified manually as needed
