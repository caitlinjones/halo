## Analyze Halo data

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
OR
```{r eval=FALSE}
counts.R -m counts_manifest.txt
```
then
```{r eval=FALSE}
pie_charts.R -m pie_charts_manifest.txt
```
