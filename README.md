## Analyze Halo data

### Step 1: Install
On MSK server, with R 
```{r eval=FALSE}
/home/byrne/R/R-3.4.3/bin/R
```
Run
```{r eval=FALSE}
library(devtools)
install_github("halo")
```
### Step 2: Create manifest with all project parameters
Example here [TO DO: link to example]

### Step 3: Run pipeline to generate counts and/or plots
Run 
```{r eval=FALSE}
halo_pipeline.R -m my_project_manifest.txt 
```
to generate counts and/or plots for Object Analysis Data from Halo.