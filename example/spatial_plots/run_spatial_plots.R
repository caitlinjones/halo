ls *manifest.txt | xargs -I {} Rscript ../../scripts/spatial_plots.R -m {}
