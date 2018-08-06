### study_config.yaml
This configuration file contains all parameters for a pipeline run. 

### cell_type_config.yaml
This configuration file contains all cell type information needed to generate counts file of all possible marker combinations that represent cell types. It can be manually created or auto-generated based on *CellTypes.xlsx file in meta data directory. 

### marker_config.yaml
A marker configuration file contains all information needed to run area/density analyses and to plot results. Info includes marker sets (e.g., exhaustion, identity, etc.), marker combinations to plot for each set, and labels to use on plots.

### plot_config.yaml
A plot configuration file is used in conjunction with a marker config file to indicate which colors to use for certain marker combinations represented by labels in marker config.
