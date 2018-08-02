##
## (1) Configure study; minimum required info:
##
Rscript ../scripts/configure_halo_pipeline.R \
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

##
## (2) Run pipeline
##
cd halodevTest
Rscript ../../scripts/final_pipeline.R -m config/study_config.yaml
