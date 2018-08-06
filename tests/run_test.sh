##
## (1) Configure study; minimum required info:
##
Rscript ../scripts/configure_halo_pipeline.R \
     --rawDataDir /home/socci/Work/Users/MellingI/Halo/Melanoma_IL2__Final/Cohort2/data/Melanoma_IL2__Final_C2_v1 \
     --studyName halodevTest \
     --dataDir $PWD/objectAnalysisData \
     --metaDir /ifs/tcga/socci/Multiomyx/HaloData/Melanoma_IL2__Final/Cohort2/MetaData \
     --driftDir /ifs/tcga/socci/Multiomyx/Cell_drift_loss_masks/melanoma_drift_result/drift_summary \
     --markerConfigFile $PWD/halodevTest/config/all_marker_sets_config.yaml \
     --cellTypeConfigFile $PWD/halodevTest/config/cell_type_config.yaml \
     --plotConfigFile $PWD/halodevTest/config/plot_config.yaml \
     --annotationsDirs /ifs/tcga/socci/Multiomyx/HaloData/Melanoma_IL2__Final/Cohort2/HaloCoordinates \
     --setDefaultDirectoryStructure

##
## (2) Run pipeline
##
cd halodevTest
Rscript ../../scripts/final_pipeline.R -m config/study_config.yaml
