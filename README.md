**Overview**

This folder contains the major scripts used in my thesis titled : 
Exploring the transpositional landscape and recent transposable element activity in beech trees using long read mobilome and genome sequencing  and with new computational tools

**Structure**



**eccDNA_SR/**\
>>├ eccdna_sr_analysis.sh     # eccDNA peaks calling (using total read coverage and split read and disordant read coverage) and intersection with TE annotation using eccDNA Illumina samples\
>>├ eccdna_sr_figures_mobilome.Rmd    # Plot eccDNA heatmap\


**eccDNA_LR/**\
    ├ eccdna_lr_analysis.sh     # eccDNA peaks calling (using total read coverage and repetitive read coverage) and intersection with TE annotation using eccDNA ONT samples (Fig 1C.)\
      ├ repeats_mobilome0.py    # Find reads containing repetitive units (Fig 1C.)\
    ├ eccdna_lr_figures_mobilome.Rmd    # Plot eccDNA heatmap, PCA and barplot of abundant TEs in samples (Fig 1BDE.)\
    ├ squirrel1_2_figures_mobilome.Rmd    # Barlot of squirrel1 and squirrel2 abundance in samples, histogram of copy size and distance to closest gene distribution (Fig 3AFG.)\


**Somatic_SV_detection/**\
    ├ filter_reads.sh     # Filtering ONT WGS reads of #354M and #354R (Fig 2A)\
    ├ somatic_sv_detection.sh   # Structurant variant calling using ONT WGS reads of #354M and #354R (Fig 2A)\


**TIPs_detection/**\
    ├ eccdna_lr_analysis.sh     # eccDNA peaks calling (using total read coverage and repetitive read coverage) and intersection with TE annotation using eccDNA ONT samples (Fig 1C.)\
      ├ repeats_mobilome0.py    # Find reads containing repetitive units (Fig 1C.)\
    ├ eccdna_lr_figures_mobilome.Rmd    # Plot eccDNA heatmap, PCA and barplot of abundant TEs in samples (Fig 1BDE.)\
