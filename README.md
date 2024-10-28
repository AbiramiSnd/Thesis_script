## Overview
This folder contains the major scripts used in my thesis titled : 
Exploring the transpositional landscape and recent transposable element activity in beech trees using long read mobilome and genome sequencing  and with new computational tools
The thesis is funded by the ED305 doctoral school of Perpignan University, in France.

## Usage
The shell script are launch as : ./script.sh
The R script are included in the shell scripts.
The Rmd files to generate figures are launched under Rstudio. 


## Structure
Here is the detailed organisation of the folders and scripts. 

### TE_library_clustering
> ├ **clustering_te_library.sh**     # eccDNA peaks calling (using total read coverage and split read and disordant read coverage) and intersection with TE annotation using eccDNA Illumina samples\


### eccDNA_SR
> ├ **eccdna_sr_analysis.sh**     # eccDNA peaks calling (using total read coverage and split read and disordant read coverage) and intersection with TE annotation using eccDNA Illumina samples (Supp Fig 1C)\
> ├ **eccdna_sr_figures_mobilome.Rmd**    # Plot eccDNA heatmap (Supp Fig 2)


### eccDNA_LR
> ├ **eccdna_lr_analysis.sh**     # eccDNA peaks calling (using total read coverage and repetitive read coverage) and intersection with TE annotation using eccDNA ONT samples (Fig 1C.)
> > ├ **repeats_mobilome0.py**    # Find reads containing repetitive units (Fig 1C.)\
> ├ **eccdna_lr_figures_mobilome.Rmd**    # Plot eccDNA heatmap, PCA and barplot of abundant TEs in samples (Fig 1BDE; Supp Fig 3)\
> ├ **squirrel1_2_figures_mobilome.Rmd**    # Barlot of squirrel1 and squirrel2 abundance in samples, histogram of copy size and distance to closest gene distribution (Fig 3AFG.)


### Somatic_SV_detection
> ├ **filter_reads.sh**     # Filtering ONT WGS reads of #354M and #354R (Fig 2A)\
> ├ **somatic_sv_detection.sh**   # Structurant variant calling using ONT WGS reads of #354M and #354R (Fig 2A)


### TIPs_detection
> ├ **launch_trackposon.sh**     # Launch trackposon on all individuals (Fig 4C)
> > ├ **TRACKPOSON.sh**    # TRACKPOOSN tool to detect insertion of a TE consensus family per individual (Fig 4C)\
> ├ **launch_local_coverage.sh**     # Launch local coverage calculation on all individuals\
> > ├ **local_window_coverage.sh**    # Using mosdepth to calculate local coverage of all genomic windows containing a TIP to estimate allelic frequency\
> ├ **launch_pipeline_post.sh**     # Generate numerical matrices per TE family by concatenating bedfiles output of trackposon with columnwise individual and rowwise windows, using a loop for all TEs (Fig 4C)\
> > ├ **all_cov_matrix.R**    # Generate a numerical matrix of local coverage with row-wise windows and column-wise individuals
> > > ├ **analysis_pipeline.sh**     # Generate read bedfiles of insertions per TE of all indvidual (Fig 4C)
> > > > ├ **analysis_pipeline.R**    # Generate numerical matrix of coverage, read total coverage matrix and estimate allelic frequency and missing data and gererate matrix with 4 factors : 0: no insertion, 1: heterozygous insertion, 2: homozygous insertion, NA: no mapping coverage (Fig 4C)\
> ├ **all_tips_analysis.sh**    # Concatenate all LTR and MITE TIPs to identify most TIPs producing windows\
> ├ **tips_figures_mobilome.Rmd** # Plot TIPs figures (Fig 4DE; Fig 5ABCD; Supp Fig 10B; Supp Fig 11; Supp Fig 12)


### TIPs_GWAS
> ├ **launch_gwas_rmvp.sh**    # Using a loop to run GWAS on all TE matrices (Fig XI)
> > ├ **run_gwas_rmvp.sh**    # Run population structure calculation, and gwas script
> > > ├ **rmvp_gwas.sh**    # Run GWAS using numerical TE matrix, phenological dataset and popualtion striuture\
> > > ├ **rmvp_pop-struct.R**    # Calculate population structure using SNP calling

## Notes
The script have been updated on 28th october and will not be modified untill 6th december. 
These programs are still under developpemnt to build a robust tool in the future.
