## Overview
This folder contains the major scripts used in my thesis titled : 
Exploring the transpositional landscape and recent transposable element activity in beech trees using long read mobilome and genome sequencing  and with new computational tools


## Usage
The shell script are launch as : ./script.sh
The R script are included in the shell scripts.
The Rmd files to generate figures are launched under Rstudio. 


## Structure
Here is the detailed organisation of the folders and scripts. 


### eccDNA_SR
> ├ eccdna_sr_analysis.sh     # eccDNA peaks calling (using total read coverage and split read and disordant read coverage) and intersection with TE annotation using eccDNA Illumina samples\
> ├ eccdna_sr_figures_mobilome.Rmd    # Plot eccDNA heatmap


### eccDNA_LR
> ├ eccdna_lr_analysis.sh     # eccDNA peaks calling (using total read coverage and repetitive read coverage) and intersection with TE annotation using eccDNA ONT samples (Fig 1C.)
> > ├ repeats_mobilome0.py    # Find reads containing repetitive units (Fig 1C.)\
> ├ eccdna_lr_figures_mobilome.Rmd    # Plot eccDNA heatmap, PCA and barplot of abundant TEs in samples (Fig 1BDE.)\
> ├ squirrel1_2_figures_mobilome.Rmd    # Barlot of squirrel1 and squirrel2 abundance in samples, histogram of copy size and distance to closest gene distribution (Fig 3AFG.)


### Somatic_SV_detection
> ├ filter_reads.sh     # Filtering ONT WGS reads of #354M and #354R (Fig 2A)\
> ├ somatic_sv_detection.sh   # Structurant variant calling using ONT WGS reads of #354M and #354R (Fig 2A)


### TIPs_detection
> ├ launch_trackposon.sh     # Launch trackposon on all individuals
> > ├ TRACKPOSON.sh    # TRACKPOOSN tool to detect insertion of a TE consensus family per individual\
> ├ launch_local_coverage.sh     # Launch local coverage calculation on all individuals\
> > ├ local_window_coverage.sh    # Using mosdepth to calculate local coverage of all genomic windows containing a TIP to estimate allelic frequency\
> ├ launch_pipeline_post.sh     # Generate numerical matrices per TE family by concatenating bedfiles output of trackposon with columnwise individual and rowwise windows, using a loop for all TEs\
> > ├ all_cov_matrix.R    # Generate a numerical matrix of local coverage with row-wise windows and column-wise individuals\
> > > ├ Analysis_pipeline.sh     # Generate read bedfiles of insertions per TE of all indvidual\
> > > > ├ Analysis_pipeline.R    # Generate numerical matrix of coverage, read total coverage matrix and estimate allelic frequency and missing data and gererate matrix with 4 factors : 0: no insertion, 1: heterozygous insertion, 2: homozygous insertion, NA: no mapping coverage\
> ├ all_tips_analysis.sh    # Concatenate all LTR and MITE TIPs to identify most TIPs producing windows\
> ├ tips_figures_mobilome.Rmd # Plot TIPs figures


### TIPs_GWAS
> ├ launch_gwas_rmvp.sh    # Using a loop to run GWAS on all TE matrices
> > ├ run_gwas_rmvp.sh    # Run population structure calculation, and gwas script\
> > > ├ rmvp_gwas.sh    # Run GWAS using numerical TE matrix, phenological dataset and popualtion striuture\
> > > ├ rmvp_pop-struct.R    # Calculate population structure using SNP calling

**Notes**
The script have been updated on 28th october and will not be modified untill 6th december. 
These scripts are still underdeveloppemnt to build a robust tool in the future.
