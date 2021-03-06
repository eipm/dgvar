# DGVar

![workflow](eg/workflow.png)

## Overview
DGVar is developed for identifying high confidence germline deleterious variants (DGVs) based on whole-exome sequencing data. The tool takes BAM files or VCF files as input, and screens DGVs by going through a series of flitering steps (please see the above workflow chart). 

## System requirements
DGVar was developed and tested on linux system. To run DGVar, you'll need a linux system and the following softwares being installed:
- Python v2.7+
- Java v1.6 & 1.7
- Samtools v0.1.19
- Genome Analysis Toolkit v2.5.2
- snpEff v4.2

## Pre-compiled database
The following databases are required:
- Human GRCh37 reference genome (b37, available at GATK resource [bundle](https://software.broadinstitute.org/gatk/download/bundle))
- NCBI dbSNP v137 (available at GATK resource [bundle](https://software.broadinstitute.org/gatk/download/bundle))
- SnpEff database GRCh37.75 (available at [SnpEff](http://snpeff.sourceforge.net/download.html#databases))
- NCBI ClinVar v20180805 (slim version available at [Box](https://wcm.box.com/s/nzzhudb6371cwuv1w2omtp84ripsfnsr))
- ExAC v0.3.1 (slim version available at [Box](https://wcm.box.com/s/8t01e5eb3vf3f1idpg7zwq71ypos8qgz))
- Exome target regions (Agilent HaloPlex bed file available at [Box](https://wcm.box.com/s/4bkw0f2rn858re30hq85lwgxxk33mes3))
- Onco/TSG annotations (available at [here](db/gene_annotations.txt))
- Common variants in EIPM cohort (available at [here](db/eipm_common_variants.txt))

## Installation
Download the codes as well as the required database; update the path in the main shell script `call_variants.sh`, then run a quick test:
- Make sure the current working directory is where dgvar is installed, and then type `cd eg`
- Run a test by typing `sh test.sh`
- Check results at `results/test/annotated/filt/test.candidates.filter_common.txt.gz`, you should see the variant in BRCA2 gene (13:32914437 GT --> G).

## Instructions for use
To run DGVar:
- Update the path in the main shell script `call_variants.sh`
- Run shell script `sh  call_variants.sh  in.bam  sampleID  target.bed`

  where in.bam is the input bam file, sampleID is an unique id (no space allowed), target.bed is the bed file listing genomic regions to check for variants

Run an example using the script in the `eg` folder.

Expected run time may vary (10-30 min) depending on input file size and your computing resources.


