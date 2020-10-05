# Metagen

Taxonomic analysis for 16s rRna single-end files

## Usage

```sh
$ ./metagen.sh <path_to_fastQ> <ref_tax_DB>
```

## Required Tools

The following tools are required

* FastQ Screen - Map filtered fastQ files with contaminant non-bacterial reference genomes hg19, mm10 (https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/)
* R version 3.6.2 libraries
    * DADA2 - Standard pipeline for taxonomy assignment (https://benjjneb.github.io/dada2/tutorial.html)
    * Phyloseq - Taxonomy visualisations (https://joey711.github.io/phyloseq/index.html)
* R Data management and formatting libraries
    * Biostrings
    * ggplot2
    * ape
    * biomformat
    * readr
    * scales
    * rmarkdown
    * dplyr
    * tibble
    * treemap
    * highcharter
    * kableExtra
    * tidyr
    * plotly
    * heatmaply
* PICRUSt2 - generate functional predictions from amplicon sequences for treemap and pathway analysis (https://github.com/picrust/picrust2/wiki)

## Important notes
1. The pipeline is currently designed and optimised for IonTorrent single-end data
2. Taxonomy libraries that are installed and can currently be used are:
    - Silva-v.132
    - GreenGenes-v.13_8
    - RefSeq-RDP16S_v2_May2018
3. Input fastQ file must have the naming format "subsample_sampleName.fastq". The string before the underscore will be used to declare the analysis' sample categories.
4. Output is stored in the metagen_report.html file that will be generated nside the analysis' directory


## Example
```sh
$ ./oncompnet.sh /home/user/16_analysis_directory greengenes
```
