#!/bin/bash

# Exit on error
set -e

path="$1"
DB="$2"
 echo $path
# Functions
printToLog() {
	date +"%F__%T ---> $1"
}

usage () {
  printf "USAGE:\n"
  printf "metagen.sh <path_to_fastq> <DB>\n\n"
}

exitWithError () {
  printf "\n***ERROR: $1\n\n"
  usage
  exit
}

# Check params
if [ "$#" -ne 2 ]; then
  exitWithError "Missing parameters"
fi

# Create log directory
rm -rf "$path/logs"; mkdir "$path/logs"
rm -rf "$path/out"; mkdir "$path/out"
mkdir "$path/out/plots"

# Base Pipeline - generate required outputs
printf "\nRunning DADA2 pipeline...\n\n"
Rscript /home/makis/tools/metagen/dada2_pipeline.R "$path" "$DB" >"$path"/logs/out 2>"$path"/logs/error

# Run picrust2 pipeline
printf "\nRunning PICRUSt2 pipeline...\n\n"
printToLog "Running functional annotation..." >>"$path"/logs/out 2>>"$path"/logs/error
source ~/tools/miniconda3/etc/profile.d/conda.sh
conda activate picrust2
picrust2_pipeline.py -s "$path"/out/asv.fasta -i "$path"/out/my.biom -o "$path"/out/functionalAnn_out -p 16 >>"$path"/logs/out 2>>"$path"/logs/error

printToLog "Adding KEGG_ID descriptions..." >>"$path"/logs/out 2>>"$path"/logs/error
add_descriptions.py -i "$path"/out/functionalAnn_out/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC -o "$path"/out/functionalAnn_out/EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz >>"$path"/logs/out 2>>"$path"/logs/error
add_descriptions.py -i "$path"/out/functionalAnn_out/pathways_out/path_abun_unstrat.tsv.gz -m METACYC -o "$path"/out/functionalAnn_out/pathways_out/path_abun_unstrat_descrip.tsv.gz >>"$path"/logs/out 2>>"$path"/logs/error

conda deactivate

# Generate rmd report
printf "\nRunning phyloseq RMD pipeline...\n\n"
printToLog "Generating visualizations and Rmarkdown report ..." >>"$path"/logs/out 2>>"$path"/logs/error
Rscript -e "path = '$path'; DB= '$DB'; rmarkdown::render('/home/makis/tools/metagen/metagen.Rmd', output_file = paste0(path,'metagen_report.html'), output_dir = path)" >>$1/logs/out 2>>$1/logs/error

printf "\nAnalysis Complete! \n\n"
printf "\nResults in $path \n\n"