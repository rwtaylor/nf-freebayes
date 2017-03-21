#!/bin/bash
#SBATCH --job-name=fb-10x
#SBATCH --time=2-0:0:0
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --partition=dpetrov,hns
#################

./nextflow run main.nf -resume -c nextflow.config --samples samples.tsv -with-timeline -with-trace -with-dag dag.png
