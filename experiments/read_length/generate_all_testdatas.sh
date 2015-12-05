#!/bin/sh
#HOW TO RUN THIS EXPERIMENT:
#- create a directory for this experiment
#- compile everything in the source directory and copy them to the experiment directory
#- copy the scripts to the experiment directory
#- acquire 4 genomes, name the files B.fasta, C.fasta, D.fasta, E.fasta, and name the genomes B, C, D, E respectively
#- also put the genomes to a file B_C_D_E.fasta
#- get DiscoSNP and make sure that the path in generate_snpsupports.sh refers to it
#-   NOTE: the path in generate_snpsupports.sh is relative to the subdirectory for that read length, not the experiment directory!
#- run ./generate_all_testdatas.sh
#- run ./generate_simulated_results.sh
#- run ./parse_results.sh
#- the results will be in the generated results_*.csv files


./generate_snpsupports.sh 100 2000
./generate_snpsupports.sh 200 1000
./generate_snpsupports.sh 300 666
./generate_snpsupports.sh 400 500
./generate_snpsupports.sh 500 400
./generate_snpsupports.sh 600 333
./generate_snpsupports.sh 700 285
./generate_snpsupports.sh 800 250
./generate_snpsupports.sh 900 222
./generate_snpsupports.sh 1000 200
./generate_snpsupports.sh 1200 166
./generate_snpsupports.sh 1400 142
./generate_snpsupports.sh 1600 125
./generate_snpsupports.sh 1800 111
./generate_snpsupports.sh 2000 100

./generate_results.sh 100
./generate_results.sh 200
./generate_results.sh 300
./generate_results.sh 400
./generate_results.sh 500
./generate_results.sh 600
./generate_results.sh 700
./generate_results.sh 800
./generate_results.sh 900
./generate_results.sh 1000
./generate_results.sh 1200
./generate_results.sh 1400
./generate_results.sh 1600
./generate_results.sh 1800
./generate_results.sh 2000
