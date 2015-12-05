#!/bin/sh
#HOW TO RUN THIS EXPERIMENT:
#- first run the read length experiment to generate the haplotyping matrices
#- create a directory for this experiment
#- compile all source files
#- create a subdirectory for EVERY read length used in this experiment
#- NOTE: do not use the same directory names as the simulated script!
#- do the following for each read length:
#- copy all compiled executables to the read length subdirectory
#- get the following preprocessing results from the read length experiment:
#-    snpsupport_discpsnp.txt, reads_all.fasta, B_C_D_E.fasta
#- copy the above files to the read length subdirectory
#- copy error_rate_experiment_discosnp.sh to the read length subdirectory
#- in the read length subdirectory, run ./error_rate_experiment_discosnp.sh snpsupport_discosnp.txt
#- after each read length experiment for both DiscoSNP and simulated method have been run:
#- in the experiment directory, run parse_results.sh
#- results will be in the generated results_*.csv files

runTest(){
	./mutate_snpsupports.exe $1 tmpsupports.tmp $3
	./remove_zero_columns.exe tmpsupports.tmp tmpsupports1.tmp tmprenumbering1.tmp
	./collapse_multiples.exe tmpsupports1.tmp tmpsupports2.tmp tmprenumbering2.tmp
	./merge_subsets.exe tmpsupports2.tmp tmpsupports3.tmp tmprenumbering3.tmp
	./merge_similars.exe tmpsupports3.tmp 11 100 tmpsupports4.tmp tmprenumbering4.tmp
	./collapse_multiples.exe tmpsupports4.tmp tmpsupports5.tmp tmprenumbering5.tmp
	./remove_zero_columns.exe tmpsupports5.tmp tmpsupports6.tmp tmprenumbering6.tmp
	./c1p_solver.exe tmpsupports6.tmp tmprenumbering7.tmp tmpsupports7.tmp
	./matrix_bander.exe tmpsupports7.tmp 50 tmpsupports8.tmp tmprenumbering8.tmp 10 0,98 200000 5 1

	./remove_outliers.exe tmpsupports8.tmp tmpsupports9.tmp 11
	./collapse_multiples.exe tmpsupports9.tmp tmpsupports10.tmp tmprenumbering9.tmp
	./merge_subsets.exe tmpsupports10.tmp tmpsupports11.tmp tmprenumbering10.tmp
	./merge_similars.exe tmpsupports11.tmp 11 100 tmpsupports12.tmp tmprenumbering11.tmp
	./collapse_multiples.exe tmpsupports12.tmp tmpsupports13.tmp tmprenumbering12.tmp
	./remove_zero_columns.exe tmpsupports13.tmp tmpsupports14.tmp tmprenumbering13.tmp
	./c1p_solver.exe tmpsupports14.tmp tmprenumbering14.tmp tmpsupports15.tmp
	./matrix_bander.exe tmpsupports15.tmp 50 tmpsupports16.tmp tmprenumbering15.tmp 10 0,98 200000 5 1

	./merge_similars.exe tmpsupports16.tmp 100 14 tmpsupports17.tmp tmprenumbering16.tmp
	./remove_zero_columns.exe tmpsupports17.tmp tmpsupports18.tmp tmprenumbering17.tmp
	./haplotyper_main.exe tmpsupports18.tmp 4 > $2
	./renumberer.exe $4 $2 tmprenumbering1.tmp tmprenumbering2.tmp tmprenumbering3.tmp tmprenumbering4.tmp tmprenumbering5.tmp tmprenumbering6.tmp tmprenumbering7.tmp tmprenumbering8.tmp tmprenumbering9.tmp tmprenumbering10.tmp tmprenumbering11.tmp tmprenumbering12.tmp tmprenumbering13.tmp tmprenumbering14.tmp tmprenumbering15.tmp tmprenumbering16.tmp tmprenumbering17.tmp

	./verify_result.exe $4 reads_all.fasta > $5
	./switch_distance_scorer.exe $4 reads_all.fasta B_C_D_E.fasta 10 4 > $6

	./evaluate_merge_success.exe reads_all.fasta tmprenumbering17.tmp tmprenumbering1.tmp tmprenumbering2.tmp tmprenumbering3.tmp tmprenumbering4.tmp tmprenumbering5.tmp tmprenumbering6.tmp tmprenumbering7.tmp tmprenumbering8.tmp tmprenumbering9.tmp tmprenumbering10.tmp tmprenumbering11.tmp tmprenumbering12.tmp tmprenumbering13.tmp tmprenumbering14.tmp tmprenumbering15.tmp tmprenumbering16.tmp > $7

}

runTest2(){
	runTest $1 result_raw_0_$2.txt 0.00 result_fixed_0_$2.txt result_score_0_$2.txt result_switch_0_$2.txt result_mergesuccess_0_$2.txt
	runTest $1 result_raw_1_$2.txt 0.01 result_fixed_1_$2.txt result_score_1_$2.txt result_switch_1_$2.txt result_mergesuccess_1_$2.txt
	runTest $1 result_raw_2_$2.txt 0.02 result_fixed_2_$2.txt result_score_2_$2.txt result_switch_2_$2.txt result_mergesuccess_2_$2.txt
	runTest $1 result_raw_3_$2.txt 0.03 result_fixed_3_$2.txt result_score_3_$2.txt result_switch_3_$2.txt result_mergesuccess_3_$2.txt
	runTest $1 result_raw_4_$2.txt 0.04 result_fixed_4_$2.txt result_score_4_$2.txt result_switch_4_$2.txt result_mergesuccess_4_$2.txt
	runTest $1 result_raw_5_$2.txt 0.05 result_fixed_5_$2.txt result_score_5_$2.txt result_switch_5_$2.txt result_mergesuccess_5_$2.txt
	runTest $1 result_raw_6_$2.txt 0.06 result_fixed_6_$2.txt result_score_6_$2.txt result_switch_6_$2.txt result_mergesuccess_6_$2.txt
	runTest $1 result_raw_7_$2.txt 0.07 result_fixed_7_$2.txt result_score_7_$2.txt result_switch_7_$2.txt result_mergesuccess_7_$2.txt
	runTest $1 result_raw_8_$2.txt 0.08 result_fixed_8_$2.txt result_score_8_$2.txt result_switch_8_$2.txt result_mergesuccess_8_$2.txt
	runTest $1 result_raw_9_$2.txt 0.09 result_fixed_9_$2.txt result_score_9_$2.txt result_switch_9_$2.txt result_mergesuccess_9_$2.txt
	runTest $1 result_raw_10_$2.txt 0.10 result_fixed_10_$2.txt result_score_10_$2.txt result_switch_10_$2.txt result_mergesuccess_10_$2.txt
	runTest $1 result_raw_11_$2.txt 0.11 result_fixed_11_$2.txt result_score_11_$2.txt result_switch_11_$2.txt result_mergesuccess_11_$2.txt
	runTest $1 result_raw_12_$2.txt 0.12 result_fixed_12_$2.txt result_score_12_$2.txt result_switch_12_$2.txt result_mergesuccess_12_$2.txt
	runTest $1 result_raw_13_$2.txt 0.13 result_fixed_13_$2.txt result_score_13_$2.txt result_switch_13_$2.txt result_mergesuccess_13_$2.txt
	runTest $1 result_raw_14_$2.txt 0.14 result_fixed_14_$2.txt result_score_14_$2.txt result_switch_14_$2.txt result_mergesuccess_14_$2.txt
	runTest $1 result_raw_15_$2.txt 0.15 result_fixed_15_$2.txt result_score_15_$2.txt result_switch_15_$2.txt result_mergesuccess_15_$2.txt
}

runTest2 $1 1
runTest2 $1 2
runTest2 $1 3
runTest2 $1 4
runTest2 $1 5
runTest2 $1 6
runTest2 $1 7
runTest2 $1 8
runTest2 $1 9
runTest2 $1 10
