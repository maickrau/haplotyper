#!/bin/sh

runTest(){
	./mutate_snpsupports.exe $1 tmpsupports.tmp $3
	./remove_zero_columns.exe tmpsupports.tmp tmpsupports2.tmp tmprenumbering1.tmp
	./collapse_multiples.exe tmpsupports2.tmp tmpsupports3.tmp tmprenumbering2.tmp
	./merge_subsets.exe tmpsupports3.tmp tmpsupports4.tmp tmprenumbering3.tmp
	./merge_similars.exe tmpsupports4.tmp 100 14 tmpsupports5.tmp tmprenumbering4.tmp
	./remove_zero_columns.exe tmpsupports5.tmp tmpsupports6.tmp tmprenumbering5.tmp
	./evaluate_merge_success.exe reads_all_repositioned.fasta tmprenumbering4.tmp tmprenumbering1.tmp tmprenumbering2.tmp tmprenumbering3.tmp > $7
	./haplotyper_main.exe tmpsupports6.tmp 4 > $2
	./renumberer.exe $4 $2 tmprenumbering1.tmp tmprenumbering2.tmp tmprenumbering3.tmp tmprenumbering4.tmp tmprenumbering5.tmp
	./verify_result.exe $4 reads_all_repositioned.fasta > $5
	./switch_distance_scorer.exe $4 reads_all_repositioned.fasta B_C_D_E.fasta 10 4 > $6
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
