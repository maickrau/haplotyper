#!/bin/sh
#note: run generate_all_testdatas.sh before running this

runTest(){
	./remove_zero_columns.exe $1/snpsupport_simulated.txt $1/snpsupport_simulated_startunzero.txt $1/renumbering_simulated_startunzero.txt
	./collapse_multiples.exe $1/snpsupport_simulated_startunzero.txt $1/snpsupport_simulated_collapsed.txt $1/renumbering_simulated_collapse.txt
	./merge_subsets.exe $1/snpsupport_simulated_collapsed.txt $1/snpsupport_simulated_merged.txt $1/renumbering_simulated_merging.txt
	./merge_similars.exe $1/snpsupport_simulated_merged.txt 100 14 $1/snpsupport_simulated_remerged.txt $1/renumbering_simulated_remerging.txt
	./remove_zero_columns.exe $1/snpsupport_simulated_remerged.txt $1/snpsupport_simulated_unzero.txt $1/renumbering_simulated_unzero.txt
	./haplotyper_main.exe $1/snpsupport_simulated_unzero.txt 4 > $1/result_simulated_raw.txt
	./renumberer.exe $1/result_simulated_fixed.txt $1/result_simulated_raw.txt $1/renumbering_simulated_startunzero.txt $1/renumbering_simulated_collapse.txt $1/renumbering_simulated_merging.txt $1/renumbering_simulated_remerging.txt $1/renumbering_simulated_unzero.txt
	./verify_result.exe $1/result_simulated_fixed.txt $1/reads_all_repositioned.fasta > $1/result_simulated_score.txt
	./switch_distance_scorer.exe $1/result_simulated_fixed.txt $1/reads_all_repositioned.fasta B_C_D_E.fasta 10 4 > $1/result_simulated_switch.txt
	./evaluate_merge_success.exe $1/reads_all_repositioned.fasta $1/renumbering_simulated_remerging.txt $1/renumbering_simulated_startunzero.txt $1/renumbering_simulated_collapse.txt $1/renumbering_simulated_merging.txt > $1/merge_success_simulated.txt

	./prettify_sort.exe $1/snpsupport_simulated_unzero.txt $1/snpsupport_simulated_print.txt $1/renumbering_simulated_prettysort.txt
	./visualize_snpsupport.exe $1/snpsupport_simulated_print.txt $1/visualize_snpsupport_simulated_final.txt
}

runTest 100
runTest 200
runTest 300
runTest 400
runTest 500
runTest 600
runTest 700
runTest 800
runTest 900
runTest 1000
runTest 1200
runTest 1400
runTest 1600
runTest 1800
runTest 2000
