#!/bin/sh

./haplotyper_main.exe $1/snpsupport_final.txt 4 > $1/result_raw.txt
./renumberer.exe $1/result_fixed.txt $1/result_raw.txt $1/renumbering_startunzero.txt $1/renumbering_collapsing.txt $1/renumbering_merging.txt $1/renumbering_remerging.txt $1/renumbering_recollapsing.txt $1/renumbering_unzero.txt $1/renumbering_c1p.txt $1/renumbering_banding.txt $1/renumbering_collapsing2.txt $1/renumbering_merging2.txt $1/renumbering_remerging2.txt $1/renumbering_recollapsing2.txt $1/renumbering_unzero2.txt $1/renumbering_c1p2.txt $1/renumbering_banding2.txt $1/renumbering_finalmerge.txt $1/renumbering_finalunzero.txt

./verify_result.exe $1/result_fixed.txt $1/reads_all.fasta > $1/result_score.txt
./switch_distance_scorer.exe $1/result_fixed.txt $1/reads_all.fasta B_C_D_E.fasta 10 4 > $1/result_switch.txt

#./evaluate_merge_success.exe $1/reads_all.fasta $1/renumbering_collapsing.txt $1/renumbering_startunzero.txt > $1/merge_success_collapsing1.txt
#./evaluate_merge_success.exe $1/reads_all.fasta $1/renumbering_merging.txt $1/renumbering_startunzero.txt $1/renumbering_collapsing.txt > $1/merge_success_subsets1.txt
#./evaluate_merge_success.exe $1/reads_all.fasta $1/renumbering_remerging.txt $1/renumbering_startunzero.txt $1/renumbering_collapsing.txt $1/renumbering_merging.txt > $1/merge_success_similars1.txt
#./evaluate_merge_success.exe $1/reads_all.fasta $1/renumbering_collapsing2.txt $1/renumbering_startunzero.txt $1/renumbering_collapsing.txt $1/renumbering_merging.txt $1/renumbering_remerging.txt $1/renumbering_recollapsing.txt $1/renumbering_unzero.txt $1/renumbering_c1p.txt $1/renumbering_banding.txt > $1/merge_success_collapsing2.txt
#./evaluate_merge_success.exe $1/reads_all.fasta $1/renumbering_merging2.txt $1/renumbering_startunzero.txt $1/renumbering_collapsing.txt $1/renumbering_merging.txt $1/renumbering_remerging.txt $1/renumbering_recollapsing.txt $1/renumbering_unzero.txt $1/renumbering_c1p.txt $1/renumbering_banding.txt $1/renumbering_collapsing2.txt > $1/merge_success_subsets2.txt
./evaluate_merge_success.exe $1/reads_all.fasta $1/renumbering_remerging2.txt $1/renumbering_startunzero.txt $1/renumbering_collapsing.txt $1/renumbering_merging.txt $1/renumbering_remerging.txt $1/renumbering_recollapsing.txt $1/renumbering_unzero.txt $1/renumbering_c1p.txt $1/renumbering_banding.txt $1/renumbering_collapsing2.txt $1/renumbering_merging2.txt > $1/merge_success_similars2.txt
