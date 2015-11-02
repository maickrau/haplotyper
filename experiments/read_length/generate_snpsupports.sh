#!/bin/sh

mkdir $1
./make_reads.exe B.fasta $1/reads_B.fasta $2 $1
./make_reads.exe C.fasta $1/reads_C.fasta $2 $1
./make_reads.exe D.fasta $1/reads_D.fasta $2 $1
./make_reads.exe E.fasta $1/reads_E.fasta $2 $1
cat $1/reads_*.fasta > $1/reads_all.fasta
./reposition_reads.exe $1/reads_all.fasta $1/reads_all_repositioned.fasta
./simulate_SNPfinding.exe $1/reads_all_repositioned.fasta B_C_D_E.fasta $1/snpsupport_simulated.txt $1
./visualize_snpsupport.exe $1/snpsupport_simulated.txt $1/visualize_snpsupport_simulated.txt
cd $1
../../discosnp_test/DiscoSNP++-2.1.6-Source/run_discoSnp++.sh -P 10 -b 2 -r "reads_all.fasta"
cd ..
./brute_force_mapper.exe $1/reads_all.fasta $1/discoRes_k_31_c_4_D_0_P_10_b_2_withlow_coherent.fa $1/snpsupport_discosnp.txt
./remove_zero_columns.exe $1/snpsupport_discosnp.txt $1/snpsupport_startunzero.txt $1/renumbering_startunzero.txt
./collapse_multiples.exe $1/snpsupport_startunzero.txt $1/snpsupport_collapsed.txt $1/renumbering_collapsing.txt
./merge_subsets.exe $1/snpsupport_collapsed.txt $1/snpsupport_merged.txt $1/renumbering_merging.txt
./merge_similars.exe $1/snpsupport_merged.txt 11 100 $1/snpsupport_similared.txt $1/renumbering_remerging.txt
./collapse_multiples.exe $1/snpsupport_similared.txt $1/snpsupport_recollapsed.txt $1/renumbering_recollapsing.txt
./remove_zero_columns.exe $1/snpsupport_recollapsed.txt $1/snpsupport_unzero.txt $1/renumbering_unzero.txt
./c1p_solver.exe $1/snpsupport_unzero.txt $1/renumbering_c1p.txt $1/snpsupport_c1p.txt
./matrix_bander.exe $1/snpsupport_c1p.txt 50 $1/snpsupport_banded.txt $1/renumbering_banding.txt 10 0,98 200000 5 1

./remove_outliers.exe $1/snpsupport_banded.txt $1/snpsupport_outliers.txt 11
./collapse_multiples.exe $1/snpsupport_outliers.txt $1/snpsupport_collapsed2.txt $1/renumbering_collapsing2.txt
./merge_subsets.exe $1/snpsupport_collapsed2.txt $1/snpsupport_merged2.txt $1/renumbering_merging2.txt
./merge_similars.exe $1/snpsupport_merged2.txt 11 100 $1/snpsupport_similared2.txt $1/renumbering_remerging2.txt
./collapse_multiples.exe $1/snpsupport_similared2.txt $1/snpsupport_recollapsed2.txt $1/renumbering_recollapsing2.txt
./remove_zero_columns.exe $1/snpsupport_recollapsed2.txt $1/snpsupport_unzero2.txt $1/renumbering_unzero2.txt
./c1p_solver.exe $1/snpsupport_unzero2.txt $1/renumbering_c1p2.txt $1/snpsupport_c1p2.txt
./matrix_bander.exe $1/snpsupport_c1p2.txt 50 $1/snpsupport_banded2.txt $1/renumbering_banding2.txt 10 0,98 200000 5 1

./merge_similars.exe $1/snpsupport_banded2.txt 100 14 $1/snpsupport_finalmerge.txt $1/renumbering_finalmerge.txt
./remove_zero_columns.exe $1/snpsupport_finalmerge.txt $1/snpsupport_final.txt $1/renumbering_finalunzero.txt

./prettify_sort.exe $1/snpsupport_final.txt $1/snpsupport_print.txt $1/renumbering_prettysort.txt
./visualize_snpsupport.exe $1/snpsupport_print.txt $1/visualize_snpsupport_final.txt
./active_rows.exe $1/snpsupport_print.txt > $1/active_rows.txt
