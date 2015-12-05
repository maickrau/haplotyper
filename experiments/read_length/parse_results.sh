#!/bin/sh

parse(){
	echo $1 | tr '\n' '\t' >> $2
	cat $1/result_score.txt | sed -n '1~5p' | tr '\.' '\,' >> $2

	echo $1 | tr '\n' '\t' >> $3
	cat $1/result_switch.txt | sed -n '11~14p' | tr ' ' '\t' >> $3

	echo $1 | tr '\n' '\t' >> $9
	cat $1/result_switch.txt | sed -n '13~14p' | tr ' ' '\t' >> $9

	echo $1 | tr '\n' '\t' >> $4
	cat $1/merge_success_similars2.txt | sed -n '7~7p' | tr '\.' ',' >> $4

	echo $1 | tr '\n' '\t' >> $5
	cat $1/result_simulated_score.txt | sed -n '1~5p' | tr '\.' '\,' >> $5

	echo $1 | tr '\n' '\t' >> $6
	cat $1/result_simulated_switch.txt | sed -n '11~14p' | tr ' ' '\t' >> $6

	echo $1 | tr '\n' '\t' >> $8
	cat $1/result_simulated_switch.txt | sed -n '13~14p' | tr ' ' '\t' >> $8

	echo $1 | tr '\n' '\t' >> $7
	cat $1/merge_success_simulated.txt | sed -n '7~7p' | tr '\.' ',' >> $7
}

rm results_score.csv
rm results_switch.csv
rm results_switch_normalized.csv
rm results_mergesuccess.csv
rm results_simulated_score.csv
rm results_simulated_switch.csv
rm results_simulated_mergesuccess.csv
rm results_simulated_switch_normalized.csv
touch results_score.csv
touch results_switch.csv
touch results_switch_normalized.csv
touch results_mergesuccess.csv
touch results_simulated_score.csv
touch results_simulated_switch.csv
touch results_simulated_mergesuccess.csv
touch results_simulated_switch_normalized.csv
parse 100 results_score.csv results_switch.csv results_mergesuccess.csv results_simulated_score.csv results_simulated_switch.csv results_simulated_mergesuccess.csv results_simulated_switch_normalized.csv results_switch_normalized.csv
parse 200 results_score.csv results_switch.csv results_mergesuccess.csv results_simulated_score.csv results_simulated_switch.csv results_simulated_mergesuccess.csv results_simulated_switch_normalized.csv results_switch_normalized.csv
parse 300 results_score.csv results_switch.csv results_mergesuccess.csv results_simulated_score.csv results_simulated_switch.csv results_simulated_mergesuccess.csv results_simulated_switch_normalized.csv results_switch_normalized.csv
parse 400 results_score.csv results_switch.csv results_mergesuccess.csv results_simulated_score.csv results_simulated_switch.csv results_simulated_mergesuccess.csv results_simulated_switch_normalized.csv results_switch_normalized.csv
parse 500 results_score.csv results_switch.csv results_mergesuccess.csv results_simulated_score.csv results_simulated_switch.csv results_simulated_mergesuccess.csv results_simulated_switch_normalized.csv results_switch_normalized.csv
parse 600 results_score.csv results_switch.csv results_mergesuccess.csv results_simulated_score.csv results_simulated_switch.csv results_simulated_mergesuccess.csv results_simulated_switch_normalized.csv results_switch_normalized.csv
parse 700 results_score.csv results_switch.csv results_mergesuccess.csv results_simulated_score.csv results_simulated_switch.csv results_simulated_mergesuccess.csv results_simulated_switch_normalized.csv results_switch_normalized.csv
parse 800 results_score.csv results_switch.csv results_mergesuccess.csv results_simulated_score.csv results_simulated_switch.csv results_simulated_mergesuccess.csv results_simulated_switch_normalized.csv results_switch_normalized.csv
parse 900 results_score.csv results_switch.csv results_mergesuccess.csv results_simulated_score.csv results_simulated_switch.csv results_simulated_mergesuccess.csv results_simulated_switch_normalized.csv results_switch_normalized.csv
parse 1000 results_score.csv results_switch.csv results_mergesuccess.csv results_simulated_score.csv results_simulated_switch.csv results_simulated_mergesuccess.csv results_simulated_switch_normalized.csv results_switch_normalized.csv
parse 1200 results_score.csv results_switch.csv results_mergesuccess.csv results_simulated_score.csv results_simulated_switch.csv results_simulated_mergesuccess.csv results_simulated_switch_normalized.csv results_switch_normalized.csv
parse 1400 results_score.csv results_switch.csv results_mergesuccess.csv results_simulated_score.csv results_simulated_switch.csv results_simulated_mergesuccess.csv results_simulated_switch_normalized.csv results_switch_normalized.csv
parse 1600 results_score.csv results_switch.csv results_mergesuccess.csv results_simulated_score.csv results_simulated_switch.csv results_simulated_mergesuccess.csv results_simulated_switch_normalized.csv results_switch_normalized.csv
parse 1800 results_score.csv results_switch.csv results_mergesuccess.csv results_simulated_score.csv results_simulated_switch.csv results_simulated_mergesuccess.csv results_simulated_switch_normalized.csv results_switch_normalized.csv
parse 2000 results_score.csv results_switch.csv results_mergesuccess.csv results_simulated_score.csv results_simulated_switch.csv results_simulated_mergesuccess.csv results_simulated_switch_normalized.csv results_switch_normalized.csv
