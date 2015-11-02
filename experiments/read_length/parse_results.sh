#!/bin/sh

parse(){
	echo $1 | tr '\n' '\t' >> $2
	cat $1/result_score.txt | sed -n '1~5p' | tr '\.' '\,' >> $2

	echo $1 | tr '\n' '\t' >> $3
	cat $1/result_switch.txt | sed -n '10~10p' | tr ' ' '\t' >> $3

	echo $1 | tr '\n' '\t' >> $4
	cat $1/merge_success_similars2.txt | sed -n '7~7p' | tr '\.' ',' >> $4
}

rm results_score.csv
rm results_switch.csv
rm results_mergesuccess.csv
touch results_score.csv
touch results_switch.csv
touch results_mergesuccess.csv
parse 100 results_score.csv results_switch.csv results_mergesuccess.csv
parse 200 results_score.csv results_switch.csv results_mergesuccess.csv
parse 300 results_score.csv results_switch.csv results_mergesuccess.csv
parse 400 results_score.csv results_switch.csv results_mergesuccess.csv
parse 500 results_score.csv results_switch.csv results_mergesuccess.csv
parse 600 results_score.csv results_switch.csv results_mergesuccess.csv
parse 700 results_score.csv results_switch.csv results_mergesuccess.csv
parse 800 results_score.csv results_switch.csv results_mergesuccess.csv
parse 900 results_score.csv results_switch.csv results_mergesuccess.csv
parse 1000 results_score.csv results_switch.csv results_mergesuccess.csv
parse 1200 results_score.csv results_switch.csv results_mergesuccess.csv
parse 1400 results_score.csv results_switch.csv results_mergesuccess.csv
parse 1600 results_score.csv results_switch.csv results_mergesuccess.csv
parse 1800 results_score.csv results_switch.csv results_mergesuccess.csv
parse 2000 results_score.csv results_switch.csv results_mergesuccess.csv
