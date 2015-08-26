#!/bin/sh

runTest(){
	echo $2 | tr '\n' '\t' >> $1/$3
	cat $1/result_score_$2_* | sed -n '1~5p' | tr '\n' '\t' | tr '\.' ',' >> $1/$3
	echo >> $1/$3
	echo $2 | tr '\n' '\t' >> $1/$4
	cat $1/result_switch_$2_* | sed -n '10~10p' | grep -Eo '^[^ ]+' | tr '\n' '\t' >> $1/$4
	echo >> $1/$4
}

runTest2(){
	rm $1/results_score.csv
	rm $1/results_switch.csv
	touch $1/results_score.csv
	touch $1/results_switch.csv
	runTest $1 0 results_score.csv results_switch.csv
	runTest $1 1 results_score.csv results_switch.csv
	runTest $1 2 results_score.csv results_switch.csv
	runTest $1 3 results_score.csv results_switch.csv
	runTest $1 4 results_score.csv results_switch.csv
	runTest $1 5 results_score.csv results_switch.csv
	runTest $1 6 results_score.csv results_switch.csv
	runTest $1 7 results_score.csv results_switch.csv
	runTest $1 8 results_score.csv results_switch.csv
	runTest $1 9 results_score.csv results_switch.csv
	runTest $1 10 results_score.csv results_switch.csv
	runTest $1 11 results_score.csv results_switch.csv
	runTest $1 12 results_score.csv results_switch.csv
	runTest $1 13 results_score.csv results_switch.csv
	runTest $1 14 results_score.csv results_switch.csv
	runTest $1 15 results_score.csv results_switch.csv
}

runTest2 simul_500
runTest2 simul_700
runTest2 simul_1000
