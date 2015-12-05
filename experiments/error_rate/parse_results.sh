#!/bin/sh

runTest(){
	echo $2 | tr '\n' '\t' >> $1/$3
	cat $1/result_score_$2_* | sed -n '1~5p' | tr '\n' '\t' | tr '\.' ',' >> $1/$3
	echo >> $1/$3
	
	echo $2 | tr '\n' '\t' >> $1/$4
	cat $1/result_switch_$2_* | sed -n '11~14p' | grep -Eo '^[^ ]+' | tr '\n' '\t' >> $1/$4
	echo >> $1/$4
	
	echo $2 | tr '\n' '\t' >> $1/$6
	cat $1/result_switch_$2_* | sed -n '13~14p' | grep -Eo '^[^ ]+' | tr '\n' '\t' >> $1/$6
	echo >> $1/$6
	
	echo $2 | tr '\n' '\t' >> $1/$5
	cat $1/result_mergesuccess_$2_*.txt | sed -n '7~7p' | tr '\n' '\t' | tr '\.' ',' >> $1/$5
	echo >> $1/$5
}

runTest2(){
	rm $1/results_score.csv
	rm $1/results_switch.csv
	rm $1/results_mergesuccess.csv
	rm $1/results_switch_normalized.csv
	touch $1/results_score.csv
	touch $1/results_switch.csv
	touch $1/results_mergesuccess.csv
	touch $1/results_switch_normalized.csv
	runTest $1 0 results_score.csv results_switch.csv results_mergesuccess.csv results_switch_normalized.csv
	runTest $1 1 results_score.csv results_switch.csv results_mergesuccess.csv results_switch_normalized.csv
	runTest $1 2 results_score.csv results_switch.csv results_mergesuccess.csv results_switch_normalized.csv
	runTest $1 3 results_score.csv results_switch.csv results_mergesuccess.csv results_switch_normalized.csv
	runTest $1 4 results_score.csv results_switch.csv results_mergesuccess.csv results_switch_normalized.csv
	runTest $1 5 results_score.csv results_switch.csv results_mergesuccess.csv results_switch_normalized.csv
	runTest $1 6 results_score.csv results_switch.csv results_mergesuccess.csv results_switch_normalized.csv
	runTest $1 7 results_score.csv results_switch.csv results_mergesuccess.csv results_switch_normalized.csv
	runTest $1 8 results_score.csv results_switch.csv results_mergesuccess.csv results_switch_normalized.csv
	runTest $1 9 results_score.csv results_switch.csv results_mergesuccess.csv results_switch_normalized.csv
	runTest $1 10 results_score.csv results_switch.csv results_mergesuccess.csv results_switch_normalized.csv
	runTest $1 11 results_score.csv results_switch.csv results_mergesuccess.csv results_switch_normalized.csv
	runTest $1 12 results_score.csv results_switch.csv results_mergesuccess.csv results_switch_normalized.csv
	runTest $1 13 results_score.csv results_switch.csv results_mergesuccess.csv results_switch_normalized.csv
	runTest $1 14 results_score.csv results_switch.csv results_mergesuccess.csv results_switch_normalized.csv
	runTest $1 15 results_score.csv results_switch.csv results_mergesuccess.csv results_switch_normalized.csv
}

runTest2 simul_500
runTest2 simul_700
runTest2 simul_1000
runTest2 discosnp_1000
runTest2 discosnp_1600
runTest2 discosnp_2000