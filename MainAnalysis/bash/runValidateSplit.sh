#!/bin/bash

DATE=`date +%Y%m%d`
mkdir -p logs/$DATE

rVals=(3 4 6 8 10)
MCDijet=(MCDijet30 MCDijet80 MCDijet170 MCDijet280 MCDijet370 MCDijet540)

topPath=(/home/cfmcginn/Samples/FullJR/akCs /home/cfmcginn/Samples/FullJR/ak)

for i in "${topPath[@]}"
do
    pbpbOrPPTag=${i#*/}   
    while [[ $pbpbOrPPTag == *"/"* ]]
    do
	pbpbOrPPTag=${pbpbOrPPTag#*/}   	
    done
    
    for r in "${rVals[@]}"
    do
	counts=$(ls $i$r*/*/*.root | wc -l)

	if [[ $counts -gt 0 ]]
	then
	    for mc in "${MCDijet[@]}"
	    do
		mainFile=""
		splitFiles=""

		for l in $i$r*/*/*$mc*.root
		do
		    if [[ $l == *"SPLIT"* ]]
		    then
			splitFiles=$splitFiles$l,
		    else
			mainFile=$l
		    fi		    
		done

		
		./bin/validateSplit.exe $mainFile $splitFiles >& logs/$DATE/validateSplit_"$mc"_"$pbpbOrPPTag$r".log &

		counts=$(ps | grep validateS | wc -l)

		while [[ $counts -ge 10 ]]
		do
		    sleep 10
		    counts=$(ps | grep validateS | wc -l)
		done

	    done
	fi
    done	     
done
