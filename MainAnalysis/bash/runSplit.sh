#!/bin/bash

DATE=`date +%Y%m%d`

topPath1=/home/cfmcginn/Samples/FullJR/ak
topPath2=PFJetAnalyzer/20190116/

rVals=(3 4)

MCDijet=(MCDijet30 MCDijet80 MCDijet170 MCDijet280 MCDijet370 MCDijet540)
splitLevel=(1 2 8 4 4 2)

mkdir -p logs/$DATE

for r in "${rVals[@]}"
do	
    rm -f $topPath1$r$topPath2/*SPLIT*.root
done
    
pos=0
for mc in "${MCDijet[@]}"
do   
    for r in "${rVals[@]}"
    do		
	for i in $topPath1$r$topPath2/*$mc*
	do
	    if [[ $i == *"SPLIT"* ]]
	    then
		continue
	    fi

	    while [[ $i == *"//"* ]]
	    do
		i=$(echo $i | sed -r "s@//@/@g")
	    done
	    	    
	    fileName=${i%*.root}
	    while [[ $fileName == *"/"* ]]
	    do
		fileName=${fileName#*/}	       
	    done
	    
	    ./bin/splitFiles.exe $i ${splitLevel[$pos]} >& logs/$DATE/$fileName.log &

	    counts=$(ps | grep splitF | wc -l)
	    while [[ $counts -ge 10 ]]
	    do
		sleep 5
		counts=$(ps | grep splitF | wc -l)
	    done
	done
    done
   
    pos=$((pos + 1))
done

wait

rm -f paths/*SPLIT*.txt
rm -f paths/*/*SPLIT*.txt

for mc in "${MCDijet[@]}"
do   
    for r in "${rVals[@]}"
    do		
	for i in $topPath1$r$topPath2/*$mc*
	do
	    if [[ $i == *"SPLIT"* ]]
	    then
		dummyVal=0
	    else
		continue
	    fi

	    while [[ $i == *"//"* ]]
	    do
		i=$(echo $i | sed -r "s@//@/@g")
	    done

	    subFile=${i%_SPLIT*}
	    appendStr=${i#$subFile}
	    appendStr=${appendStr%.root}
	    subFile=$subFile.root

	    txtFile=$(grep -nr $subFile paths/)
	    while [[ $txtFile == *":"* ]]
	    do
		txtFile=${txtFile%:*}
	    done
	    txtFile=${txtFile%.txt}

	    while read line;
	    do
		if [[ $line == *"/"* ]]
		then
		    continue
		fi

		echo $line >> $txtFile$appendStr.txt		
	    done < $txtFile.txt

	    echo $i >> $txtFile$appendStr.txt
	done
    done
done


echo "bash runSplit.sh complete!"
