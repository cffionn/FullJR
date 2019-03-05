#!/bin/bash

DATE=`date +%Y%m%d`

nSplits=10
if [[ -d /data/cmcginn ]]
then
    nSplits=20
fi

topPath1PbPb=/data/cmcginn/Forests/Pythia6HydjetDijet/LargeConeRAA_20180601/akCs
topPath1PP=/data/cmcginn/Forests/Pythia6Dijet/LargeConeRAA_20180712/ak

topPath2PbPb=PU3PFFlowJetAnalyzer/20190116/
topPath2PP=PFJetAnalyzer/20190116/

#rVals=(3)
rVals=(3 4 6 8 10)

MCDijet=(MCDijet30 MCDijet80 MCDijet170 MCDijet280 MCDijet370 MCDijet540)
#MCDijet=(MCDijet540)
splitLevel=(1 2 8 8 8 4)
#splitLevel=(2)

mkdir -p logs/$DATE

#for r in "${rVals[@]}"
#do	
#    rm -f $topPath1PbPb$r$topPath2PbPb/*SPLIT*.root
#    rm -f $topPath1PP$r$topPath2PP/*SPLIT*.root
#done
#    
#pos=0
#for mc in "${MCDijet[@]}"
#do   
#    for r in "${rVals[@]}"
#    do		
#	for i in $topPath1PbPb$r$topPath2PbPb/*$mc*
#	do
#	    if [[ $i == *"SPLIT"* ]]
#	    then
#		continue
#	    fi
#
#	    while [[ $i == *"//"* ]]
#	    do
#		i=$(echo $i | sed -r "s@//@/@g")
#	    done
#	    	    
#	    fileName=${i%*.root}
#	    while [[ $fileName == *"/"* ]]
#	    do
#		fileName=${fileName#*/}	       
#	    done
#	    
#	    ./bin/splitFiles.exe $i ${splitLevel[$pos]} >& logs/$DATE/$fileName.log &
#
#	    counts=$(ps | grep splitF | wc -l)
#	    while [[ $counts -ge $nSplits ]]
#	    do
#		sleep 5
#		counts=$(ps | grep splitF | wc -l)
#	    done
#	done
#
#	for i in $topPath1PP$r$topPath2PP/*$mc*
#	do
#	    if [[ $i == *"SPLIT"* ]]
#	    then
#		continue
#	    fi
#
#	    while [[ $i == *"//"* ]]
#	    do
#		i=$(echo $i | sed -r "s@//@/@g")
#	    done
#	    	    
#	    fileName=${i%*.root}
#	    while [[ $fileName == *"/"* ]]
#	    do
#		fileName=${fileName#*/}	       
#	    done
#	    
#	    ./bin/splitFiles.exe $i ${splitLevel[$pos]} >& logs/$DATE/$fileName.log &
#
#	    counts=$(ps | grep splitF | wc -l)
#	    while [[ $counts -ge $nSplits ]]
#	    do
#		sleep 5
#		counts=$(ps | grep splitF | wc -l)
#	    done
#	done
#    done
#   
#    pos=$((pos + 1))
#done
#
#wait

rm -f paths/*SPLIT*.txt
rm -f paths/*/*SPLIT*.txt

for mc in "${MCDijet[@]}"
do   
    for r in "${rVals[@]}"
    do		
	for i in $topPath1PbPb$r$topPath2PbPb/*$mc*
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

		echo $line >> $txtFile$mc$appendStr.txt		
	    done < $txtFile.txt

	    echo $i >> $txtFile$mc$appendStr.txt
	done

	for i in $topPath1PP$r$topPath2PP/*$mc*
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

		echo $line >> $txtFile$mc$appendStr.txt		
	    done < $txtFile.txt

	    echo $i >> $txtFile$mc$appendStr.txt
	done

    done
done


echo "bash runSplit.sh complete!"
