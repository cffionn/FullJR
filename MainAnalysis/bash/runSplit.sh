#!/bin/bash

DATE=`date +%Y%m%d`

nSplits=10
if [[ -d /data/cmcginn ]]
then
    nSplits=20
fi

#topPath1PbPb=/data/cmcginn/Forests/Pythia6HydjetDijet/LargeConeRAA_20180601/akCs
topPath1PP=/data/cmcginn/Forests/Pythia6Dijet/LargeConeRAA_20180712/ak

if [[ -d /data/cmcginn ]]
then
    dummyVal=0
else
#    topPath1PbPb=/home/cfmcginn/Samples/FullJR/akCs
    topPath1PP=/home/cfmcginn/Samples/FullJR/ak
fi

#topPath2PbPb=PU3PFFlowJetAnalyzer/20190116/
topPath2PP=PFJetAnalyzer/20190116/

#rVals=(2)
rVals=(2 3 4 6 8 10)

MCDijet=(MCDijet30 MCDijet80 MCDijet170 MCDijet280 MCDijet370 MCDijet540)
#MCDijet=(MCDijet540)
splitLevel=(1 2 8 8 8 4)
#splitLevel=(2)

mkdir -p logs/$DATE

for r in "${rVals[@]}"
do	
#    rm -f $topPath1PbPb$r$topPath2PbPb/*SPLIT*.root
    rm -f $topPath1PP$r$topPath2PP/*SPLIT*.root
done
    
pos=0
for mc in "${MCDijet[@]}"
do   
    for r in "${rVals[@]}"
    do

	echo $mc $r

#	counts=$(ls $topPath1PbPb$r$topPath2PbPb/*$mc* | wc -l)

#	if [[ $counts -ge 1 ]]
#	then
#	    for i in $topPath1PbPb$r$topPath2PbPb/*$mc*
#	    do
#		if [[ $i == *"SPLIT"* ]]
#		then
#		    continue
#		fi
#		
#		while [[ $i == *"//"* ]]
#		do
#		    i=$(echo $i | sed -r "s@//@/@g")
#		done
#	    	
#		fileName=${i%*.root}
#		while [[ $fileName == *"/"* ]]
#		do
#		    fileName=${fileName#*/}	       
#		done
#		
#		./bin/splitFiles.exe $i ${splitLevel[$pos]} >& logs/$DATE/$fileName.log &
#		
#		counts=$(ps | grep splitF | wc -l)
#		while [[ $counts -ge $nSplits ]]
#		do
#		    sleep 5
#		    counts=$(ps | grep splitF | wc -l)
#		done
#	    done
#	fi
	
	counts=$(ls $topPath1PP$r$topPath2PP/*$mc* | wc -l)
	if [[ $counts -ge 1 ]]
	then
	    for i in $topPath1PP$r$topPath2PP/*$mc*
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
		while [[ $counts -ge $nSplits ]]
		do
		    sleep 5
		    counts=$(ps | grep splitF | wc -l)
		done
	    done
	fi
    done
   
    pos=$((pos + 1))
done

wait

#rm -f paths/*SPLIT*.txt
#rm -f paths/*/*SPLIT*.txt

for mc in "${MCDijet[@]}"
do   
    for r in "${rVals[@]}"
    do		

#	counts=$(ls $topPath1PbPb$r$topPath2PbPb/*$mc* | wc -l)
#	if [[ $counts -ge 1 ]]
#	then	    
#	    for i in $topPath1PbPb$r$topPath2PbPb/*$mc*
#	    do
#		if [[ $i == *"SPLIT"* ]]
#		then
#		    dummyVal=0
#		else
#		    continue
#		fi
#		
#		while [[ $i == *"//"* ]]
#		do
#		    i=$(echo $i | sed -r "s@//@/@g")
#		done
#		
#		subFile=${i%_SPLIT*}
#		appendStr=${i#$subFile}
#		appendStr=${appendStr%.root}
#		subFile=$subFile.root
#		
#		txtFile=$(grep -nr $subFile paths/)
#		while [[ $txtFile == *":"* ]]
#		do
#		    txtFile=${txtFile%:*}
#		done
#		txtFile=${txtFile%.txt}
#		
#		while read line;
#		do
#		    if [[ $line == *"/"* ]]
#		    then
#			continue
#		    fi
#		    
#		    echo $line >> $txtFile$mc$appendStr.txt		
#		done < $txtFile.txt
#		
#		echo $i >> $txtFile$mc$appendStr.txt
#	    done
#	fi

#	echo $r $mc
	counts=$(ls $topPath1PP$r$topPath2PP/*$mc* | wc -l)
	echo " COuNTS $counts"
#	echo $(ls $topPath1PP$r$topPath2PP/*$mc*)


	if [[ $counts -ge 1 ]]
	then	    	    
	    for i in $topPath1PP$r$topPath2PP/*$mc*
	    do
#		echo $i

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

#		echo "APPEND SUB " $appendStr $subFile
		
		txtFile=$(grep -nr $subFile paths/)
		
		while [[ $txtFile == *":"* ]]
		do
		    txtFile=${txtFile%:*}
		done
		txtFile=${txtFile%.txt}

		echo "TXT " $txtFile
		
		if [[ -f $txtFile$mc$appendStr.txt ]]
		then
		    rm $txtFile$mc$appendStr.txt
		fi

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
	fi
    done
done


echo "bash runSplit.sh complete!"
