#!/bin/bash

runDir=$FULLJRDIR/MainAnalysis
curDir=$PWD

if [[ "$runDir" == "$curDir" ]]
then
    dummyVal=0
else
    echo "Please run from $PWD, as bash bash/runResponse.sh. exit 1"
    exit 1
fi

if [[ -f $runDir/bin/makeJetResponseTree.exe ]]
then
    dummyVal=0
else
    echo "Necessary executable ./bin/makeJetResponseTree.exe is missing. run make. exit 1"
    exit 1
fi

jobNumber=16
if [[ -d /home/cfmcginn ]]
then
    jobNumber=10
fi

DATE=`date +%Y%m%d`
TIMESTART=`date +%H%M%S`

mkdir -p logs
mkdir -p logs/$DATE

filesToProcess=()
counts=$(ls paths/*SPLIT*.txt | wc -l)
if [[ $counts -ge 1 ]]
then
    for i in paths/*SPLIT*.txt
    do
	filesToProcess+=($i)
    done
fi    

counts=$(ls paths/*/*SPLIT*.txt | wc -l)
if [[ $counts -ge 1 ]]
then
    for i in paths/*/*SPLIT*.txt
    do
	filesToProcess+=($i)
    done
fi


for i in "${filesToProcess[@]}"
do
    pbpbBool=1
    if [[ $i == *"akCs"* ]]
    then
	pbpbBool=0
    fi

    logName=${i%.txt}
    while [[ $logName == *"/"* ]]
    do
	logName=${logName#*/}
    done

    ./bin/makeJetResponseTree.exe $i $pbpbBool 1.0 >& logs/$DATE/response_$logName.log &

    counts=$(ps | grep makeJet | wc -l)
    while [[ $counts -ge $jobNumber ]]
    do
	sleep 10
	counts=$(ps | grep makeJet | wc -l)
    done
done

wait

DATEEND=`date +%Y%m%d`

rVals=(ak3 ak4 ak6 ak8 ak10 akCs3 akCs4 akCs6 akCs8 akCs10)
rPaths=()

for r in "${rVals[@]}"
do
    rPath="NA"
    for i in "${filesToProcess[@]}"
    do
	val=$(grep $r $i)
	if [[ $val == "/"* ]]
	then
	    echo $val, $r
	    rPath=${val%/*}
	    break
	fi
    done

    rPaths+=($rPath)
done

pos=0
for r in "${rVals[@]}"
do
    echo $r, ${rPaths[$pos]}
    pos=$((pos + 1))
done


pos=0
for r in "${rVals[@]}"
do
    files=""

    countLocalDate=0
    if [[ -d output/$DATE ]]
    then
	countLocalDate=$(ls output/$DATE/*$r*SPLIT*.root | wc -l)
    fi

    countLocalDateEnd=0
    if [[ $DATE -ne $DATEEND ]]
    then
	if [[ -d output/$DATEEND ]]
	then
	    countLocalDateEnd=$(ls output/$DATEEND/*$r*SPLIT*.root | wc -l)
	fi
    fi
    
    countGlobalDate=0
    if [[ -d ${rPaths[$pos]}/output/$DATE ]]
    then
	countGlobalDate=$(ls ${rPaths[$pos]}/output/$DATE/*$r*SPLIT*.root | wc -l)
    fi

    countGlobalDateEnd=0
    if [[ $DATE -ne $DATEEND ]]
    then
	if [[ -d ${rPaths[$pos]}/output/$DATEEND ]]
	then
	    countGlobalDateEnd=$(ls ${rPaths[$pos]}/output/$DATEEND/*$r*SPLIT*.root | wc -l)
	fi
    fi
    count=$((countLocalDate + $countLocalDateEnd + $countGlobalDate + $countGlobalDateEnd))
    

    if [[ $count -eq 0 ]]
    then
	pos=$((pos + 1))
	continue
    fi

    if [[ $countLocalDate -ne 0 ]]
    then
	for i in output/$DATE/*$r*SPLIT*.root
	do
	    files=$files$i,
	done
    fi
    
    if [[ $countLocalDateEnd -ne 0 ]]
    then
	for i in output/$DATEEND/*$r*SPLIT*.root
	do
	    files=$files$i,
	done	
    fi

    if [[ $countGlobalDate -ne 0 ]]
    then
	for i in ${rPaths[$pos]}/output/$DATE/*$r*SPLIT*.root
	do
	    files=$files$i,
	done
    fi
    
    if [[ $countGlobalDateEnd -ne 0 ]]
    then
	for i in ${rPaths[$pos]}/output/$DATEEND/*$r*SPLIT*.root
	do
	    files=$files$i,
	done	
    fi

    ./bin/combineResponse.exe $files $r >& logs/$DATE/combineResponse_$r.log &
    jobCounts=$(ps | grep combineR | wc -l)
    if [[ $jobCounts -ge $jobNumber ]]
    then
	sleep 10
	jobCounts=$(ps | grep combineR | wc -l)
    fi

    pos=$((pos + 1))
done

wait
	     
TIMEEND=`date +%H%M%S`

echo "Time start: $TIMESTART"
echo "Time end: $TIMEEND"

echo "WARNING: DATE START \'$DATE\' not same as DATE END \'$DATEEND\'"

echo "bash/runSplitResponse.sh Complete!"
