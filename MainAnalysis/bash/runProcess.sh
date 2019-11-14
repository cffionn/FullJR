#!/bin/bash

runDir=$FULLJRDIR/MainAnalysis
curDir=$PWD

if [[ "$runDir" == "$curDir" ]]
then
    dummyVal=0
else
    echo "Please run from $PWD, as bash bash/runProcess.sh. exit 1"
    exit 1
fi

if [[ -f $runDir/bin/processRawData.exe ]]
then
    dummyVal=0
else
    echo "Necessary executable ./bin/processRawData.exe is missing. run make. exit 1"
    exit 1
fi


dateStrPP=20190604
dateStrPbPb=20190604

DATE=`date +%Y%m%d`

mkdir -p logs
mkdir -p logs/$DATE


jetsPP=(ak2PFJetAnalyzer ak3PFJetAnalyzer ak4PFJetAnalyzer ak6PFJetAnalyzer ak8PFJetAnalyzer ak10PFJetAnalyzer)
jetsPbPb=(akCs2PU3PFFlowJetAnalyzer akCs3PU3PFFlowJetAnalyzer akCs4PU3PFFlowJetAnalyzer akCs6PU3PFFlowJetAnalyzer akCs8PU3PFFlowJetAnalyzer akCs10PU3PFFlowJetAnalyzer)

#jetsPP=(ak8PFJetAnalyzer ak10PFJetAnalyzer)
#jetsPbPb=(akCs8PU3PFFlowJetAnalyzer akCs10PU3PFFlowJetAnalyzer)

#PbPbDataFile=/data/cmcginn/Forests/PbPb2015Data/HIHardProbes/HiForestAOD_HIHardProbes_HLTJet100_AllR_PtCut140_AbsEta3_20180626_21LumiPer_180626_152510_1050_OutOf1050_MERGED.root
PbPbDataFile=/data/cmcginn/Forests/PbPb2015Data/HIHardProbes/HiForestAOD_HIHardProbes_HLTJet100_AllR_PtCut140_AbsEta3_20190422_21LumiPer_190422_210102_MERGED_1049FilesOutOf1050_MERGED_TTreeSkim_nEvtAll_nEvtStart0_NoRLESkim_NoCut_20190501
PbPbMCFilePre=output/"$dateStrPbPb"/combinedResponse_
PbPbMCFilePost=_"$dateStrPbPb".root

#BAD JSON
#PPDataFile=/data/cmcginn/Forests/pp2015Data/HIHardProbes/HiForestAOD_HighPtJet80_HLTJet80_LargeROR_PtCut110_AbsEta5_20180627_11Lumi_180627_122059_355_OutOf355_MERGED.root
#PPDataFile=/data/cmcginn/Forests/pp2015Data/HIHardProbes/HiForestAOD_HighPtJet80_HLTJet80_LargeROR_PtCut110_AbsEta5_20190220_11Lumi_190220_221659_561_OutOf561_MERGED
PPDataFile=/data/cmcginn/Forests/pp2015Data/HIHardProbes/HiForestAOD_HighPtJet80_HLTJet80_LargeROR_PtCut110_AbsEta5_20190220_11Lumi_190220_221659_561_OutOf561_MERGED_TTreeSkim_nEvtAll_nEvtStart0_NoRLESkim_NoCut_20190524
PPMCFilePre=output/"$dateStrPP"/combinedResponse_
PPMCFilePost=_"$dateStrPP".root


pos=0
for i in "${jetsPP[@]}"
do
    echo "Processing $i, ${jetsPbPb[$pos]}..."

    newPPDataFile="$PPDataFile"_$i.root
    if [[ -f $newPPDataFile ]]
    then
	dummyVal=0
    else
	echo "$newPPDataFile is not found, please update. exit 1"
	exit 1
    fi
    
    fileNamePP=$PPMCFilePre${i%PF*}$PPMCFilePost
    
    if [[ -f $fileNamePP ]]
    then
	./bin/processRawData.exe $newPPDataFile $fileNamePP 1 $i >& logs/$DATE/processPP_$i.log &
    else
	echo " $fileNamePP is not found, continue on $i"
    fi

    fileNamePbPb=${jetsPbPb[$pos]}
    fileNamePbPb=${fileNamePbPb%PU*}
    fileNamePbPb=$PbPbMCFilePre$fileNamePbPb$PbPbMCFilePost

    newPbPbDataFile="$PbPbDataFile"_${jetsPbPb[$pos]}.root
    if [[ -f $newPbPbDataFile ]]
    then
	dummyVal=0
    else
	echo "$newPbPbDataFile is not found, please update. exit 1"
	exit 1
    fi


    if [[ -f $fileNamePbPb ]]
    then
	dummyVal=0
	./bin/processRawData.exe $newPbPbDataFile $fileNamePbPb 0 ${jetsPbPb[$pos]} >& logs/$DATE/processPbPb_${jetsPbPb[$pos]}.log &
    else
	echo " $fileNamePbPb is not found, continue on ${jetsPbPb[$pos]}"
    fi

    pos=$((pos+1))
##    wait
done

wait

echo "runProcess.sh Complete"
