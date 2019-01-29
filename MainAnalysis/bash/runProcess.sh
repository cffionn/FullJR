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

dateStrPP=20190129
dateStrPbPb=20190129

DATE=`date +%Y%m%d`

mkdir -p logs
mkdir -p logs/$DATE


jetsPP=(ak3PFJetAnalyzer ak4PFJetAnalyzer ak6PFJetAnalyzer ak8PFJetAnalyzer ak10PFJetAnalyzer)
jetsPbPb=(akCs3PU3PFFlowJetAnalyzer akCs4PU3PFFlowJetAnalyzer akCs6PU3PFFlowJetAnalyzer akCs8PU3PFFlowJetAnalyzer akCs10PU3PFFlowJetAnalyzer)

PbPbDataFile=/data/cmcginn/Forests/PbPb2015Data/HIHardProbes/HiForestAOD_HIHardProbes_HLTJet100_AllR_PtCut140_AbsEta3_20180626_21LumiPer_180626_152510_1050_OutOf1050_MERGED.root
PbPbMCFilePre=output/"$dateStrPbPb"/Pythia6_Dijet_pp502_Hydjet_Cymbal_MB_PbPb_MCDijet_20180521_ExcludeTop4_ExcludeToFrac_Frac0p7_Full_5Sigma_20180608_SVM_
PbPbMCFilePost=_FracNEntries1p00_JetResponse_"$dateStrPbPb".root

PPDataFile=/data/cmcginn/Forests/pp2015Data/HIHardProbes/HiForestAOD_HighPtJet80_HLTJet80_LargeROR_PtCut110_AbsEta5_20180627_11Lumi_180627_122059_355_OutOf355_MERGED.root
PPMCFilePre=output/"$dateStrPP"/Pythia6_Dijet_pp502_MCDijet_20180712_ExcludeTop4_ExcludeToFrac_Frac0p7_Full_5Sigma_20180712_SVM_
PPMCFilePost=_FracNEntries1p00_JetResponse_"$dateStrPP".root

if [[ -f $PbPbDataFile ]]
then
    dummyVal=0
else
    echo "$PbPbDataFile is not found, please update. exit 1"
    exit 1
fi

if [[ -f $PPDataFile ]]
then
    dummyVal=0
else
    echo "$PPDataFile is not found, please update. exit 1"
    exit 1
fi

pos=0
for i in "${jetsPP[@]}"
do
    echo "Processing $i, ${jetsPbPb[$pos]}..."

    fileNamePP=$PPMCFilePre$i$PPMCFilePost
    
    if [[ -f $fileNamePP ]]
    then
	./bin/processRawData.exe $PPDataFile $fileNamePP 1 $i >& logs/$DATE/processPP_$i.log &
    else
	echo " $fileNamePP is not found, continue on $i"
    fi

    fileNamePbPb=$PbPbMCFilePre${jetsPbPb[$pos]}$PbPbMCFilePost

    if [[ -f $fileNamePbPb ]]
    then
	./bin/processRawData.exe $PbPbDataFile $fileNamePbPb 0 ${jetsPbPb[$pos]} >& logs/$DATE/processPbPb_${jetsPbPb[$pos]}.log &
    else
	echo " $fileNamePbPb is not found, continue on ${jetsPbPb[$pos]}"
    fi

    pos=$((pos+1))
    wait
done

wait

echo "runProcess.sh Complete"
