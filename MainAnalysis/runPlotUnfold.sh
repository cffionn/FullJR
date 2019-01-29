#!/bin/bash

pbpbVals=(akCs3PU3PFFlow akCs4PU3PFFlow akCs6PU3PFFlow akCs8PU3PFFlow akCs10PU3PFFlow)
ppVals=(ak3PF ak4PF ak6PF ak8PF ak10PF)

fileOutPbPbPre=output/HiForestAOD_HIHardProbes_HLTJet100_AllR_Pythia6_Dijet_pp502_Hydjet_Cymbal_MB_PbP_UnfoldRawData_NSuperBayes0_
fileOutPbPbPost=_20180913.root

fileOutPPPre=output/HiForestAOD_HighPtJet80_HLTJet80_LargeRO_Pythia6_Dijet_pp502_MCDijet_20180712_Exc_UnfoldRawData_NSuperBayes0_
fileOutPPPost=_20180913.root

pos=0
for i in "${pbpbVals[@]}"
do
    echo "Processing $i, ${ppVals[$pos]}"
    ./bin/plotUnfoldedSpectra.exe $fileOutPPPre${ppVals[$pos]}$fileOutPPPost $fileOutPbPbPre$i$fileOutPbPbPost >& logs/plotUnfold_${ppVals[$pos]}_$i.log &
    pos=$((pos + 1))
    
done

wait

echo "runPlotUnfold.sh complete!"