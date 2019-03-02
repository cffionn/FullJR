//cpp dependencies
#include <iostream>
#include <map>
#include <string>
#include <vector>

//ROOT dependencies
#include "TFile.h"
#include "TTree.h"

//Non-local dependencies
#include "Utility/include/fileUtilities.h"
#include "Utility/include/returnRootFileContentsList.h"

//Intended for cross check against ntuples from yeonju - reduced content ttrees

ULong64_t getRunLumiKey(ULong64_t run, ULong64_t lumi){return run*100000 + lumi;}
ULong64_t runLumiKeyToRun(ULong64_t key){return key/100000;}
ULong64_t runLumiKeyToLumi(ULong64_t key){return key%100000;}

int checkOverlap(const std::string fileName1, const std::string subsetTrig1, const std::string fileName2)
{
  if(!fileIsGood(fileName1, ".root")) return 1;
  if(!fileIsGood(fileName2, ".root")) return 1;

  const std::string hltTreeStr = "hltanalysis/HltTree";
  const std::string hltTreeAltStr = "hltTree";
  const std::string hltTreeAlt2Str = "hltanalysisMOD/HltTree";
  const std::string hiTreeStr = "hiEvtAnalyzer/HiTree";
  const std::string hiTreeAltStr = "HiEvt";

  bool file1HasHLTTree = false;
  bool file1HasHLTTreeAlt = false;
  bool file1HasHLTTreeAlt2 = false;

  bool file2HasHITree = false;
  bool file2HasHITreeAlt = false;

  TFile* inFile_p = new TFile(fileName1.c_str(), "READ");
  std::vector<std::string> treeList1 = returnRootFileContentsList(inFile_p, "TTree", "");
  inFile_p->Close();
  delete inFile_p;

  inFile_p = new TFile(fileName2.c_str(), "READ");
  std::vector<std::string> treeList2 = returnRootFileContentsList(inFile_p, "TTree", "");
  inFile_p->Close();
  delete inFile_p;

  for(auto const & tree : treeList1){
    if(isStrSame(tree, hltTreeStr)) file1HasHLTTree = true;
    else if(isStrSame(tree, hltTreeAltStr)) file1HasHLTTreeAlt = true;
    else if(isStrSame(tree, hltTreeAlt2Str)) file1HasHLTTreeAlt2 = true;
  }
  for(auto const & tree : treeList2){
    if(isStrSame(tree, hiTreeAltStr)) file2HasHITreeAlt = true;
    else if(isStrSame(tree, hiTreeAltStr)) file2HasHITreeAlt = true;
  }

  if(!file1HasHLTTree && !file1HasHLTTreeAlt && !file1HasHLTTreeAlt2){
    std::cout << "File \'" << fileName1 << "\' doesn't contain hlt tree. return 1" << std::endl;
    return 1;
  }

  if(!file2HasHITree && !file2HasHITreeAlt){
    std::cout << "File \'" << fileName2 << "\' doesn't contain hi tree. return 1" << std::endl;
    return 1;
  }

  std::string hltTreeStrFinal = "";
  if(file1HasHLTTree) hltTreeStrFinal = hltTreeStr;
  else if(file1HasHLTTreeAlt) hltTreeStrFinal = hltTreeAltStr;
  else if(file1HasHLTTreeAlt2) hltTreeStrFinal = hltTreeAlt2Str;

  Int_t RunHLT, LumiBlockHLT, trigger;
  Long64_t EventHLT;

  std::map<ULong64_t, std::vector<ULong64_t> > runLumiMapToEvt;
  inFile_p = new TFile(fileName1.c_str(), "READ");
  TTree* hltTree_p = (TTree*)inFile_p->Get(hltTreeStrFinal.c_str());
  hltTree_p->SetBranchStatus("*", 0);
  hltTree_p->SetBranchStatus("Run", 1);
  hltTree_p->SetBranchStatus("LumiBlock", 1);
  hltTree_p->SetBranchStatus("Event", 1);
  hltTree_p->SetBranchStatus(subsetTrig1.c_str(), 1);

  hltTree_p->SetBranchAddress("Run", &RunHLT);
  hltTree_p->SetBranchAddress("LumiBlock", &LumiBlockHLT);
  hltTree_p->SetBranchAddress("Event", &EventHLT);
  hltTree_p->SetBranchAddress(subsetTrig1.c_str(), &trigger);
  
  const ULong64_t nEntries1 = hltTree_p->GetEntries();
  const ULong64_t nDiv1 = TMath::Max(nEntries1/20, (ULong64_t)1);
  std::cout << "Processing file \'" << fileName1 << "\'..." << std::endl;
  for(ULong64_t entry = 0; entry < nEntries1; ++entry){
    if(entry%nDiv1 == 0) std::cout << " Entry " << entry << "/" << nEntries1 << std::endl;
    hltTree_p->GetEntry(entry);
    if(trigger == 0) continue;
    
    ULong64_t key = getRunLumiKey(RunHLT, LumiBlockHLT);
    if(runLumiMapToEvt.count(key) == 0) runLumiMapToEvt[key] = {};
    runLumiMapToEvt[key].push_back(EventHLT);
  }

  inFile_p->Close();
  delete inFile_p;

  std::string hiTreeStrFinal = "";
  if(file2HasHITree) hiTreeStrFinal = hiTreeStr;
  else if(file2HasHITreeAlt) hiTreeStrFinal = hiTreeAltStr;
  
  UInt_t runHI, lumiHI;
  ULong64_t evtHI;

  inFile_p = new TFile(fileName2.c_str(), "READ");
  TTree* hiTree_p = (TTree*)inFile_p->Get(hiTreeStrFinal.c_str());
  hiTree_p->SetBranchStatus("*", 0);
  hiTree_p->SetBranchStatus("run", 1);
  hiTree_p->SetBranchStatus("lumi", 1);
  hiTree_p->SetBranchStatus("evt", 1);

  hiTree_p->SetBranchAddress("run", &runHI);
  hiTree_p->SetBranchAddress("lumi", &lumiHI);
  hiTree_p->SetBranchAddress("evt", &evtHI);

  const ULong64_t nEntries2 = hiTree_p->GetEntries();
  const ULong64_t nDiv2 = TMath::Max(nEntries2/20, (ULong64_t)1);
  std::cout << "Processing file \'" << fileName2 << "\'..." << std::endl;

  for(ULong64_t entry = 0; entry < nEntries2; ++entry){
    if(entry%nDiv2 == 0) std::cout << " Entry " << entry << "/" << nEntries2 << std::endl;
    hiTree_p->GetEntry(entry);

    ULong64_t key = getRunLumiKey(runHI, lumiHI);
    std::vector<ULong64_t>* tempVect = &(runLumiMapToEvt[key]);

    bool isFound = false;
    for(auto const & evt : (*tempVect)){
      if(evt == evtHI){
	isFound = true;
	break;
      }
    }

    if(!isFound) std::cout << "EVENT NOT FOUND: " << runHI << ", " << lumiHI << ", " << evtHI << std::endl;
  }

  inFile_p->Close();
  delete inFile_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 4){
    std::cout << "Usage: ./bin/checkOverlap.exe <fileName1> <subsetTrig1> <fileName2>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += checkOverlap(argv[1], argv[2], argv[3]);
  return retVal;
}
