//cpp dependencies
#include <iostream>
#include <map>
#include <string>
#include <vector>

//ROOT dependencies
#include "TFile.h"
#include "TTree.h"

//Non-local dependencies
#include "Utility/include/checkMakeDir.h"
#include "Utility/include/returnRootFileContentsList.h"

ULong64_t runLumiEvtToKey(UInt_t run, UInt_t lumi, ULong64_t evt)
{
  ULong64_t key = ((ULong64_t)lumi) + ((ULong64_t)(run*10000)) + ((ULong64_t)(evt*100000000000));
  return key;
}

int validateSplit(std::string mainFileName, std::string splitFileNames)
{
  if(!checkFile(mainFileName) || mainFileName.find(".root") == std::string::npos){
    std::cout << "mainFileName \'" << mainFileName << "\' is invalid. return 1" << std::endl;
    return 1;
  }

  if(splitFileNames.size() != 0){
    if(splitFileNames.substr(splitFileNames.size()-1, 1).find(",") == std::string::npos) splitFileNames = splitFileNames + ",";
  }

  const std::string splitFileNamesPerma = splitFileNames;
  std::vector<std::string> fileList;
  while(splitFileNames.find(",") != std::string::npos){
    std::string tempName = splitFileNames.substr(0, splitFileNames.find(","));
    std::cout << tempName << std::endl;
    if(checkFile(tempName) && tempName.find(".root") != std::string::npos) fileList.push_back(tempName);
    splitFileNames.replace(0, splitFileNames.find(",")+1, "");
  }

  if(fileList.size() == 0){
    std::cout << "None of given fileNames in splitFileNames \'" << splitFileNamesPerma << "\' are valid. return 1" << std::endl;
    return 1;
  }

  TFile* inFile_p = new TFile(mainFileName.c_str(), "READ");
  std::vector<std::string> treeListMain = returnRootFileContentsList(inFile_p, "TTree"); 
  std::vector<ULong64_t> countsMain;
  std::vector<ULong64_t> countsSplit;

  std::map<ULong64_t, UInt_t> runLumiEvtMapMain;
  
  UInt_t run, lumi;
  ULong64_t evt;

  for(auto const & tree : treeListMain){
    TTree* inTree_p = (TTree*)inFile_p->Get(tree.c_str());
    countsMain.push_back(inTree_p->GetEntries());
    countsSplit.push_back(0);

    if(tree.find("HiTree") != std::string::npos){
      inTree_p->SetBranchStatus("*", 0);
      inTree_p->SetBranchStatus("run", 1);
      inTree_p->SetBranchStatus("lumi", 1);
      inTree_p->SetBranchStatus("evt", 1);

      inTree_p->SetBranchAddress("run", &run);
      inTree_p->SetBranchAddress("lumi", &lumi);
      inTree_p->SetBranchAddress("evt", &evt);

      const ULong64_t nEntries = inTree_p->GetEntries();

      for(ULong64_t entry = 0; entry < nEntries; ++entry){
	inTree_p->GetEntry(entry);
	ULong64_t key = runLumiEvtToKey(run, lumi, evt);
	if(runLumiEvtMapMain.count(key) != 0){
	  std::cout << "WARNING: key is reused: " << run << ", " << lumi << ", " << evt << std::endl;
	}
	runLumiEvtMapMain[key] = 1;
      }
    }
  }  
  
  inFile_p->Close();
  delete inFile_p;
  
  for(auto const & file : fileList){
    inFile_p = new TFile(file.c_str(), "READ");
    std::vector<std::string> treeListSplit = returnRootFileContentsList(inFile_p, "TTree");
    if(treeListSplit.size() != treeListMain.size()){
      std::cout << "Warning: TTree list mismatch. return 1" << std::endl;
      return 1;
    }

    unsigned int pos = 0;
    for(auto const & tree1 : treeListMain){
      bool treeFound = false;

      for(auto const & tree2 : treeListSplit){
	if(isStrSame(tree1, tree2)){
	  treeFound = true;
	  break;
	}
      }

      if(!treeFound){
	std::cout << "Missing tree: " << tree1 << std::endl;
	return 1;
      }

      TTree* inTree_p = (TTree*)inFile_p->Get(tree1.c_str());
      countsSplit[pos] += inTree_p->GetEntries();
      
      if(tree1.find("HiTree") != std::string::npos){
	inTree_p->SetBranchStatus("*", 0);
	inTree_p->SetBranchStatus("run", 1);
	inTree_p->SetBranchStatus("lumi", 1);
	inTree_p->SetBranchStatus("evt", 1);
	
	inTree_p->SetBranchAddress("run", &run);
	inTree_p->SetBranchAddress("lumi", &lumi);
	inTree_p->SetBranchAddress("evt", &evt);

	const ULong64_t nEntries = inTree_p->GetEntries();
	for(ULong64_t entry = 0; entry < nEntries; ++entry){
	  inTree_p->GetEntry(entry);
	  ULong64_t key = runLumiEvtToKey(run, lumi, evt);

	  if(runLumiEvtMapMain.count(key) != 1){
	    std::cout << "ENTRY MISSING: " << run << ", " << lumi << ", " << evt << std::endl;
	  }
	}
      }

      ++pos;
    }
    
    inFile_p->Close();
    delete inFile_p;
  }

  for(unsigned int tI = 0; tI < treeListMain.size(); ++tI){
    std::cout << treeListMain[tI] << ": " << countsMain[tI] << "/" << countsSplit[tI] << std::endl;
  }
  
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 3){
    std::cout << "Usage: ./bin/validateSplit.exe <mainFileName> <splitFileName>. return 1" << std::endl;
    return 1;
  }
  
  int retVal = 0;
  retVal += validateSplit(argv[1], argv[2]);
  return retVal;
}
