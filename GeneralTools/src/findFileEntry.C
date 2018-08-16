#include <iostream>
#include <string>
#include <vector>
#include <fstream>

#include "TFile.h"
#include "TTree.h"

int findFileEntry(const std::string inName, const unsigned int inRun, const unsigned int inLumi, const unsigned long long inEvt)
{
  std::vector<std::string> fileList;
  if(inName.find(".root") != std::string::npos) fileList.push_back(inName);
  else if(inName.find(".txt") != std::string::npos){
    std::ifstream file(inName.c_str());
    std::string tempStr;
    while(std::getline(file, tempStr)){
      if(tempStr.size() == 0) continue;

      if(tempStr.find(".root") != std::string::npos) fileList.push_back(tempStr);
    }
    file.close();
  }
  else{
    std::cout << "inName \'" << inName << "\' is invalid. return 1" << std::endl;
    return 1;
  }

  if(fileList.size() == 0){
    std::cout << "inName \'" << inName << "\' gives no valid files. return 1" << std::endl;
    return 1;
  }

  std::cout << "Begin processing " << fileList.size() << " files..." << std::endl;
  
  bool isFound = false;
  std::string foundFile = "";
  int foundEntry = -1;

  unsigned int run_, lumi_;
  unsigned long long evt_;

  for(unsigned int fI = 0; fI < fileList.size(); ++fI){
    std::cout << " File: " << fI << "/" << fileList.size() << ", " << fileList.at(fI) << std::endl;

    TFile* inFile_p = new TFile(fileList.at(fI).c_str(), "READ");
    TTree* hiTree_p = (TTree*)inFile_p->Get("hiEvtAnalyzer/HiTree");
    hiTree_p->SetBranchStatus("*", 0);
    hiTree_p->SetBranchStatus("run", 1);
    hiTree_p->SetBranchStatus("lumi", 1);
    hiTree_p->SetBranchStatus("evt", 1);

    hiTree_p->SetBranchAddress("run", &run_);
    hiTree_p->SetBranchAddress("lumi", &lumi_);
    hiTree_p->SetBranchAddress("evt", &evt_);

    const Int_t nEntries = hiTree_p->GetEntries();

    for(Int_t entry = 0; entry < nEntries; ++entry){
      hiTree_p->GetEntry(entry);

      if(run_ != inRun) continue;
      if(lumi_ != inLumi) continue;
      if(evt_ != inEvt) continue;

      foundFile = fileList.at(fI);
      foundEntry = entry;
      
      isFound = true;
      break;
    }

    if(isFound) break;

    inFile_p->Close();
    delete inFile_p;
  }

  std::cout << "Run,lumi,evt: " << inRun << ", " << inLumi << ", " << inEvt << std::endl;

  if(isFound) std::cout << "Is found in file \'" << foundFile << "\', entry " << foundEntry << std::endl;
  else std::cout << "Is not found..." << std::endl;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 5){
    std::cout << "Usage: ./bin/findFileEntry.exe <inName> <run> <lumi> <evt>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += findFileEntry(argv[1], std::stoi(argv[2]), std::stoi(argv[3]), std::stol(argv[4]));
  return retVal;
}
