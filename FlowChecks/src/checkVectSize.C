#include <iostream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"

int checkVectSize(const std::string inFileName)
{
  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  TTree* tree_p = (TTree*)inFile_p->Get("v2V3Tree");

  std::vector<float>* eByEPt_p=NULL;

  tree_p->SetBranchStatus("*", 0);
  tree_p->SetBranchStatus("eByEPt", 1);

  tree_p->SetBranchAddress("eByEPt", &eByEPt_p);

  const Int_t nEntries = tree_p->GetEntries();

  std::cout << "Processing nEntries: " << nEntries << std::endl;
  for(Int_t entry = 0; entry < nEntries; ++entry){
    if(entry%10000 == 0) std::cout << " Entry " << entry << "/" << nEntries << std::endl;
    tree_p->GetEntry(entry);
					       
    if(eByEPt_p->size() == 0){
      std::cout << "Zero entries i event: " << entry << std::endl;
      break;
    }
  }
  
  inFile_p->Close();
  delete inFile_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/checkVectSize.exe <inFileName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += checkVectSize(argv[1]);
  return retVal;
}
