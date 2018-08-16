//cpp dependencies
#include <iostream>
#include <string>

//ROOT dependencies
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TDatime.h"
#include "TNamed.h"

#include "Utility/include/checkMakeDir.h"
#include "Utility/include/etaPhiFunc.h"

int makeJetRRatio(const std::string inFileName)
{
  if(!checkFile(inFileName)){
    std::cout << "makePartonRRatio: inFileName \'" << inFileName << "\' is not a valid file. return 1" << std::endl;
    return 1;
  }
  else if(inFileName.find(".root") == std::string::npos){
    std::cout << "makePartonRRatio: inFileName \'" << inFileName << "\' is not a .root file. return 1" << std::endl;
    return 1;
  }

 TDatime* date = new TDatime();
  std::string outFileName = inFileName;
  while(outFileName.find("/") != std::string::npos){outFileName.replace(0, outFileName.find("/")+1, "");}
  outFileName.replace(outFileName.find(".root"), 5, "");
  outFileName = "output/" + outFileName + "_JetRRatio_" + std::to_string(date->GetDate()) + ".root";
  delete date;

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  TH2D* r0p4Over0p8vsR0p8Pt_h = new TH2D("r0p4Over0p8vsR0p8Pt_h", ";Gen. Jet p_{T} R=0.8;R=0.4/R=0.8", 100, 30, 330, 100, 0, 2);


  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  TTree* inTreeR0p4_p = (TTree*)inFile_p->Get("recoAndGenTreeR0p4");
  TTree* inTreeR0p8_p = (TTree*)inFile_p->Get("recoAndGenTreeR0p8");

  Float_t pthatWeight_;

  const Int_t nPartMax = 2000;
  Int_t nJtR0p4_;
  Float_t genJtPtR0p4_[nPartMax];
  Float_t genJtPhiR0p4_[nPartMax];
  Float_t genJtEtaR0p4_[nPartMax];

  Int_t nJtR0p8_;
  Float_t genJtPtR0p8_[nPartMax];
  Float_t genJtPhiR0p8_[nPartMax];
  Float_t genJtEtaR0p8_[nPartMax];

  inTreeR0p4_p->SetBranchStatus("*", 0);
  inTreeR0p4_p->SetBranchStatus("pthatWeight", 1);
  inTreeR0p4_p->SetBranchStatus("nJt", 1);
  inTreeR0p4_p->SetBranchStatus("genJtPt", 1);
  inTreeR0p4_p->SetBranchStatus("genJtPhi", 1);
  inTreeR0p4_p->SetBranchStatus("genJtEta", 1);

  inTreeR0p4_p->SetBranchAddress("pthatWeight", &pthatWeight_);
  inTreeR0p4_p->SetBranchAddress("nJt", &nJtR0p4_);
  inTreeR0p4_p->SetBranchAddress("genJtPt", genJtPtR0p4_);
  inTreeR0p4_p->SetBranchAddress("genJtPhi", genJtPhiR0p4_);
  inTreeR0p4_p->SetBranchAddress("genJtEta", genJtEtaR0p4_);

  inTreeR0p8_p->SetBranchStatus("*", 0);
  inTreeR0p8_p->SetBranchStatus("nJt", 1);
  inTreeR0p8_p->SetBranchStatus("genJtPt", 1);
  inTreeR0p8_p->SetBranchStatus("genJtPhi", 1);
  inTreeR0p8_p->SetBranchStatus("genJtEta", 1);

  inTreeR0p8_p->SetBranchAddress("nJt", &nJtR0p8_);
  inTreeR0p8_p->SetBranchAddress("genJtPt", genJtPtR0p8_);
  inTreeR0p8_p->SetBranchAddress("genJtPhi", genJtPhiR0p8_);
  inTreeR0p8_p->SetBranchAddress("genJtEta", genJtEtaR0p8_);

  const Int_t nEntries = inTreeR0p4_p->GetEntries();

  std::cout << "Processing " << nEntries << " entries..." << std::endl;
  for(Int_t entry = 0; entry < nEntries; ++entry){
    if(entry%100000 == 0) std::cout << " Entry " << entry << "/" << nEntries << std::endl;
    inTreeR0p4_p->GetEntry(entry);
    inTreeR0p8_p->GetEntry(entry);

    const Int_t nUsed = nJtR0p4_;
    bool isUsed[nUsed];
    for(Int_t rI = 0; rI < nJtR0p4_; ++rI){
      isUsed[rI] = false;
    }

    for(Int_t jI = 0; jI < nJtR0p8_; ++jI){
      for(Int_t jI2 = 0; jI2 < nJtR0p4_; ++jI2){
	if(isUsed[jI2]) continue;

	if(getDR(genJtEtaR0p4_[jI2], genJtPhiR0p4_[jI2], genJtEtaR0p8_[jI], genJtPhiR0p8_[jI]) < 0.8){
	  isUsed[jI2] = true;

	  if(genJtPtR0p4_[jI2]/genJtPtR0p8_[jI] > 1.3){
	    std::cout << "Greater than 1 (R=0.4,R=0.8): " << genJtPtR0p4_[jI2]/genJtPtR0p8_[jI] << " (" << genJtPtR0p4_[jI2] << ", " << genJtPtR0p8_[jI] << "), Entry=" << entry << std::endl;
	  }

	  r0p4Over0p8vsR0p8Pt_h->Fill(genJtPtR0p8_[jI], genJtPtR0p4_[jI2]/genJtPtR0p8_[jI], pthatWeight_);
	}
      }
    }

  }

  inFile_p->Close();
  delete inFile_p;

  outFile_p->cd();
  
  r0p4Over0p8vsR0p8Pt_h->Write("", TObject::kOverwrite);
  delete r0p4Over0p8vsR0p8Pt_h;

  outFile_p->Close();
  delete outFile_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/makeJetRRatio.exe <inFileName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += makeJetRRatio(argv[1]);
  return retVal;
}
