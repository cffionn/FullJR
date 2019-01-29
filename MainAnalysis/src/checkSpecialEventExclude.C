//cpp dependencies
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

//ROOT dependencies
#include "TFile.h"
#include "TTree.h"

//Non-local dependencies
#include "Utility/include/goodGlobalSelection.h"
#include "Utility/include/mntToXRootdFileString.h"
#include "Utility/include/specialHYDJETEventExclude.h"

int checkSpecialEventExclude(const std::string inName)
{
  std::vector<std::string> fileList;

  if(inName.find(".root") != std::string::npos) fileList.push_back(inName);
  else if(inName.find(".txt") != std::string::npos){
    std::ifstream file(inName.c_str());
    std::string tempStr;

    while(std::getline(file, tempStr)){
      while(tempStr.find(" ") != std::string::npos){tempStr.replace(tempStr.find(" "), 1, "");}
      if(tempStr.size() == 0) continue;
      if(tempStr.find(".root") != std::string::npos) fileList.push_back(tempStr);
    }

    file.close();
  }
  else{
    std::cout << "Given inName \'" << inName << "\' is invalid. return 1" << std::endl;
    return 1;
  }

  if(fileList.size() == 0){
    std::cout << "Given inName \'" << inName << "\' is gives no valid root files. return 1" << std::endl;
    return 1;
  }

  const Int_t nMaxJet_ = 500;
  Int_t ngen_;
  Float_t genpt_[nMaxJet_];
  Float_t geneta_[nMaxJet_];
  Float_t genphi_[nMaxJet_];
  Int_t gensubid_[nMaxJet_];

  unsigned int run_, lumi_;
  unsigned long long evt_;
  Int_t hiBin_;
  Float_t vz_;
  Float_t hiHF_;

  Int_t HBHENoiseFilterResultRun2Loose_ = -1;
  Int_t pprimaryVertexFilter_ = -1;
  Int_t pBeamScrapingFilter_ = -1;
  Int_t phfCoincFilter3_ = -1;
  Int_t pclusterCompatibilityFilter_ = -1;

  goodGlobalSelection globalSel;
  globalSel.setIsPbPb(true);

  specialHYDJETEventExclude specialSel;

  for(unsigned int fI = 0; fI < fileList.size(); ++fI){
    std::cout << "Processing file " << fI << "/" << fileList.size() << ": \'" << fileList.at(fI) << "\'" << std::endl;

    TFile* inFile_p = TFile::Open(mntToXRootdFileString(fileList.at(fI)).c_str(), "READ");
    TTree* jetTree_p = (TTree*)inFile_p->Get("akCs4PU3PFFlowJetAnalyzer/t");
    jetTree_p->SetBranchStatus("*", 0);
    jetTree_p->SetBranchStatus("ngen", 1);
    jetTree_p->SetBranchStatus("genpt", 1);
    jetTree_p->SetBranchStatus("geneta", 1);
    jetTree_p->SetBranchStatus("genphi", 1);
    jetTree_p->SetBranchStatus("gensubid", 1);

    jetTree_p->SetBranchAddress("ngen", &ngen_);
    jetTree_p->SetBranchAddress("genpt", genpt_);
    jetTree_p->SetBranchAddress("geneta", geneta_);
    jetTree_p->SetBranchAddress("genphi", genphi_);
    jetTree_p->SetBranchAddress("gensubid", gensubid_);

    TTree* hiTree_p = (TTree*)inFile_p->Get("hiEvtAnalyzer/HiTree");

    hiTree_p->SetBranchStatus("*", 0);
    hiTree_p->SetBranchStatus("run", 1);
    hiTree_p->SetBranchStatus("lumi", 1);
    hiTree_p->SetBranchStatus("evt", 1);
    hiTree_p->SetBranchStatus("hiBin", 1);
    hiTree_p->SetBranchStatus("hiHF", 1);
    hiTree_p->SetBranchStatus("vz", 1);

    hiTree_p->SetBranchAddress("run", &run_);
    hiTree_p->SetBranchAddress("lumi", &lumi_);
    hiTree_p->SetBranchAddress("evt", &evt_);
    hiTree_p->SetBranchAddress("hiBin", &hiBin_);
    hiTree_p->SetBranchAddress("hiHF", &hiHF_);
    hiTree_p->SetBranchAddress("vz", &vz_);

    TTree* skimTree_p = (TTree*)inFile_p->Get("skimanalysis/HltTree");

    skimTree_p->SetBranchStatus("*", 0);
    skimTree_p->SetBranchStatus("HBHENoiseFilterResultRun2Loose", 1);
    skimTree_p->SetBranchAddress("HBHENoiseFilterResultRun2Loose", &HBHENoiseFilterResultRun2Loose_);

    if(true){//use this with pbpb
      skimTree_p->SetBranchStatus("pprimaryVertexFilter", 1);
      skimTree_p->SetBranchStatus("phfCoincFilter3", 1);
      skimTree_p->SetBranchStatus("pclusterCompatibilityFilter", 1);

      skimTree_p->SetBranchAddress("pprimaryVertexFilter", &pprimaryVertexFilter_);
      skimTree_p->SetBranchAddress("phfCoincFilter3", &phfCoincFilter3_);
      skimTree_p->SetBranchAddress("pclusterCompatibilityFilter", &pclusterCompatibilityFilter_);
    }
    else{
      skimTree_p->SetBranchStatus("pBeamScrapingFilter", 1);
      skimTree_p->SetBranchStatus("pPAprimaryVertexFilter", 1);

      skimTree_p->SetBranchAddress("pBeamScrapingFilter", &pBeamScrapingFilter_);
      skimTree_p->SetBranchAddress("pPAprimaryVertexFilter", &pprimaryVertexFilter_);
    }


    const Int_t nEntries = jetTree_p->GetEntries();
    const Int_t printInterval = TMath::Max(1, nEntries/20);

    Int_t totalEvents = 0;
    Int_t eventsPassingGoodGlobalSel = 0;
    Int_t eventsPassingSpecialSel = 0;

    for(Int_t entry = 0; entry < nEntries; ++entry){
      if(nEntries >= 50000 && entry%printInterval == 0) std::cout << " Entry: " << entry << "/" << nEntries << std::endl;

      jetTree_p->GetEntry(entry);
      hiTree_p->GetEntry(entry);
      skimTree_p->GetEntry(entry);

      totalEvents++;

      globalSel.setVz(vz_);
      globalSel.setHiHF(hiHF_);
      globalSel.setPprimaryVertexFilter(pprimaryVertexFilter_);
      globalSel.setPBeamScrapingFilter(pBeamScrapingFilter_);
      globalSel.setPhfCoincFilter3(phfCoincFilter3_);
      globalSel.setHBHENoiseFilterResultRun2Loose(HBHENoiseFilterResultRun2Loose_);
      globalSel.setPclusterCompatibilityFilter(pclusterCompatibilityFilter_);

      if(!globalSel.isGood()) continue;

      eventsPassingGoodGlobalSel++;

      bool badJetSpecialSel = specialSel.CheckEventBadJet(ngen_, genpt_, genphi_, geneta_, gensubid_);

      if(badJetSpecialSel) continue;
      
      eventsPassingSpecialSel++;
    }

    inFile_p->Close();
    delete inFile_p;

    std::cout << " Total Evt: " << totalEvents << std::endl;
    std::cout << "  Passing goodGlobal: " << eventsPassingGoodGlobalSel << ", " << 100.*((Double_t)eventsPassingGoodGlobalSel)/((Double_t)totalEvents) << "%" << std::endl;
    std::cout << "  Passing specialSel: " << eventsPassingSpecialSel<< ", " << 100.*((Double_t)eventsPassingSpecialSel)/((Double_t)totalEvents) << "%" << std::endl;
  }

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/checkSpecialEventExclude.exe <inFileName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += checkSpecialEventExclude(argv[1]);

  std::cout << "Job complete. Return " << retVal << "." << std::endl;
  return retVal;
}
