#include <iostream>
#include <string>
#include <fstream>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TDatime.h"

#include "Utility/include/checkMakeDir.h"
#include "Utility/include/goodGlobalSelection.h"
#include "Utility/include/mntToXRootdFileString.h"
#include "Utility/include/ncollFunctions_5TeV.h"

int makeMuonFakeCheck(const std::string inName, bool isPP = false)
{
  std::vector<std::string> fileList;
  std::vector<double> pthats;
  std::vector<double> pthatWeights;

  if(inName.find(".root") != std::string::npos){
    fileList.push_back(inName);
    pthats.push_back(1.);
    pthatWeights.push_back(1.);
  }
  else if(inName.find(".txt") != std::string::npos){
    std::ifstream file(inName.c_str());
    std::string tempStr;

    while(std::getline(file, tempStr)){
      while(tempStr.find(" ") != std::string::npos){tempStr.replace(tempStr.find(" "), 1, "");}
      if(tempStr.size() == 0) continue;
      if(tempStr.find(".root") != std::string::npos) fileList.push_back(tempStr);
      else{
        if(tempStr.substr(0, std::string("PTHAT=").size()).find("PTHAT=") != std::string::npos){
          tempStr.replace(0, std::string("PTHAT=").size(), "");
          while(tempStr.find(",") != std::string::npos){
            pthats.push_back(std::stod(tempStr.substr(0, tempStr.find(","))));
            tempStr.replace(0, tempStr.find(",")+1, "");
          }
          if(tempStr.size() != 0) pthats.push_back(std::stod(tempStr));
        }
        else if(tempStr.substr(0, std::string("PTHATWEIGHTS=").size()).find("PTHATWEIGHTS=") != std::string::npos){
          tempStr.replace(0, std::string("PTHATWEIGHTS=").size(), "");
          while(tempStr.find(",") != std::string::npos){
            pthatWeights.push_back(std::stod(tempStr.substr(0, tempStr.find(","))));
            tempStr.replace(0, tempStr.find(",")+1, "");
          }
          if(tempStr.size() != 0) pthatWeights.push_back(std::stod(tempStr));
        }
        else if(tempStr.substr(0, std::string("ISPP=").size()).find("ISPP=") != std::string::npos){
          tempStr.replace(0,tempStr.find("=")+1, "");
          isPP = std::stoi(tempStr);
        }
        else std::cout << "WARNING: Line in \'" << inName << "\', \'" << tempStr << "\' is invalid. check input" << std::endl;
      }
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
  else if(pthats.size() == 0){
    std::cout << "Given inName \'" << inName << "\' contains no pthat list. return 1" << std::endl;
    return 1;
  }
  else if(pthatWeights.size() == 0){
    std::cout << "Given inName \'" << inName << "\' contains no pthatWeights list. return 1" << std::endl;
    return 1;
  }

  std::cout << "Pthats and weights: " << std::endl;
  for(unsigned int pI = 0; pI < pthats.size(); ++pI){
    std::cout << " " << pI << "/" << pthats.size() << ": " << pthats.at(pI) << ", " << pthatWeights.at(pI) << std::endl;
  }

  std::cout << "IsPP: " << isPP << std::endl;

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  std::string outFileName = inName;
  while(outFileName.find("/") != std::string::npos){outFileName.replace(0, outFileName.find("/")+1, "");}
  if(outFileName.find(".root") != std::string::npos) outFileName.replace(outFileName.find(".root"), std::string(".root").size(), "");
  else if(outFileName.find(".txt") != std::string::npos) outFileName.replace(outFileName.find(".txt"), std::string(".txt").size(), "");
  outFileName = "output/" + outFileName + "_MuonFakeCheck_" + dateStr + ".root";
  checkMakeDir("output");

  const Double_t jtAbsEtaMax = 2.;
  const Int_t nJtPtBins = 8;
  const Double_t jtPtBins[nJtPtBins+1] = {300., 400., 500., 600., 700., 800., 900., 1000., 1100.};

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  TH1D* recoJtPt_Raw_h = new TH1D("recoJtPt_Raw_h", ";Reco. Jet p_{T};Counts (Weighted)", nJtPtBins, jtPtBins);
  TH1D* recoJtPt_RawWithRef_h = new TH1D("recoJtPt_RawWithRef_h", ";Reco. Jet p_{T};Counts (Weighted)", nJtPtBins, jtPtBins);

  TH1D* recoJtPt_Raw_Loose_h = new TH1D("recoJtPt_Raw_Loose_h", ";Reco. Jet p_{T};Counts (Weighted)", nJtPtBins, jtPtBins);
  TH1D* recoJtPt_RawWithRef_Loose_h = new TH1D("recoJtPt_RawWithRef_Loose_h", ";Reco. Jet p_{T};Counts (Weighted)", nJtPtBins, jtPtBins);

  TH1D* recoJtPt_Raw_Tight_h = new TH1D("recoJtPt_Raw_Tight_h", ";Reco. Jet p_{T};Counts (Weighted)", nJtPtBins, jtPtBins);
  TH1D* recoJtPt_RawWithRef_Tight_h = new TH1D("recoJtPt_RawWithRef_Tight_h", ";Reco. Jet p_{T};Counts (Weighted)", nJtPtBins, jtPtBins);

  TH1D* recoJtPt_Raw_TightLepVeto_h = new TH1D("recoJtPt_Raw_TightLepVeto_h", ";Reco. Jet p_{T};Counts (Weighted)", nJtPtBins, jtPtBins);
  TH1D* recoJtPt_RawWithRef_TightLepVeto_h = new TH1D("recoJtPt_RawWithRef_TightLepVeto_h", ";Reco. Jet p_{T};Counts (Weighted)", nJtPtBins, jtPtBins);

  const Int_t nMaxJet_ = 500;
  Float_t pthat_;
  Int_t nref_;
  Float_t jtpt_[nMaxJet_];
  Float_t refpt_[nMaxJet_];
  Float_t jteta_[nMaxJet_];
  Float_t jtphi_[nMaxJet_];
  Float_t jtPfCHF_[nMaxJet_];
  Float_t jtPfCEF_[nMaxJet_];
  Float_t jtPfNHF_[nMaxJet_];
  Float_t jtPfNEF_[nMaxJet_];
  Float_t jtPfMUF_[nMaxJet_];
  Int_t jtPfCHM_[nMaxJet_];
  Int_t jtPfCEM_[nMaxJet_];
  Int_t jtPfNHM_[nMaxJet_];
  Int_t jtPfNEM_[nMaxJet_];
  Int_t jtPfMUM_[nMaxJet_];

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
  globalSel.setIsPbPb(!isPP);

  for(unsigned int fI = 0; fI < fileList.size(); ++fI){
    std::cout << "Processing file " << fI << "/" << fileList.size() << ": \'" << fileList.at(fI) << "\'" << std::endl;

    TFile* inFile_p = TFile::Open(mntToXRootdFileString(fileList.at(fI)).c_str(), "READ");
    TTree* jetTree_p = (TTree*)inFile_p->Get("akCs4PU3PFFlowJetAnalyzer/t");
    jetTree_p->SetBranchStatus("*", 0);
    jetTree_p->SetBranchStatus("pthat", 1);
    jetTree_p->SetBranchStatus("nref", 1);
    jetTree_p->SetBranchStatus("jtpt", 1);
    jetTree_p->SetBranchStatus("refpt", 1);
    jetTree_p->SetBranchStatus("jteta", 1);
    jetTree_p->SetBranchStatus("jtphi", 1);
    jetTree_p->SetBranchStatus("jtPfCHF", 1);
    jetTree_p->SetBranchStatus("jtPfCEF", 1);
    jetTree_p->SetBranchStatus("jtPfNHF", 1);
    jetTree_p->SetBranchStatus("jtPfNEF", 1);
    jetTree_p->SetBranchStatus("jtPfMUF", 1);
    jetTree_p->SetBranchStatus("jtPfCHM", 1);
    jetTree_p->SetBranchStatus("jtPfCEM", 1);
    jetTree_p->SetBranchStatus("jtPfNHM", 1);
    jetTree_p->SetBranchStatus("jtPfNEM", 1);
    jetTree_p->SetBranchStatus("jtPfMUM", 1);

    jetTree_p->SetBranchAddress("pthat", &pthat_);
    jetTree_p->SetBranchAddress("nref", &nref_);
    jetTree_p->SetBranchAddress("jtpt", jtpt_);
    jetTree_p->SetBranchAddress("refpt", refpt_);
    jetTree_p->SetBranchAddress("jteta", jteta_);
    jetTree_p->SetBranchAddress("jtphi", jtphi_);
    jetTree_p->SetBranchAddress("jtPfCHF", jtPfCHF_);
    jetTree_p->SetBranchAddress("jtPfCEF", jtPfCEF_);
    jetTree_p->SetBranchAddress("jtPfNHF", jtPfNHF_);
    jetTree_p->SetBranchAddress("jtPfNEF", jtPfNEF_);
    jetTree_p->SetBranchAddress("jtPfMUF", jtPfMUF_);
    jetTree_p->SetBranchAddress("jtPfCHM", jtPfCHM_);
    jetTree_p->SetBranchAddress("jtPfCEM", jtPfCEM_);
    jetTree_p->SetBranchAddress("jtPfNHM", jtPfNHM_);
    jetTree_p->SetBranchAddress("jtPfNEM", jtPfNEM_);
    jetTree_p->SetBranchAddress("jtPfMUM", jtPfMUM_);

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

    if(!isPP){
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

    for(Int_t entry = 0; entry < nEntries; ++entry){
      if(nEntries >= 50000 && entry%printInterval == 0) std::cout << " Entry: " << entry << "/" << nEntries << std::endl;

      jetTree_p->GetEntry(entry);
      hiTree_p->GetEntry(entry);
      skimTree_p->GetEntry(entry);

      globalSel.setVz(vz_);
      globalSel.setHiHF(hiHF_);
      globalSel.setPprimaryVertexFilter(pprimaryVertexFilter_);
      globalSel.setPBeamScrapingFilter(pBeamScrapingFilter_);
      globalSel.setPhfCoincFilter3(phfCoincFilter3_);
      globalSel.setHBHENoiseFilterResultRun2Loose(HBHENoiseFilterResultRun2Loose_);
      globalSel.setPclusterCompatibilityFilter(pclusterCompatibilityFilter_);

      if(!globalSel.isGood()) continue;
      if(hiBin_ > 20) continue;

      Double_t pthatWeight_ = -1;
      for(unsigned int pI = 0; pI < pthats.size()-1; ++pI){
        if(pthats.at(pI) <= pthat_ && pthat_ < pthats.at(pI+1)){
          pthatWeight_ = pthatWeights.at(pI);
          break;
        }
      }
      if(pthat_ > pthats.at(pthats.size()-1)) pthatWeight_ = pthatWeights.at(pthatWeights.size()-1);

      if(pthatWeight_ < 0){
	std::cout << "WARNING - NO WEIGHT FOR pthat \'" << pthat_ << "\'. Set to 1" << std::endl;
        pthatWeight_ = 1.;
      }

      Double_t ncollWeight_ = findNcoll_Renorm(hiBin_);
      Double_t fullWeight_ = ncollWeight_*pthatWeight_;

      for(Int_t jI = 0; jI < nref_; ++jI){
	if(TMath::Abs(jteta_[jI]) > jtAbsEtaMax) continue;
	if(jtpt_[jI] < 300.) continue;

	recoJtPt_Raw_h->Fill(jtpt_[jI], fullWeight_);
	if(refpt_[jI] >= 150.) recoJtPt_RawWithRef_h->Fill(jtpt_[jI], fullWeight_);
	else if(jtpt_[jI] > 600.){
	  std::cout << " Bad jetpt, phi, eta: " << jtpt_[jI] << ", " << jtphi_[jI] << ", " << jteta_[jI] << std::endl;
	  std::cout << "  run, lumi, evt, entry: " << run_ << ", " << lumi_ << ", " << evt_ << ", " << entry << std::endl;
	  std::cout << " check if cut by loose, tight, lepveto..." << std::endl;
	}

	//via https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2016

	bool passesLoose = jtPfNHF_[jI] < 0.99 && jtPfNEF_[jI] < 0.99;
	passesLoose = passesLoose && jtPfCHM_[jI] + jtPfCEM_[jI] + jtPfNHM_[jI] + jtPfNEM_[jI] + jtPfMUM_[jI] > 1;
	passesLoose = passesLoose && jtPfCHF_[jI] > 0.0;
	passesLoose = passesLoose && jtPfCHM_[jI] > 0;
	passesLoose = passesLoose && jtPfCEF_[jI] < 0.99;

	if(!passesLoose) continue;

        recoJtPt_Raw_Loose_h->Fill(jtpt_[jI], fullWeight_);
        if(refpt_[jI] >= 150.) recoJtPt_RawWithRef_Loose_h->Fill(jtpt_[jI], fullWeight_);
	else if(jtpt_[jI] > 600.){
	  std::cout << "  NOT CUT BY LOOSE" << std::endl;
	}

	bool passesTight = jtPfNHF_[jI] < 0.9 && jtPfNEF_[jI] < 0.9;
	passesTight = passesTight && jtPfCHM_[jI] + jtPfCEM_[jI] + jtPfNHM_[jI] + jtPfNEM_[jI] + jtPfMUM_[jI] > 1;
	passesTight = passesTight && jtPfCHF_[jI] > 0.0;
	passesTight = passesTight && jtPfCHM_[jI] > 0;
	passesTight = passesTight && jtPfCEF_[jI] < 0.99;

	if(!passesTight) continue;

	recoJtPt_Raw_Tight_h->Fill(jtpt_[jI], fullWeight_);
        if(refpt_[jI] >= 150.) recoJtPt_RawWithRef_Tight_h->Fill(jtpt_[jI], fullWeight_);
	else if(jtpt_[jI] > 600.){
	  std::cout << "  NOT CUT BY TIGHT" << std::endl;
	}

	bool passesTightLepVeto = jtPfNHF_[jI] < 0.9 && jtPfNEF_[jI] < 0.9;
	passesTightLepVeto = passesTightLepVeto && jtPfCHM_[jI] + jtPfCEM_[jI] + jtPfNHM_[jI] + jtPfNEM_[jI] + jtPfMUM_[jI] > 1;
	passesTightLepVeto = passesTightLepVeto && jtPfMUF_[jI] < 0.8;
	passesTightLepVeto = passesTightLepVeto && jtPfCHF_[jI] > 0.0;
	passesTightLepVeto = passesTightLepVeto && jtPfCHM_[jI] > 0;
	passesTightLepVeto = passesTightLepVeto && jtPfCEF_[jI] < 0.9;

	if(!passesTightLepVeto) continue;

	recoJtPt_Raw_TightLepVeto_h->Fill(jtpt_[jI], fullWeight_);
	if(refpt_[jI] >= 150.) recoJtPt_RawWithRef_TightLepVeto_h->Fill(jtpt_[jI], fullWeight_);
	else if(jtpt_[jI] > 600.){
	  std::cout << "  NOT CUT BY LEPVETO" << std::endl;
	}
      }
    }

    inFile_p->Close();
    delete inFile_p;
  }

  outFile_p->cd();

  recoJtPt_Raw_h->Write("", TObject::kOverwrite);
  delete recoJtPt_Raw_h;

  recoJtPt_RawWithRef_h->Write("", TObject::kOverwrite);
  delete recoJtPt_RawWithRef_h;

  recoJtPt_Raw_Loose_h->Write("", TObject::kOverwrite);
  delete recoJtPt_Raw_Loose_h;

  recoJtPt_RawWithRef_Loose_h->Write("", TObject::kOverwrite);
  delete recoJtPt_RawWithRef_Loose_h;

  recoJtPt_Raw_Tight_h->Write("", TObject::kOverwrite);
  delete recoJtPt_Raw_Tight_h;

  recoJtPt_RawWithRef_Tight_h->Write("", TObject::kOverwrite);
  delete recoJtPt_RawWithRef_Tight_h;

  recoJtPt_Raw_TightLepVeto_h->Write("", TObject::kOverwrite);
  delete recoJtPt_Raw_TightLepVeto_h;

  recoJtPt_RawWithRef_TightLepVeto_h->Write("", TObject::kOverwrite);
  delete recoJtPt_RawWithRef_TightLepVeto_h;

  outFile_p->Close();
  delete outFile_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2 && argc != 3){
    std::cout << "Usage: ./bin/makeMuonFakeCheck.exe <inFileName> <isPP-Opt>" << std::endl;
    return 1;
  }

  int retVal = 0; 
  if(argc == 2) retVal += makeMuonFakeCheck(argv[1]);
  else if(argc == 3) retVal += makeMuonFakeCheck(argv[1], std::stoi(argv[2]));
  return retVal;
}
