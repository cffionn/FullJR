//cpp dependencies
#include <iostream>
#include <string>

//ROOT dependencies
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TDatime.h"
#include "TMath.h"
#include "TNamed.h"

//Non-Local FullJR (Utility, etc.) dependencies
#include "Utility/include/checkMakeDir.h"
#include "Utility/include/etaPhiFunc.h"
#include "Utility/include/returnRootFileContentsList.h"
#include "Utility/include/histDefUtility.h"
#include "Utility/include/plotUtilities.h"
#include "Utility/include/getLinBins.h"

int makePartonRRatio(const std::string inFileName)
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
  outFileName = "output/" + outFileName + "_PartonRRatio_" + std::to_string(date->GetDate()) + ".root";
  delete date;

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");

  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  TTree* inTree_p = (TTree*)inFile_p->Get("partonTree");
  TTree* jetTreeR0p8_p = (TTree*)inFile_p->Get("recoAndGenTreeR0p8");

  std::vector<std::string> listOfTNamed = returnRootFileContentsList(inFile_p, "TNamed", "");

  unsigned int pos = 0;
  while(pos < listOfTNamed.size()){   
    bool isUnique = true;
    for(unsigned int i = pos+1; i < listOfTNamed.size(); ++i){
      if(listOfTNamed.at(i).size() == listOfTNamed.at(pos).size() && listOfTNamed.at(i).find(listOfTNamed.at(pos)) != std::string::npos){
	isUnique = false;
	listOfTNamed.erase(listOfTNamed.begin()+i);
      }
    }

    if(isUnique) ++pos;
  }

  Int_t nEnergyFracBinsPos = -1;
  Int_t energyFracLowPos = -1; 
  Int_t energyFracHiPos = -1;
  Int_t energyFracBinsPos = -1;
 
  std::cout << "Check tnamed: " << std::endl;
  for(unsigned int i = 0; i < listOfTNamed.size(); ++i){
    std::cout << " " << i << "/" << listOfTNamed.size() << ": " << listOfTNamed.at(i) << std::endl;

    std::string tempStr = listOfTNamed.at(i);
    while(tempStr.find("/") != std::string::npos){tempStr.replace(0, tempStr.find("/")+1, "");}

    if(tempStr.find("nEnergyFracBins") != std::string::npos && tempStr.size() == std::string("nEnergyFracBins").size()) nEnergyFracBinsPos = i;

    if(tempStr.find("energyFracLow") != std::string::npos && tempStr.size() == std::string("energyFracLow").size()) energyFracLowPos = i;

    if(tempStr.find("energyFracHi") != std::string::npos && tempStr.size() == std::string("energyFracHi").size()) energyFracHiPos = i;

    if(tempStr.find("energyFracBins") != std::string::npos && tempStr.size() == std::string("energyFracBins").size()) energyFracBinsPos = i;
 }

  std::cout << "Gets here: " << __LINE__ << std::endl;
  std::cout << nEnergyFracBinsPos << std::endl;
  std::cout << energyFracLowPos << std::endl;
  std::cout << energyFracHiPos << std::endl;
  std::cout << energyFracBinsPos << std::endl;

  const Int_t nEnergyFracBins = std::stoi(std::string(((TNamed*)inFile_p->Get(listOfTNamed.at(nEnergyFracBinsPos).c_str()))->GetTitle()));
  const Float_t energyFracLow = std::stof(std::string(((TNamed*)inFile_p->Get(listOfTNamed.at(energyFracLowPos).c_str()))->GetTitle()));
  const Float_t energyFracHi = std::stof(std::string(((TNamed*)inFile_p->Get(listOfTNamed.at(energyFracHiPos).c_str()))->GetTitle()));
  Double_t energyFracBins[nEnergyFracBins+1];
  
  std::string energyFracBinsStr = ((TNamed*)inFile_p->Get(listOfTNamed.at(energyFracBinsPos).c_str()))->GetTitle();

  std::cout << "Gets here: " << __LINE__ << std::endl;

  pos = 0;
  while(energyFracBinsStr.find(",") != std::string::npos){
    energyFracBins[pos] = std::stod(energyFracBinsStr.substr(0, energyFracBinsStr.find(",")));
    energyFracBinsStr.replace(0, energyFracBinsStr.find(",")+1, "");
    ++pos;
  }
  if(energyFracBinsStr.size() != 0) energyFracBins[nEnergyFracBins] = std::stod(energyFracBinsStr);

  std::cout << "nEnergyFracBins, energyFracLow, energyFracHi: " << nEnergyFracBins << ", " << energyFracLow << ", " << energyFracHi << std::endl;
  for(Int_t i = 0; i < nEnergyFracBins+1; ++i){
    std::cout << energyFracBins[i] << ",";
  }
  std::cout << std::endl;

  const Int_t nPartonPtBins = 4;
  const Float_t partonPtLowBins[nPartonPtBins] = {30, 50, 100, 200};
  const Float_t partonPtHiBins[nPartonPtBins] = {50, 100, 200, 400};

  const Int_t nJtPtBins = 4;
  const Float_t jtPtLowBins[nJtPtBins] = {30, 50, 100, 200};
  const Float_t jtPtHiBins[nJtPtBins] = {50, 100, 200, 400};

  const Int_t nPartonFracBins = 300;
  const Float_t partonFracLow = 0.;
  const Float_t partonFracHi = 2.;
  Double_t partonFracBins[nPartonFracBins+1];
  getLinBins(partonFracLow, partonFracHi, nPartonFracBins, partonFracBins);

  const Int_t nJtFracBins = 300;
  const Float_t jtFracLow = 0.;
  const Float_t jtFracHi = 2.;
  Double_t jtFracBins[nJtFracBins+1];
  getLinBins(jtFracLow, jtFracHi, nJtFracBins, jtFracBins);

  outFile_p->cd();
  TH2D* partonRRatio_h[nPartonPtBins];
  TH1D* partonRRatio_Bins_h[nPartonPtBins][nEnergyFracBins];
  TH1D* partonRRatio_Mean_h[nPartonPtBins];

  TH2D* jtRRatio_h[nJtPtBins];
  TH1D* jtRRatio_Bins_h[nJtPtBins][nEnergyFracBins];
  TH1D* jtRRatio_Mean_h[nJtPtBins];

  for(Int_t pI = 0; pI < nPartonPtBins; ++pI){
    const std::string partStr = "PartonPt" + prettyString(partonPtLowBins[pI],1,true) + "to" + prettyString(partonPtHiBins[pI],1,true);
    const std::string jtStr = "JtPt" + prettyString(jtPtLowBins[pI],1,true) + "to" + prettyString(jtPtHiBins[pI],1,true);
    
    partonRRatio_h[pI] = new TH2D(("partonRRatio_" + partStr + "_h").c_str(), ";#DeltaR_{Hadrons,Parton};Energy Fraction", nEnergyFracBins, energyFracBins, nPartonFracBins, partonFracBins);
    jtRRatio_h[pI] = new TH2D(("jtRRatio_" + partStr + "_h").c_str(), ";#DeltaR_{Hadrons,Jt};Energy Fraction", nEnergyFracBins, energyFracBins, nJtFracBins, jtFracBins);
    
    for(Int_t fI = 0; fI < nEnergyFracBins; ++fI){
      partonRRatio_Bins_h[pI][fI] = new TH1D(("partonRRatio_Bins_" + partStr + "_Bin" + std::to_string(fI) + "_h").c_str(), "", 100, partonFracLow, partonFracHi);
      jtRRatio_Bins_h[pI][fI] = new TH1D(("jtRRatio_Bins_" + jtStr + "_Bin" + std::to_string(fI) + "_h").c_str(), "", 100, jtFracLow, jtFracHi);
      centerTitles({partonRRatio_Bins_h[pI][fI], jtRRatio_Bins_h[pI][fI]});
      setSumW2({partonRRatio_Bins_h[pI][fI], jtRRatio_Bins_h[pI][fI]});
    }

    partonRRatio_Mean_h[pI] = new TH1D(("partonRRatio_Mean_" + partStr + "_h").c_str(), ";#DeltaR_{Hadrons,Parton};#LTEnergy Fraction#GT", nEnergyFracBins, energyFracBins);
    jtRRatio_Mean_h[pI] = new TH1D(("jtRRatio_Mean_" + jtStr + "_h").c_str(), ";#DeltaR_{Hadrons,Jt};#LTEnergy Fraction#GT", nEnergyFracBins, energyFracBins);

    centerTitles({partonRRatio_h[pI], jtRRatio_h[pI]});
    setSumW2({partonRRatio_h[pI], jtRRatio_h[pI]});

    centerTitles({partonRRatio_Mean_h[pI], jtRRatio_Mean_h[pI]});
    setSumW2({partonRRatio_Mean_h[pI], jtRRatio_Mean_h[pI]});

    partonRRatio_Mean_h[pI]->SetMarkerColor(1);
    partonRRatio_Mean_h[pI]->SetMarkerSize(1);
    partonRRatio_Mean_h[pI]->SetMarkerStyle(20);
    partonRRatio_Mean_h[pI]->SetLineColor(1);

    jtRRatio_Mean_h[pI]->SetMarkerColor(1);
    jtRRatio_Mean_h[pI]->SetMarkerSize(1);
    jtRRatio_Mean_h[pI]->SetMarkerStyle(20);
    jtRRatio_Mean_h[pI]->SetLineColor(1);
  }

  Float_t pthatWeight_;

  const Int_t nQGMax_ = 2;
  Int_t nQG_;
  Float_t qgPt_[nQGMax_];
  Float_t qgPhi_[nQGMax_];
  Float_t qgEta_[nQGMax_];
  Int_t qgID_[nQGMax_];
  Float_t energyFrac_[nQGMax_][nEnergyFracBins];

  inTree_p->SetBranchStatus("*", 0);
  inTree_p->SetBranchStatus("pthatWeight", 1);
  inTree_p->SetBranchStatus("nQG", 1);
  inTree_p->SetBranchStatus("qgPt", 1);
  inTree_p->SetBranchStatus("qgPhi", 1);
  inTree_p->SetBranchStatus("qgEta", 1);
  inTree_p->SetBranchStatus("qgID", 1);
  inTree_p->SetBranchStatus("energyFrac", 1);

  inTree_p->SetBranchAddress("pthatWeight", &pthatWeight_);
  inTree_p->SetBranchAddress("nQG", &nQG_);
  inTree_p->SetBranchAddress("qgPt", qgPt_);
  inTree_p->SetBranchAddress("qgPhi", qgPhi_);
  inTree_p->SetBranchAddress("qgEta", qgEta_);
  inTree_p->SetBranchAddress("qgID", qgID_);
  inTree_p->SetBranchAddress("energyFrac", energyFrac_);

  const Int_t nMaxJet = 500;
  Int_t nJt_;
  Float_t genJtPt_[nMaxJet];
  Float_t genJtEta_[nMaxJet];
  Float_t jtEnergyFrac_[nMaxJet][nEnergyFracBins];

  jetTreeR0p8_p->SetBranchStatus("*", 0);
  jetTreeR0p8_p->SetBranchStatus("nJt", 1);
  jetTreeR0p8_p->SetBranchStatus("genJtPt", 1);
  jetTreeR0p8_p->SetBranchStatus("genJtEta", 1);
  jetTreeR0p8_p->SetBranchStatus("jtEnergyFrac", 1);

  jetTreeR0p8_p->SetBranchAddress("nJt", &nJt_);
  jetTreeR0p8_p->SetBranchAddress("genJtPt", genJtPt_);
  jetTreeR0p8_p->SetBranchAddress("genJtEta", genJtEta_);
  jetTreeR0p8_p->SetBranchAddress("jtEnergyFrac", jtEnergyFrac_);

  const Int_t nEntries = inTree_p->GetEntries();

  std::cout << "Processing " << nEntries << " entries..." << std::endl;
  for(Int_t entry = 0; entry < nEntries; ++entry){
    if(entry%100000 == 0) std::cout << " Entry " << entry << "/" << nEntries << std::endl;
    inTree_p->GetEntry(entry);
    jetTreeR0p8_p->GetEntry(entry);

    for(Int_t qI = 0; qI < nQG_; ++qI){
      if(TMath::Abs(qgEta_[qI]) > 1.) continue;

      Int_t partonPos = -1;
      for(Int_t pI = 0; pI < nPartonPtBins; ++pI){
	if(qgPt_[qI] >= partonPtLowBins[pI] && qgPt_[qI] < partonPtHiBins[pI]){
	  partonPos = pI;
	  break;
	}
      }
     
      if(partonPos < 0) continue;

      for(Int_t pI = 0; pI < nEnergyFracBins; ++pI){
	partonRRatio_h[partonPos]->Fill((energyFracBins[pI] + energyFracBins[pI+1])/2., energyFrac_[qI][pI], pthatWeight_);
	partonRRatio_Bins_h[partonPos][pI]->Fill(energyFrac_[qI][pI], pthatWeight_);
      }
    }

    for(Int_t jI = 0; jI < nJt_; ++jI){
      if(TMath::Abs(genJtEta_[jI]) > 1.) continue;

      Int_t jtPos = -1;
      for(Int_t pI = 0; pI < nJtPtBins; ++pI){
        if(genJtPt_[jI] >= jtPtLowBins[pI] && genJtPt_[jI] < jtPtHiBins[pI]){
          jtPos = pI;
          break;
        }
      }

      if(jtPos < 0) continue;

      for(Int_t pI = 0; pI < nEnergyFracBins; ++pI){
        jtRRatio_h[jtPos]->Fill((energyFracBins[pI] + energyFracBins[pI+1])/2., jtEnergyFrac_[jI][pI], pthatWeight_);
        jtRRatio_Bins_h[jtPos][pI]->Fill(jtEnergyFrac_[jI][pI], pthatWeight_);
      }
    }
  }

  inFile_p->Close();
  delete inFile_p;

  outFile_p->cd();

  for(Int_t pI = 0; pI < nPartonPtBins; ++pI){  
    partonRRatio_h[pI]->Write("", TObject::kOverwrite);

    for(Int_t fI = 0; fI < nEnergyFracBins; ++fI){
      //      partonRRatio_Bins_h[pI][fI]->Write("", TObject::kOverwrite);
      partonRRatio_Mean_h[pI]->SetBinContent(fI+1, partonRRatio_Bins_h[pI][fI]->GetMean());
      partonRRatio_Mean_h[pI]->SetBinError(fI+1, partonRRatio_Bins_h[pI][fI]->GetMeanError());
    }

    partonRRatio_Mean_h[pI]->Write("", TObject::kOverwrite);
  }

  for(Int_t pI = 0; pI < nJtPtBins; ++pI){  
    jtRRatio_h[pI]->Write("", TObject::kOverwrite);

    for(Int_t fI = 0; fI < nEnergyFracBins; ++fI){
      //      jtRRatio_Bins_h[pI][fI]->Write("", TObject::kOverwrite);
      jtRRatio_Mean_h[pI]->SetBinContent(fI+1, jtRRatio_Bins_h[pI][fI]->GetMean());
      jtRRatio_Mean_h[pI]->SetBinError(fI+1, jtRRatio_Bins_h[pI][fI]->GetMeanError());
    }

    jtRRatio_Mean_h[pI]->Write("", TObject::kOverwrite);
  }

  for(Int_t pI = 0; pI < nPartonPtBins; ++pI){  
    delete partonRRatio_h[pI];

    for(Int_t fI = 0; fI < nEnergyFracBins; ++fI){
      delete partonRRatio_Bins_h[pI][fI];
    }

    delete partonRRatio_Mean_h[pI];
  }

  for(Int_t pI = 0; pI < nJtPtBins; ++pI){  
    delete jtRRatio_h[pI];

    for(Int_t fI = 0; fI < nEnergyFracBins; ++fI){
      delete jtRRatio_Bins_h[pI][fI];
    }

    delete jtRRatio_Mean_h[pI];
  }

  outFile_p->Close();
  delete outFile_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/makePartonRRatio.exe <inFileName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += makePartonRRatio(argv[1]);
  return retVal;
}
