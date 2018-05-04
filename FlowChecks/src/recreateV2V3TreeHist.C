#include <iostream>
#include <string>
#include <vector>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TMath.h"
#include "TDatime.h"
#include "TF1.h"

#include "Utility/include/doGlobalDebug.h"

int recreateV2V3TreeHist(const std::string inFileName)
{
  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  TFile* outFile_p = new TFile(("output/v2CrossCheck_TreeHist_" + dateStr + ".root").c_str(), "RECREATE");

  const Int_t nCentBins = 11;
  const Int_t centBinsLow[nCentBins] = {5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55};
  const Int_t centBinsHi[nCentBins] = {10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60};
  const Double_t centBins[nCentBins+1] = {5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60};

  
  TH1F* v2Raw_h[nCentBins];
  TH1F* v2RawCorr_h[nCentBins];

  TH1F* v2Obs_h[nCentBins];
  TH1F* v2ObsCorr_h[nCentBins];

  TH1F* v2Raw_Mean_h = new TH1F("v2Raw_Mean_h", ";Centrality (%);#LTv_{2}^{raw}#GT", nCentBins, centBins);
  TH1F* v2Raw_Sigma_h = new TH1F("v2Raw_Sigma_h", ";Centrality (%);#sigma(v_{2}^{raw})", nCentBins, centBins);
  TH1F* v2RawCorr_Mean_h = new TH1F("v2RawCorr_Mean_h", ";Centrality (%);#LTv_{2}^{raw}#GT", nCentBins, centBins);
  TH1F* v2RawCorr_Sigma_h = new TH1F("v2RawCorr_Sigma_h", ";Centrality (%);#sigma(v_{2}^{raw})", nCentBins, centBins);

  TH1F* v2Obs_Mean_h = new TH1F("v2Obs_Mean_h", ";Centrality (%);#LTv_{2}^{obs}#GT", nCentBins, centBins);
  TH1F* v2Obs_Sigma_h = new TH1F("v2Obs_Sigma_h", ";Centrality (%);#sigma(v_{2}^{obs})", nCentBins, centBins);
  TH1F* v2ObsCorr_Mean_h = new TH1F("v2ObsCorr_Mean_h", ";Centrality (%);#LTv_{2}^{obs}#GT", nCentBins, centBins);
  TH1F* v2ObsCorr_Sigma_h = new TH1F("v2ObsCorr_Sigma_h", ";Centrality (%);#sigma(v_{2}^{obs})", nCentBins, centBins);

  for(Int_t cI = 0; cI < nCentBins; ++cI){
    v2Raw_h[cI] = new TH1F(("v2Raw_Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]) + "_h").c_str(), ";v_{2};Counts", 150, 0.0, 0.6);
    v2RawCorr_h[cI] = new TH1F(("v2RawCorr_Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]) + "_h").c_str(), ";v_{2};Counts", 150, 0.0, 0.6);
    v2Obs_h[cI] = new TH1F(("v2Obs_Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]) + "_h").c_str(), ";v_{2};Counts", 150, 0.0, 0.6);
    v2ObsCorr_h[cI] = new TH1F(("v2ObsCorr_Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]) + "_h").c_str(), ";v_{2};Counts", 150, 0.0, 0.6);
  }

  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  TTree* inTree_p = (TTree*)inFile_p->Get("v2V3Tree");

  Int_t hiBin_;
  Float_t hiEvt2Plane_;
  Float_t hiEvt3Plane_;
  Float_t v2FromTree_;
  std::vector<float>* pfPhi_p=NULL;
  std::vector<float>* pfWeight_p=NULL;

  inTree_p->SetBranchAddress("hiBin", &hiBin_);
  inTree_p->SetBranchAddress("hiEvt2Plane", &hiEvt2Plane_);
  inTree_p->SetBranchAddress("hiEvt3Plane", &hiEvt3Plane_);
  inTree_p->SetBranchAddress("v2FromTree", &v2FromTree_);
  inTree_p->SetBranchAddress("pfPhi", &pfPhi_p);
  inTree_p->SetBranchAddress("pfWeight", &pfWeight_p);

  const Int_t nEntries = TMath::Min((Int_t)inTree_p->GetEntries(), (Int_t)100000000);

  double globalN[nCentBins];
  double globalV2XRaw[nCentBins];
  double globalV2YRaw[nCentBins];
  double globalV2XRawCorr[nCentBins];
  double globalV2YRawCorr[nCentBins];
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    globalN[cI] = 0.;
    globalV2XRaw[cI] = 0.;
    globalV2YRaw[cI] = 0.;
    globalV2XRawCorr[cI] = 0.;
    globalV2YRawCorr[cI] = 0.;
  }

  for(Int_t entry = 0; entry < nEntries; ++entry){
    if(entry%10000 == 0) std::cout << " Entry " << entry << "/" << nEntries << std::endl;
    inTree_p->GetEntry(entry);

    int centPos = -1;
    for(Int_t cI = 0; cI < nCentBins; ++cI){
      if(centBinsLow[cI] <= hiBin_/2 && centBinsHi[cI] > hiBin_/2){
	centPos = cI;
	break;
      }
    }

    if(centPos < 0) continue;

    double v2xRaw = 0.;
    double v2yRaw = 0.;

    double v2xRawCorr = 0.;
    double v2yRawCorr = 0.;

    double weight = 0.;
  
    for(unsigned pfI = 0; pfI < pfPhi_p->size(); pfI++){
      double tempWeight = pfWeight_p->at(pfI);
      double deltaEventPhi = pfPhi_p->at(pfI) - hiEvt2Plane_;

      v2xRaw += TMath::Cos(2*(deltaEventPhi));
      v2yRaw += TMath::Sin(2*(deltaEventPhi));

      v2xRawCorr += tempWeight*TMath::Cos(2*(deltaEventPhi));
      v2yRawCorr += tempWeight*TMath::Sin(2*(deltaEventPhi));

      weight += tempWeight;
    }

    v2xRaw /= (double)pfPhi_p->size();
    v2yRaw /= (double)pfPhi_p->size();

    v2xRawCorr /= weight;
    v2yRawCorr /= weight;

    double v2Raw = TMath::Sqrt(v2xRaw*v2xRaw + v2yRaw*v2yRaw);
    double v2RawCorr = TMath::Sqrt(v2xRawCorr*v2xRawCorr + v2yRawCorr*v2yRawCorr);

    v2Raw_h[centPos]->Fill(v2Raw);
    v2RawCorr_h[centPos]->Fill(v2RawCorr);

    globalN[centPos] += 1;
    globalV2XRaw[centPos] += v2xRaw;
    globalV2YRaw[centPos] += v2yRaw;
    globalV2XRawCorr[centPos] += v2xRawCorr;
    globalV2YRawCorr[centPos] += v2yRawCorr;
  }

  for(Int_t cI = 0; cI < nCentBins; ++cI){
    globalV2XRaw[cI] /= globalN[cI];
    globalV2YRaw[cI] /= globalN[cI];
    globalV2XRawCorr[cI] /= globalN[cI];
    globalV2YRawCorr[cI] /= globalN[cI];
    std::cout << "Cent, N: " << centBinsLow[cI] << "-" << centBinsHi[cI] << "%, " << globalN[cI] << std::endl;
    std::cout << "  Global v2xRaw, v2yRaw, N: " << globalV2XRaw[cI] << ", " << globalV2YRaw[cI] << std::endl;
    std::cout << "  Global v2xRawCorr, v2yRawCorr, N: " << globalV2XRawCorr[cI] << ", " << globalV2YRawCorr[cI] << std::endl;
  }



  for(Int_t entry = 0; entry < nEntries; ++entry){
    if(entry%10000 == 0) std::cout << " Entry " << entry << "/" << nEntries << std::endl;
    inTree_p->GetEntry(entry);

    int centPos = -1;
    for(Int_t cI = 0; cI < nCentBins; ++cI){
      if(centBinsLow[cI] <= hiBin_/2 && centBinsHi[cI] > hiBin_/2){
	centPos = cI;
	break;
      }
    }

    if(centPos < 0) continue;

    double v2xObs = 0.;
    double v2yObs = 0.;

    double v2xObsCorr = 0.;
    double v2yObsCorr = 0.;

    double weight = 0.;
  
    for(unsigned pfI = 0; pfI < pfPhi_p->size(); pfI++){
      double tempWeight = pfWeight_p->at(pfI);
      double deltaEventPhi = pfPhi_p->at(pfI) - hiEvt2Plane_;

      v2xObs += TMath::Cos(2*(deltaEventPhi));
      v2yObs += TMath::Sin(2*(deltaEventPhi));

      v2xObsCorr += tempWeight*TMath::Cos(2*(deltaEventPhi));
      v2yObsCorr += tempWeight*TMath::Sin(2*(deltaEventPhi));

      weight += tempWeight;
    }

    v2xObs /= (double)pfPhi_p->size();
    v2yObs /= (double)pfPhi_p->size();

    v2xObsCorr /= weight;
    v2yObsCorr /= weight;

    v2xObs -= globalV2XRaw[centPos];
    v2yObs -= globalV2YRaw[centPos];

    v2xObsCorr -= globalV2XRawCorr[centPos];
    v2yObsCorr -= globalV2YRawCorr[centPos];

    double v2Obs = TMath::Sqrt(v2xObs*v2xObs + v2yObs*v2yObs);
    double v2ObsCorr = TMath::Sqrt(v2xObsCorr*v2xObsCorr + v2yObsCorr*v2yObsCorr);

    v2Obs_h[centPos]->Fill(v2Obs);
    v2ObsCorr_h[centPos]->Fill(v2ObsCorr);
  }


  inFile_p->Close();
  delete inFile_p;

  outFile_p->cd();

  for(Int_t cI = 0; cI < nCentBins; ++cI){  
    v2Raw_h[cI]->Write("", TObject::kOverwrite);
    v2RawCorr_h[cI]->Write("", TObject::kOverwrite);

    v2Obs_h[cI]->Write("", TObject::kOverwrite);
    v2ObsCorr_h[cI]->Write("", TObject::kOverwrite);

    v2Raw_Mean_h->SetBinContent(cI+1, v2Raw_h[cI]->GetMean());
    v2Raw_Mean_h->SetBinError(cI+1, v2Raw_h[cI]->GetMeanError());
    v2Raw_Sigma_h->SetBinContent(cI+1, v2Raw_h[cI]->GetStdDev());
    v2Raw_Sigma_h->SetBinError(cI+1, v2Raw_h[cI]->GetStdDevError());

    v2RawCorr_Mean_h->SetBinContent(cI+1, v2RawCorr_h[cI]->GetMean());
    v2RawCorr_Mean_h->SetBinError(cI+1, v2RawCorr_h[cI]->GetMeanError());
    v2RawCorr_Sigma_h->SetBinContent(cI+1, v2RawCorr_h[cI]->GetStdDev());
    v2RawCorr_Sigma_h->SetBinError(cI+1, v2RawCorr_h[cI]->GetStdDevError());

    v2Obs_Mean_h->SetBinContent(cI+1, v2Obs_h[cI]->GetMean());
    v2Obs_Mean_h->SetBinError(cI+1, v2Obs_h[cI]->GetMeanError());
    v2Obs_Sigma_h->SetBinContent(cI+1, v2Obs_h[cI]->GetStdDev());
    v2Obs_Sigma_h->SetBinError(cI+1, v2Obs_h[cI]->GetStdDevError());

    v2ObsCorr_Mean_h->SetBinContent(cI+1, v2ObsCorr_h[cI]->GetMean());
    v2ObsCorr_Mean_h->SetBinError(cI+1, v2ObsCorr_h[cI]->GetMeanError());
    v2ObsCorr_Sigma_h->SetBinContent(cI+1, v2ObsCorr_h[cI]->GetStdDev());
    v2ObsCorr_Sigma_h->SetBinError(cI+1, v2ObsCorr_h[cI]->GetStdDevError());

    delete v2Raw_h[cI];
    delete v2RawCorr_h[cI];

    delete v2Obs_h[cI];
    delete v2ObsCorr_h[cI];
  }

  v2Raw_Mean_h->Write("", TObject::kOverwrite);
  v2Raw_Sigma_h->Write("", TObject::kOverwrite);

  v2RawCorr_Mean_h->Write("", TObject::kOverwrite);
  v2RawCorr_Sigma_h->Write("", TObject::kOverwrite);

  v2Obs_Mean_h->Write("", TObject::kOverwrite);
  v2Obs_Sigma_h->Write("", TObject::kOverwrite);

  v2ObsCorr_Mean_h->Write("", TObject::kOverwrite);
  v2ObsCorr_Sigma_h->Write("", TObject::kOverwrite);

  outFile_p->Close();
  delete outFile_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "USAGE: ./recreateV2V3TreeHist.exe <inFileName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += recreateV2V3TreeHist(argv[1]);
  return retVal;
}
