//cpp dependencies
#include <iostream>
#include <string>
#include <vector>

//ROOT dependencies
#include "TCanvas.h"
#include "TDatime.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TTree.h"

//RooUnfold dependencies
#include "src/RooUnfoldBayes.h"
#include "src/RooUnfoldResponse.h"

//Non-local dependencies
#include "Utility/include/checkMakeDir.h"

int unfoldTest(const std::string inFileName)
{
  if(!checkFile(inFileName) || inFileName.find(".root") == std::string::npos){
    std::cout << "inFileName \'" << inFileName << "\' is invalid. return 1" << std::endl;
    return 1;
  }

  const std::string responseStr = "response_SmallBins_ak3PFJetAnalyzer_PP_Cent0to10_LightMUAndCHID_ResponseMod0p10_AbsEta0p0to2p0_RecoGenAsymm_h";
  const std::string responsePriorStr = "response_SmallBins_ak3PFJetAnalyzer_PP_Cent0to10_LightMUAndCHID_ResponseMod0p10_AbsEta0p0to2p0_PriorFlat_RecoGenAsymm_h";
  const std::string genStr = "genJtPt_SmallBins_ak3PFJetAnalyzer_PP_Cent0to10_LightMUAndCHID_ResponseMod0p10_AbsEta0p0to2p0_h";
  const std::string recoStr = "recoJtPt_ak3PFJetAnalyzer_PP_Cent0to10_LightMUAndCHID_ResponseMod0p10_AbsEta0p0to2p0_RecoTrunc_h";

  const Int_t nBayes = 6;
  
  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");

  TH2D* response_p = (TH2D*)inFile_p->Get(responseStr.c_str());
  TH2D* responsePrior_p = (TH2D*)inFile_p->Get(responsePriorStr.c_str());
  TH1D* gen_p = (TH1D*)inFile_p->Get(genStr.c_str());
  TH1D* reco_p = (TH1D*)inFile_p->Get(recoStr.c_str());
  
  RooUnfoldResponse* rooRes_p = new RooUnfoldResponse(reco_p, gen_p, response_p);
  RooUnfoldResponse* rooResPrior_p = new RooUnfoldResponse(reco_p, gen_p, responsePrior_p);

  TCanvas* canv_p = new TCanvas("canv_p", "", 400*4, 400*2);
  canv_p->SetTopMargin(0.01);
  canv_p->SetRightMargin(0.01);
  canv_p->SetLeftMargin(0.01);
  canv_p->SetBottomMargin(0.01);
  canv_p->Divide(4, 2);

  canv_p->cd();
  canv_p->cd(1);

  reco_p->SetMarkerSize(0);
  reco_p->DrawCopy("HIST E1");
  gStyle->SetOptStat(0);
  gPad->SetLogy();
  
  for(Int_t bI = 0; bI < nBayes; ++bI){
    TH1D* genClone_p = (TH1D*)gen_p->Clone((std::string(gen_p->GetName()) + "_CLONE").c_str());
    TH1D* recoClone_p = (TH1D*)reco_p->Clone((std::string(reco_p->GetName()) + "_CLONE").c_str());

    RooUnfoldBayes bayes(rooRes_p, recoClone_p, bI+1, false, "name");
    bayes.SetVerbose(-1);
    TH1D* unfold_h = (TH1D*)bayes.Hreco(RooUnfold::kCovToy);

    
    
    delete unfold_h;
    
    delete genClone_p;
    delete recoClone_p;
  }

  canv_p->SaveAs("temp.pdf");
  delete canv_p;
  
  delete rooRes_p;
  delete rooResPrior_p;

  inFile_p->Close();
  delete inFile_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/unfoldTest.exe <inFileName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += unfoldTest(argv[1]);
  return retVal;
}
