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
#include "TPad.h"
#include "TStyle.h"
#include "TTree.h"

//RooUnfold dependencies
#include "src/RooUnfoldBayes.h"
#include "src/RooUnfoldResponse.h"

//Non-local dependencies
#include "Utility/include/checkMakeDir.h"
#include "Utility/include/histDefUtility.h"
#include "Utility/include/kirchnerPalette.h"
#include "Utility/include/plotUtilities.h"

void setLabels(TH1D* hist_p)
{
  hist_p->GetXaxis()->SetTitleFont(43);
  hist_p->GetYaxis()->SetTitleFont(43);
  hist_p->GetXaxis()->SetLabelFont(43);
  hist_p->GetYaxis()->SetLabelFont(43);

  hist_p->GetXaxis()->SetTitleSize(14);
  hist_p->GetYaxis()->SetTitleSize(14);
  hist_p->GetXaxis()->SetLabelSize(14);
  hist_p->GetYaxis()->SetLabelSize(14);

  return;
}

void setLabels(TH2D* hist_p)
{
  hist_p->GetXaxis()->SetTitleFont(43);
  hist_p->GetYaxis()->SetTitleFont(43);
  hist_p->GetXaxis()->SetLabelFont(43);
  hist_p->GetYaxis()->SetLabelFont(43);

  hist_p->GetXaxis()->SetTitleSize(14);
  hist_p->GetYaxis()->SetTitleSize(14);
  hist_p->GetXaxis()->SetLabelSize(14);
  hist_p->GetYaxis()->SetLabelSize(14);
  
  return;
}

void oldToNewPrior(TH1D* distToUnfold_p, TH1D* reco_p, TH2D* response_p, TH1D* newPrior_p)
{
  newPrior_p->Scale(distToUnfold_p->Integral()/newPrior_p->Integral());  
  for(Int_t bIY = 0; bIY < response_p->GetYaxis()->GetNbins(); ++bIY){
    Double_t scale = 0.0;
    
    for(Int_t bIX = 0; bIX < response_p->GetXaxis()->GetNbins(); ++bIX){
      scale += response_p->GetBinContent(bIX+1, bIY+1);
    }

    scale = newPrior_p->GetBinContent(bIY+1)/scale;

    for(Int_t bIX = 0; bIX < response_p->GetXaxis()->GetNbins(); ++bIX){
      response_p->SetBinContent(bIX+1, bIY+1, response_p->GetBinContent(bIX+1, bIY+1)*scale);
      response_p->SetBinError(bIX+1, bIY+1, response_p->GetBinError(bIX+1, bIY+1)*scale);
    }    
  }

  for(Int_t bIX = 0; bIX < response_p->GetXaxis()->GetNbins(); ++bIX){
    Double_t scale = 0.0;

    for(Int_t bIY = 0; bIY < response_p->GetYaxis()->GetNbins(); ++bIY){
      scale += response_p->GetBinContent(bIX+1, bIY+1);
    }

    Double_t val = reco_p->GetBinContent(bIX+1);
    Double_t err = reco_p->GetBinError(bIX+1);
    reco_p->SetBinContent(bIX+1, scale);
    reco_p->SetBinError(bIX+1, err*scale/val);
  }

  return;
}


int unfoldTest(const std::string inFileName)
{
  if(!checkFile(inFileName) || inFileName.find(".root") == std::string::npos){
    std::cout << "inFileName \'" << inFileName << "\' is invalid. return 1" << std::endl;
    return 1;
  }

  kirchnerPalette kPal;
  const Int_t nStyles = 5;
  const Int_t styles[nStyles] = {24, 25, 27, 28, 46};
  const Int_t nBayes = 5;
  
  const std::string responseStr = "response_SmallBins_ak3PFJetAnalyzer_PP_Cent0to10_LightMUAndCHID_ResponseMod0p10_AbsEta0p0to2p0_RecoGenAsymm_h";
  const std::string responsePriorStr = "response_SmallBins_ak3PFJetAnalyzer_PP_Cent0to10_LightMUAndCHID_ResponseMod0p10_AbsEta0p0to2p0_PriorFlat_RecoGenAsymm_h";

  const std::string genStr = "genJtPt_SmallBins_ak3PFJetAnalyzer_PP_Cent0to10_LightMUAndCHID_ResponseMod0p10_AbsEta0p0to2p0_GoodReco_h";
  const std::string recoStr = "recoJtPt_ak3PFJetAnalyzer_PP_Cent0to10_LightMUAndCHID_ResponseMod0p10_AbsEta0p0to2p0_GoodGen_h";

  const std::string genPriorStr = "genJtPt_SmallBins_ak3PFJetAnalyzer_PP_Cent0to10_LightMUAndCHID_ResponseMod0p10_AbsEta0p0to2p0_PriorFlat_GoodReco_h";
  const std::string recoPriorStr = "recoJtPt_ak3PFJetAnalyzer_PP_Cent0to10_LightMUAndCHID_ResponseMod0p10_AbsEta0p0to2p0_PriorFlat_GoodGen_h";

  const std::string genParaStr = "genJtPt_SmallBins_ak3PFJetAnalyzer_PP_Cent0to10_LightMUAndCHID_ResponseMod0p10_AbsEta0p0to2p0_GoodReco_h";
  const std::string recoParaStr = "recoJtPt_ak3PFJetAnalyzer_PP_Cent0to10_LightMUAndCHID_ResponseMod0p10_AbsEta0p0to2p0_RecoTrunc_h";

  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");

  TH2D* response_p = (TH2D*)inFile_p->Get(responseStr.c_str());
  TH2D* responsePrior_p = (TH2D*)inFile_p->Get(responsePriorStr.c_str());

  TH1D* gen_p = (TH1D*)inFile_p->Get(genStr.c_str());
  TH1D* genPrior_p = (TH1D*)inFile_p->Get(genPriorStr.c_str());

  TH1D* reco_p = (TH1D*)inFile_p->Get(recoStr.c_str());
  TH1D* recoPrior_p = (TH1D*)inFile_p->Get(recoPriorStr.c_str());

  TH1D* genPara_p = (TH1D*)inFile_p->Get(genParaStr.c_str());
  TH1D* recoPara_p = (TH1D*)inFile_p->Get(recoParaStr.c_str());
  
  const Int_t nPads = 10;
  const Int_t nXWidths = 4;
  const Int_t nYWidths = 2;
  TCanvas* canv_p = new TCanvas("canv_p", "", 400*nXWidths, 400*nYWidths);
  canv_p->SetTopMargin(0.01);
  canv_p->SetRightMargin(0.01);
  canv_p->SetLeftMargin(0.01);
  canv_p->SetBottomMargin(0.01);

  TPad* pads_p[nPads];

  const Double_t xLow[nPads] = {0.0, 0.25, 0.5, 0.75, 0.0, 0.25, 0.5, 0.75, 0.0, 0.25};
  const Double_t xUp[nPads] = {0.25, 0.5, 0.75, 1.0, 0.25, 0.5, 0.75, 1.0, 0.25, 0.5};
  const Double_t yLow[nPads] = {0.7, 0.7 , 0.5, 0.5, 0.0, 0.0, 0.0, 0.0, 0.5, 0.5};
  const Double_t yUp[nPads] = {1.0, 1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 0.5, 0.7, 0.7};

  const Double_t lMargin[nPads] = {0.1, 0.1, 0.14, 0.14, 0.14, 0.14, 0.14, 0.1, 0.1, 0.1};
  const Double_t rMargin[nPads] = {0.01, 0.01, 0.14, 0.14, 0.14, 0.14, 0.14, 0.01, 0.01, 0.01};
  const Double_t tMargin[nPads] = {0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};
  const Double_t bMargin[nPads] = {0.01, 0.01, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 0.2};

  for(unsigned int pI = 0; pI < nPads; ++pI){        
    std::cout << pI << ": " << xLow[pI] << ", " << xUp[pI] << ", " << yLow[pI] << ", " << yUp[pI] << std::endl;

    pads_p[pI] = new TPad(("pads_" + std::to_string(pI)).c_str(), "", xLow[pI], yLow[pI], xUp[pI], yUp[pI]);
    pads_p[pI]->SetTopMargin(tMargin[pI]);
    pads_p[pI]->SetBottomMargin(bMargin[pI]);
    pads_p[pI]->SetLeftMargin(lMargin[pI]);
    pads_p[pI]->SetRightMargin(rMargin[pI]);

    pads_p[pI]->Draw();
  }

  canv_p->cd();
  pads_p[2]->cd();
  centerTitles(responsePrior_p);
  responsePrior_p->SetMaximum(100);
  responsePrior_p->SetMinimum(0.00001);
  setLabels(responsePrior_p);
  responsePrior_p->DrawCopy("COLZ");
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->SetLogz();

  canv_p->cd();
  pads_p[3]->cd();
  centerTitles(response_p);
  response_p->SetMaximum(100);
  response_p->SetMinimum(0.00001);
  setLabels(response_p);
  response_p->DrawCopy("COLZ");  
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->SetLogz();

  canv_p->cd();
  pads_p[0]->cd();
  
  recoPara_p->SetMarkerSize(0);
  recoPara_p->SetMarkerColor(1);
  recoPara_p->SetLineColor(1);

  genPara_p->SetMarkerStyle(20);
  genPara_p->SetMarkerSize(1);
  genPara_p->SetMarkerColor(1);
  genPara_p->SetLineColor(1);

  centerTitles(recoPara_p);
  
  Double_t maxVal = getNearestFactor10Up(recoPara_p->GetMaximum(), 2);
  Double_t minVal = getNearestFactor10Down(recoPara_p->GetMinimum(), 2);

  recoPara_p->SetMaximum(maxVal);
  recoPara_p->SetMinimum(minVal);

  setLabels(recoPara_p);
  recoPara_p->DrawCopy("HIST E1");
  genPara_p->DrawCopy("HIST E1 P SAME");
  gStyle->SetOptStat(0);
  gPad->SetLogy();
  gPad->SetLogx();

  canv_p->cd();
  pads_p[1]->cd();
  
  recoPara_p->DrawCopy("HIST E1");
  gStyle->SetOptStat(0);
  gPad->SetLogy();
  gPad->SetLogx();
    
  for(Int_t bI = 0; bI < nBayes; ++bI){
    TH1D* preGenClone_p = (TH1D*)genPrior_p->Clone((std::string(gen_p->GetName()) + "_PRECLONE").c_str());
    TH1D* preRecoClone_p = (TH1D*)recoPrior_p->Clone((std::string(reco_p->GetName()) + "_PRECLONE").c_str());
    TH2D* preResponseClone_p = (TH2D*)responsePrior_p->Clone((std::string(responsePrior_p->GetName()) + "_PRECLONE").c_str());
    TH1D* preRecoToUnfoldClone_p = (TH1D*)recoPara_p->Clone((std::string(recoPara_p->GetName()) + "_PRECLONE").c_str());
    RooUnfoldResponse* preRooRes_p = new RooUnfoldResponse(preRecoClone_p, preGenClone_p, preResponseClone_p);

    if(bI == 0){
      canv_p->cd();
      pads_p[4]->cd();
      preResponseClone_p->SetMaximum(100);
      preResponseClone_p->SetMinimum(0.00001);
      setLabels(preResponseClone_p);
      preResponseClone_p->DrawCopy("COLZ");
      gPad->SetLogx();
      gPad->SetLogy();
      gPad->SetLogz();
    }    
    
    RooUnfoldBayes preBayes(preRooRes_p, preRecoToUnfoldClone_p, 2, false, "name");
    preBayes.SetVerbose(-1);
    TH1D* preUnfold_h = (TH1D*)preBayes.Hreco(RooUnfold::kCovToy);      
    
    //    TH1D* postGenClone_p = (TH1D*)gen_p->Clone((std::string(gen_p->GetName()) + "_POSTCLONE").c_str());
    TH1D* postGenClone_p = (TH1D*)preUnfold_h->Clone((std::string(gen_p->GetName()) + "_POSTCLONE").c_str());
    TH1D* postRecoClone_p = (TH1D*)recoPrior_p->Clone((std::string(reco_p->GetName()) + "_POSTCLONE").c_str());
    TH2D* postResponseClone_p = (TH2D*)responsePrior_p->Clone((std::string(responsePrior_p->GetName()) + "_POSTCLONE").c_str());
    TH1D* postRecoToUnfoldClone_p = (TH1D*)recoPara_p->Clone((std::string(recoPara_p->GetName()) + "_POSTCLONE").c_str());

    std::cout << "TH2 DIMENSIONS: " << std::endl;
    std::cout << "X: ";
    for(Int_t bIX = 0; bIX < postResponseClone_p->GetXaxis()->GetNbins(); ++bIX){
      std::cout << postResponseClone_p->GetXaxis()->GetBinLowEdge(bIX+1) << ", ";
    }
    std::cout << postResponseClone_p->GetXaxis()->GetBinLowEdge(postResponseClone_p->GetXaxis()->GetNbins()+1) << std::endl;

    std::cout << "Y: ";
    for(Int_t bIY = 0; bIY < postResponseClone_p->GetYaxis()->GetNbins(); ++bIY){
      std::cout << postResponseClone_p->GetYaxis()->GetBinLowEdge(bIY+1) << ", ";
    }
    std::cout << postResponseClone_p->GetYaxis()->GetBinLowEdge(postResponseClone_p->GetYaxis()->GetNbins()+1) << std::endl;

    oldToNewPrior(postRecoToUnfoldClone_p, postRecoClone_p, postResponseClone_p, postGenClone_p);
    
    if(bI == 4){
      canv_p->cd();
      pads_p[5]->cd();
      postResponseClone_p->SetMaximum(100);
      postResponseClone_p->SetMinimum(0.00001);
      setLabels(postResponseClone_p);
      postResponseClone_p->DrawCopy("COLZ");
      gPad->SetLogx();
      gPad->SetLogy();
      gPad->SetLogz();
    }

    RooUnfoldResponse* postRooRes_p = new RooUnfoldResponse(postRecoClone_p, postGenClone_p, postResponseClone_p);
         
    RooUnfoldBayes bayes(postRooRes_p, postRecoToUnfoldClone_p, bI+1, false, "name");
    bayes.SetVerbose(-1);

    TH1D* unfold_h = (TH1D*)bayes.Hreco(RooUnfold::kCovToy);
    TH1D* unfoldClone_h = (TH1D*)unfold_h->Clone((std::string(unfold_h->GetName()) + "_CLONE").c_str());
    
    canv_p->cd();
    pads_p[0]->cd();

    unfold_h->SetMarkerSize(1);
    unfold_h->SetMarkerStyle(styles[bI%nStyles]);
    unfold_h->SetMarkerColor(kPal.getColor(bI%nBayes));
    unfold_h->SetLineColor(kPal.getColor(bI%nBayes));
    setLabels(unfold_h);
    unfold_h->DrawCopy("HIST E1 P SAME");
    
    canv_p->cd();
    pads_p[8]->cd();
    unfold_h->Divide(genPara_p);

    if(bI == 0){
      unfold_h->SetMaximum(1.5);
      unfold_h->SetMinimum(0.5);
      centerTitles(unfold_h);
      unfold_h->SetTitle("");
      unfold_h->Print("ALL");
      unfold_h->DrawCopy("HIST E1 P");
      gPad->SetLogx();
    }
    else unfold_h->DrawCopy("HIST E1 P SAME");

    canv_p->cd();
    pads_p[1]->cd();

    TH1D* refold_h = (TH1D*)postRooRes_p->ApplyToTruth(unfoldClone_h);
    refold_h->SetMarkerSize(1);
    refold_h->SetMarkerStyle(styles[bI%nStyles]);
    refold_h->SetMarkerColor(kPal.getColor(bI%nBayes));
    refold_h->SetLineColor(kPal.getColor(bI%nBayes));  
    
    refold_h->DrawCopy("HIST E1 P SAME");
    refold_h->Divide(recoPara_p);

    canv_p->cd();
    pads_p[9]->cd();

    refold_h->SetMinimum(0.5);
    refold_h->SetMaximum(1.5);
    refold_h->SetTitle("");
    centerTitles(refold_h);

    setLabels(refold_h);
    if(bI == 0) refold_h->DrawCopy("HIST E1 P");
    else refold_h->DrawCopy("HIST E1 P SAME");
    gPad->SetLogx();


    if(bI == 4){
      canv_p->cd();
      pads_p[6]->cd();      
      postResponseClone_p->Divide(response_p);
      postResponseClone_p->SetMaximum(2.);
      postResponseClone_p->SetMinimum(0);
      postResponseClone_p->DrawCopy("COLZ");
      gPad->SetLogx();
      gPad->SetLogy();
    }

    delete refold_h;
    delete unfold_h;
    delete unfoldClone_h;

    delete postRooRes_p;
    
    delete postResponseClone_p;
    delete postGenClone_p;
    delete postRecoClone_p;
    
    delete preUnfold_h;

    delete preRooRes_p;
    delete preRecoToUnfoldClone_p;
    
    delete preResponseClone_p;
    delete preGenClone_p;
    delete preRecoClone_p;
  }

  canv_p->SaveAs("temp.pdf");

  for(Int_t pI = 0; pI < nPads; ++pI){
    delete pads_p[pI];
  }
  
  delete canv_p;

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
