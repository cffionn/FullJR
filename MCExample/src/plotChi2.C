//cpp dependencies
#include <iostream>
#include <string>

//ROOT dependencies
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TDatime.h"
#include "TStyle.h"
#include "TLatex.h"

//Non-local (Utility) dependencies
#include "Utility/include/checkMakeDir.h"
#include "Utility/include/histDefUtility.h"
#include "Utility/include/plotUtilities.h"

int plotChi2(const std::string inFileName)
{
  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  TH1D* atlasRAAChi2_LogLoss_h = (TH1D*)inFile_p->Get("atlasRAAChi2_LogLoss_h");
  
  TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);
  
  canv_p->SetTopMargin(0.01);
  canv_p->SetRightMargin(0.01);
  canv_p->SetLeftMargin(0.14);
  canv_p->SetBottomMargin(0.14);
  
  atlasRAAChi2_LogLoss_h->SetMarkerStyle(20);
  atlasRAAChi2_LogLoss_h->SetMarkerSize(0.8);
  atlasRAAChi2_LogLoss_h->SetMarkerColor(1);
  atlasRAAChi2_LogLoss_h->SetLineColor(1);

  atlasRAAChi2_LogLoss_h->GetXaxis()->SetLabelFont(43);
  atlasRAAChi2_LogLoss_h->GetYaxis()->SetLabelFont(43);
  atlasRAAChi2_LogLoss_h->GetXaxis()->SetTitleFont(43);
  atlasRAAChi2_LogLoss_h->GetYaxis()->SetTitleFont(43);

  atlasRAAChi2_LogLoss_h->GetXaxis()->SetLabelSize(14);
  atlasRAAChi2_LogLoss_h->GetYaxis()->SetLabelSize(14);
  atlasRAAChi2_LogLoss_h->GetXaxis()->SetTitleSize(14);
  atlasRAAChi2_LogLoss_h->GetYaxis()->SetTitleSize(14);

  atlasRAAChi2_LogLoss_h->GetXaxis()->SetTitle("Energy Loss Param [GeV]");
  atlasRAAChi2_LogLoss_h->GetYaxis()->SetTitle("RAA Chi2");

  Double_t minChi2 = 1000;
  Int_t minChi2Pos = -1;

  for(Int_t bIX = 0; bIX < atlasRAAChi2_LogLoss_h->GetNbinsX(); ++bIX){
    if(atlasRAAChi2_LogLoss_h->GetBinContent(bIX+1) < minChi2){
      minChi2 = atlasRAAChi2_LogLoss_h->GetBinContent(bIX+1);
      minChi2Pos = bIX+1;
    }
  }

  canv_p->cd();
  atlasRAAChi2_LogLoss_h->DrawCopy("HIST E1 P");
  gStyle->SetOptStat(0);

  TLatex* label_p = new TLatex();
  label_p->SetTextFont(43);
  label_p->SetTextSize(14);
  label_p->SetNDC();

  label_p->DrawLatex(0.3, 0.8, ("Minimum at " + prettyString(atlasRAAChi2_LogLoss_h->GetBinCenter(minChi2Pos), 4, false)).c_str());

  delete label_p;

  checkMakeDir("pdfDir");
  const std::string saveName = "pdfDir/atlasRAAChi2_LogLoss_" + dateStr + ".pdf";
  quietSaveAs(canv_p, saveName);

  delete canv_p;

  inFile_p->Close();
  delete inFile_p;
    
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/plotChi2.exe <inFileName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += plotChi2(argv[1]);
  return retVal;
}
