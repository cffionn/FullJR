#include <iostream>
#include <string>

#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TDatime.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TMath.h"

#include "Utility/include/checkMakeDir.h"
#include "Utility/include/histDefUtility.h"
#include "Utility/include/plotUtilities.h"
#include "Utility/include/vanGoghPalette.h"

#include "MCExcample/include/atlasRAA.h"

int plotRAAAndDijet(const std::string inFileName, const std::string raaDijet, const std::string logLossStr)
{
  atlasRAA raa;
  vanGoghPalette vg;

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  TH1D* atlasFig_h = NULL;
  if(raaDijet.find("raa") != std::string::npos) atlasFig_h = (TH1D*)inFile_p->Get("atlasRAA_R0p4_h");
  else if(raaDijet.find("dijet") != std::string::npos) atlasFig_h = (TH1D*)inFile_p->Get("atlasDijetXJ_R0p4_h");
  raa.StyleHistogram(atlasFig_h);
  TH1D* logLoss_h = NULL;
  if(raaDijet.find("raa") != std::string::npos) logLoss_h = (TH1D*)inFile_p->Get(("raa_AtlasBinnedR0p4_LogLoss" + logLossStr + "_h").c_str());
  else if(raaDijet.find("dijet") != std::string::npos) logLoss_h = (TH1D*)inFile_p->Get(("dijetXJ_AtlasBinnedR0p4_LogLoss" + logLossStr + "_h").c_str());
  
  TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);
  
  canv_p->SetTopMargin(0.01);
  canv_p->SetRightMargin(0.01);
  canv_p->SetLeftMargin(0.14);
  canv_p->SetBottomMargin(0.14);
  
  logLoss_h->SetMarkerStyle(20);
  logLoss_h->SetMarkerSize(0.8);
  logLoss_h->SetMarkerColor(vg.getColor(0));
  logLoss_h->SetLineColor(vg.getColor(0));

  logLoss_h->GetXaxis()->SetLabelFont(43);
  logLoss_h->GetYaxis()->SetLabelFont(43);
  logLoss_h->GetXaxis()->SetTitleFont(43);
  logLoss_h->GetYaxis()->SetTitleFont(43);

  logLoss_h->GetXaxis()->SetLabelSize(14);
  logLoss_h->GetYaxis()->SetLabelSize(14);
  logLoss_h->GetXaxis()->SetTitleSize(14);
  logLoss_h->GetYaxis()->SetTitleSize(14);

  atlasFig_h->GetXaxis()->SetLabelFont(43);
  atlasFig_h->GetYaxis()->SetLabelFont(43);
  atlasFig_h->GetXaxis()->SetTitleFont(43);
  atlasFig_h->GetYaxis()->SetTitleFont(43);

  atlasFig_h->GetXaxis()->SetLabelSize(14);
  atlasFig_h->GetYaxis()->SetLabelSize(14);
  atlasFig_h->GetXaxis()->SetTitleSize(14);
  atlasFig_h->GetYaxis()->SetTitleSize(14);

  atlasFig_h->SetMaximum(1.5*TMath::Max(atlasFig_h->GetMaximum(), logLoss_h->GetMaximum()));
  atlasFig_h->SetMinimum(0.0);

  if(raaDijet.find("raa") != std::string::npos){
    logLoss_h->GetXaxis()->SetTitle("Jet p_{T} [GeV]");
    logLoss_h->GetYaxis()->SetTitle("R_{AA}");
    
    atlasFig_h->GetXaxis()->SetTitle("Jet p_{T} [GeV]");
    atlasFig_h->GetYaxis()->SetTitle("R_{AA}");
  }
  else{
    logLoss_h->GetXaxis()->SetTitle("A_{J}");
    logLoss_h->GetYaxis()->SetTitle("#frac{1}{N_{J}}#frac{dN_{J}}{dA_{J}}");

    atlasFig_h->GetXaxis()->SetTitle("A_{J}");
    atlasFig_h->GetYaxis()->SetTitle("#frac{1}{N_{J}}#frac{dN_{J}}{dA_{J}}");
  }

  centerTitles({logLoss_h, atlasFig_h});
  
  canv_p->cd();
  atlasFig_h->DrawCopy("HIST E1 P");
  logLoss_h->DrawCopy("HIST E1 P SAME");
  gStyle->SetOptStat(0);
  if(raaDijet.find("raa") != std::string::npos) gPad->SetLogx();

  TLatex* label_p = new TLatex();
  label_p->SetTextFont(43);
  label_p->SetTextSize(14);
  label_p->SetNDC();

  delete label_p;

  checkMakeDir("pdfDir");
  const std::string saveName = "pdfDir/atlas" + raaDijet + "Comp_LogLoss_" + dateStr + ".pdf";
  canv_p->SaveAs(saveName.c_str());

  delete canv_p;

  inFile_p->Close();
  delete inFile_p;
    
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 3){
    std::cout << "Usage: ./bin/plotRAAAndDijet.exe <inFileName> <logLosStr>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += plotRAAAndDijet(argv[1], "raa", argv[2]);
  retVal += plotRAAAndDijet(argv[1], "dijet", argv[2]);
  return retVal;
}
