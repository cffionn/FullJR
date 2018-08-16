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
#include "TLegend.h"
#include "TMath.h"

//Non-local (Utility) dependencies
#include "Utility/include/checkMakeDir.h"
#include "Utility/include/histDefUtility.h"
#include "Utility/include/plotUtilities.h"
#include "Utility/include/vanGoghPalette.h"

//Local dependencies
#include "MCExample/include/atlasRAA.h"

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

  TH1D* genJtPt_AtlasBinnedR0p4_Total_h = (TH1D*)inFile_p->Get("genJtPt_AtlasBinnedR0p4_Total_h");
  TH1D* genJtPt_AtlasBinnedR0p4_Quark_h = (TH1D*)inFile_p->Get("genJtPt_AtlasBinnedR0p4_Quark_h");
  //  TH1D* genJtPt_AtlasBinnedR0p4_Gluon_h = (TH1D*)inFile_p->Get("genJtPt_AtlasBinnedR0p4_Gluon_h");
  
  TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);
  
  canv_p->SetTopMargin(0.07);
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

  logLoss_h->GetXaxis()->SetTitleOffset(logLoss_h->GetXaxis()->GetTitleOffset()*1.5);
  atlasFig_h->GetXaxis()->SetTitleOffset(atlasFig_h->GetXaxis()->GetTitleOffset()*1.5);

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

  label_p->DrawLatex(0.2, 0.95, "PYTHIA 8, 5.02 TeV + ATLAS RAA");

  delete label_p;

  TLatex* num_p = new TLatex();
  num_p->SetTextFont(43);
  num_p->SetTextSize(14);

  num_p->DrawLatex(200, -0.05, "200");
  num_p->DrawLatex(400, -0.05, "400");
  num_p->DrawLatex(800, -0.05, "800");

  delete num_p;

  TLegend* leg_p = new TLegend(0.3, 0.7, 0.6, 0.9);
  leg_p->SetBorderSize(0);
  leg_p->SetFillStyle(0);
  leg_p->SetFillColor(0);
  leg_p->SetTextFont(43);
  leg_p->SetTextSize(14);

  std::string tempStr = logLossStr;
  tempStr.replace(tempStr.find("p"), 1, ".");
  const double logLossDouble = std::stod(tempStr);

  std::cout << tempStr << ", " << logLossStr << ", " << logLossDouble << std::endl;

  leg_p->AddEntry(logLoss_h, ("Log Loss s=" + tempStr/*prettyString(logLossDouble, false, 3)*/).c_str(), "P L");
  leg_p->AddEntry(atlasFig_h, "ATLAS RAA", "P L");

  leg_p->Draw("SAME");

  checkMakeDir("pdfDir");
  const std::string saveName = "pdfDir/atlas" + raaDijet + "Comp_LogLoss_" + dateStr + ".pdf";
  quietSaveAs(canv_p, saveName);

  delete leg_p;
  delete canv_p;

  if(raaDijet.find("raa") == std::string::npos) return 0;

  canv_p = new TCanvas("canv_p", "", 450, 450);
  
  canv_p->SetTopMargin(0.07);
  canv_p->SetRightMargin(0.01);
  canv_p->SetLeftMargin(0.14);
  canv_p->SetBottomMargin(0.14);
  
  genJtPt_AtlasBinnedR0p4_Quark_h->Sumw2();
  genJtPt_AtlasBinnedR0p4_Total_h->Sumw2();

  genJtPt_AtlasBinnedR0p4_Quark_h->Divide(genJtPt_AtlasBinnedR0p4_Total_h);
  genJtPt_AtlasBinnedR0p4_Total_h->Divide(genJtPt_AtlasBinnedR0p4_Total_h);

  genJtPt_AtlasBinnedR0p4_Total_h->SetFillColor(vg.getColor(2));
  genJtPt_AtlasBinnedR0p4_Total_h->SetMarkerSize(0.01);

  genJtPt_AtlasBinnedR0p4_Quark_h->SetFillColor(vg.getColor(1));
  genJtPt_AtlasBinnedR0p4_Quark_h->SetMarkerSize(0.01);

  

  genJtPt_AtlasBinnedR0p4_Total_h->GetXaxis()->SetLabelFont(43);
  genJtPt_AtlasBinnedR0p4_Total_h->GetYaxis()->SetLabelFont(43);
  genJtPt_AtlasBinnedR0p4_Total_h->GetXaxis()->SetTitleFont(43);
  genJtPt_AtlasBinnedR0p4_Total_h->GetYaxis()->SetTitleFont(43);

  genJtPt_AtlasBinnedR0p4_Total_h->GetXaxis()->SetLabelSize(14);
  genJtPt_AtlasBinnedR0p4_Total_h->GetYaxis()->SetLabelSize(14);
  genJtPt_AtlasBinnedR0p4_Total_h->GetXaxis()->SetTitleSize(14);
  genJtPt_AtlasBinnedR0p4_Total_h->GetYaxis()->SetTitleSize(14);

  genJtPt_AtlasBinnedR0p4_Quark_h->GetXaxis()->SetLabelFont(43);
  genJtPt_AtlasBinnedR0p4_Quark_h->GetYaxis()->SetLabelFont(43);
  genJtPt_AtlasBinnedR0p4_Quark_h->GetXaxis()->SetTitleFont(43);
  genJtPt_AtlasBinnedR0p4_Quark_h->GetYaxis()->SetTitleFont(43);

  genJtPt_AtlasBinnedR0p4_Quark_h->GetXaxis()->SetLabelSize(14);
  genJtPt_AtlasBinnedR0p4_Quark_h->GetYaxis()->SetLabelSize(14);
  genJtPt_AtlasBinnedR0p4_Quark_h->GetXaxis()->SetTitleSize(14);
  genJtPt_AtlasBinnedR0p4_Quark_h->GetYaxis()->SetTitleSize(14);

  genJtPt_AtlasBinnedR0p4_Total_h->SetMaximum(1.2);
  genJtPt_AtlasBinnedR0p4_Total_h->SetMinimum(0.0);

  genJtPt_AtlasBinnedR0p4_Total_h->GetXaxis()->SetTitle("Gen. Jet p_{T} [GeV]");
  genJtPt_AtlasBinnedR0p4_Total_h->GetYaxis()->SetTitle("Flavor Fraction");

  centerTitles({genJtPt_AtlasBinnedR0p4_Total_h, genJtPt_AtlasBinnedR0p4_Quark_h});
  
  genJtPt_AtlasBinnedR0p4_Total_h->GetXaxis()->SetTitleOffset(genJtPt_AtlasBinnedR0p4_Total_h->GetXaxis()->GetTitleOffset()*1.5);
  genJtPt_AtlasBinnedR0p4_Quark_h->GetXaxis()->SetTitleOffset(genJtPt_AtlasBinnedR0p4_Quark_h->GetXaxis()->GetTitleOffset()*1.5);

  canv_p->cd();
  genJtPt_AtlasBinnedR0p4_Total_h->DrawCopy("HIST E1");
  genJtPt_AtlasBinnedR0p4_Quark_h->DrawCopy("HIST E1 SAME");
  gStyle->SetOptStat(0);
  gPad->RedrawAxis();
  if(raaDijet.find("raa") != std::string::npos) gPad->SetLogx();

  label_p = new TLatex();
  label_p->SetTextFont(43);
  label_p->SetTextSize(14);
  label_p->SetNDC();

  label_p->DrawLatex(0.2, 0.95, "PYTHIA 8, 5.02 TeV");

  delete label_p;

  num_p = new TLatex();
  num_p->SetTextFont(43);
  num_p->SetTextSize(14);

  num_p->DrawLatex(200, -0.05, "200");
  num_p->DrawLatex(400, -0.05, "400");
  num_p->DrawLatex(800, -0.05, "800");

  delete num_p;

  leg_p = new TLegend(0.3, 0.8, 0.6, 0.9);
  leg_p->SetBorderSize(0);
  leg_p->SetFillStyle(0);
  leg_p->SetFillColor(0);
  leg_p->SetTextFont(43);
  leg_p->SetTextSize(14);

  leg_p->AddEntry(genJtPt_AtlasBinnedR0p4_Total_h, "Gluon", "F");
  leg_p->AddEntry(genJtPt_AtlasBinnedR0p4_Quark_h, "Quark", "F");

  leg_p->Draw("SAME");


  const std::string saveNameQG = "pdfDir/qgFrac_LogLoss_" + dateStr + ".pdf";
  quietSaveAs(canv_p, saveNameQG);

  delete leg_p;
  delete canv_p;

  canv_p = new TCanvas("canv_p", "", 450, 450);
  
  canv_p->SetTopMargin(0.07);
  canv_p->SetRightMargin(0.01);
  canv_p->SetLeftMargin(0.14);
  canv_p->SetBottomMargin(0.14);

  const Int_t nBins = genJtPt_AtlasBinnedR0p4_Quark_h->GetNbinsX();
  Double_t bins[nBins+1];
  for(Int_t bI = 0; bI < nBins+1; ++bI){
    bins[bI] = genJtPt_AtlasBinnedR0p4_Quark_h->GetBinLowEdge(bI+1);
  }

  TH1D* meanEnergyLoss_h = new TH1D("meanEnergyLoss_h", ";Gen. Jet p_{T} [GeV];#LTEnergy Loss#GT", nBins, bins);
  
  for(Int_t bI = 0; bI < nBins; ++bI){
    Double_t loss = TMath::Log(TMath::E()*meanEnergyLoss_h->GetBinCenter(bI+1)/25.)*1.*logLossDouble*genJtPt_AtlasBinnedR0p4_Quark_h->GetBinContent(bI+1);
    loss += TMath::Log(TMath::E()*meanEnergyLoss_h->GetBinCenter(bI+1)/25.)*9./5.*logLossDouble*(1.-genJtPt_AtlasBinnedR0p4_Quark_h->GetBinContent(bI+1));

    meanEnergyLoss_h->SetBinContent(bI+1, loss);
    meanEnergyLoss_h->SetBinError(bI+1, 0);
  }


  meanEnergyLoss_h->Sumw2();
  meanEnergyLoss_h->SetMarkerColor(vg.getColor(1));
  meanEnergyLoss_h->SetMarkerStyle(20);
  meanEnergyLoss_h->SetMarkerSize(0.8);
  meanEnergyLoss_h->SetLineColor(vg.getColor(1));
  
  meanEnergyLoss_h->GetXaxis()->SetLabelFont(43);
  meanEnergyLoss_h->GetYaxis()->SetLabelFont(43);
  meanEnergyLoss_h->GetXaxis()->SetTitleFont(43);
  meanEnergyLoss_h->GetYaxis()->SetTitleFont(43);

  meanEnergyLoss_h->GetXaxis()->SetLabelSize(14);
  meanEnergyLoss_h->GetYaxis()->SetLabelSize(14);
  meanEnergyLoss_h->GetXaxis()->SetTitleSize(14);
  meanEnergyLoss_h->GetYaxis()->SetTitleSize(14);

  centerTitles({meanEnergyLoss_h});
  
  meanEnergyLoss_h->GetXaxis()->SetTitleOffset(meanEnergyLoss_h->GetXaxis()->GetTitleOffset()*1.5);


  canv_p->cd();

  meanEnergyLoss_h->SetMinimum(0);
  meanEnergyLoss_h->SetMaximum(45);
  meanEnergyLoss_h->DrawCopy("HIST E1 P");
  gStyle->SetOptStat(0);
  gPad->RedrawAxis();
  if(raaDijet.find("raa") != std::string::npos) gPad->SetLogx();

  label_p = new TLatex();
  label_p->SetTextFont(43);
  label_p->SetTextSize(14);
  label_p->SetNDC();

  label_p->DrawLatex(0.2, 0.95, "PYTHIA 8, 5.02 TeV");

  delete label_p;

  num_p = new TLatex();
  num_p->SetTextFont(43);
  num_p->SetTextSize(14);

  num_p->DrawLatex(200, -2, "200");
  num_p->DrawLatex(400, -2, "400");
  num_p->DrawLatex(800, -2, "800");

  delete num_p;

  const std::string saveNameELoss = "pdfDir/eLoss_LogLoss_" + dateStr + ".pdf";
  quietSaveAs(canv_p, saveNameELoss);

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
