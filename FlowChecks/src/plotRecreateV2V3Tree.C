#include <iostream>
#include <string>

#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TDatime.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"

#include "Utility/include/histDefUtility.h"
#include "Utility/include/kirchnerPalette.h"
#include "Utility/include/vanGoghPalette.h"


void setColorStyleLabelTitle(TH1F* inHist_p, const Int_t color, const Int_t style, const Int_t labelSize, const Int_t titleSize)
{
  inHist_p->SetLineColor(color);
  inHist_p->SetMarkerColor(color);
  inHist_p->SetMarkerStyle(style);

  inHist_p->GetXaxis()->SetTitleFont(43);
  inHist_p->GetYaxis()->SetTitleFont(43);
  inHist_p->GetXaxis()->SetLabelFont(43);
  inHist_p->GetYaxis()->SetLabelFont(43);

  inHist_p->GetXaxis()->SetTitleSize(titleSize);
  inHist_p->GetYaxis()->SetTitleSize(titleSize);
  inHist_p->GetXaxis()->SetLabelSize(labelSize);
  inHist_p->GetYaxis()->SetLabelSize(labelSize);

  return;
}


int plotRecreateV2V3(const std::string inFileName, const std::string inJamesFileName)
{
  vanGoghPalette vg;
  const Int_t flowStyleSet[3] = {33, 34, 21};
  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->SetTextFont(43);
  label_p->SetTextSize(18);

  TLegend* leg_p = new TLegend(0.6, 0.55, 0.9, 0.85);
  leg_p->SetTextFont(43);
  leg_p->SetTextSize(18);
  leg_p->SetBorderSize(0);
  leg_p->SetFillColor(0);
  leg_p->SetFillStyle(0);
  
  const Int_t nCentBins = 11;
  const Int_t centBinsLow[nCentBins] = {5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55};
  const Int_t centBinsHi[nCentBins] = {10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60};
  const Double_t centBins[nCentBins+1] = {5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60};

  TFile* jamesFile_p = new TFile(inJamesFileName.c_str(), "READ");
  TH1F* jamesHist_p[nCentBins];
  TH1F* jamesMean_p = new TH1F("jamesMean_p", ";Centrality (%);#LTv_{2}^{obs}#GT", nCentBins, centBins);
  TH1F* jamesSigma_p = new TH1F("jamesMean_p", ";Centrality (%);#sigma(v_{2}^{obs}", nCentBins, centBins);

  for(Int_t cI = 0; cI < nCentBins; ++cI){
    jamesHist_p[cI] = (TH1F*)jamesFile_p->Get(("qwebye/hVnFull_c" + std::to_string(cI+1)).c_str());

    jamesMean_p->SetBinContent(cI+1, jamesHist_p[cI]->GetMean());
    jamesMean_p->SetBinError(cI+1, jamesHist_p[cI]->GetMeanError());

    jamesSigma_p->SetBinContent(cI+1, jamesHist_p[cI]->GetStdDev());
    jamesSigma_p->SetBinError(cI+1, jamesHist_p[cI]->GetStdDevError());

    centerTitles(jamesHist_p[cI]);
    setSumW2(jamesHist_p[cI]);


    jamesHist_p[cI]->Scale(1./jamesHist_p[cI]->Integral());
    jamesHist_p[cI]->SetMaximum(jamesHist_p[cI]->GetMaximum()*1.5);

    setColorStyleLabelTitle(jamesHist_p[cI], 1, 20, 12, 14);
    jamesHist_p[cI]->GetYaxis()->SetTitleOffset(jamesHist_p[cI]->GetYaxis()->GetTitleOffset()*1.5);
    jamesHist_p[cI]->SetTitle("");

    if(cI == 0) leg_p->AddEntry(jamesHist_p[cI], "HIN-16-019", "P L");
  }

  centerTitles({jamesMean_p, jamesSigma_p});
  setSumW2({jamesMean_p, jamesSigma_p});
  jamesMean_p->SetMaximum(jamesMean_p->GetMaximum()*1.5);
  jamesMean_p->SetMinimum(0.0);

  setColorStyleLabelTitle(jamesMean_p, 1, 20, 12, 14);
  jamesMean_p->GetYaxis()->SetTitleOffset(jamesHist_p[0]->GetYaxis()->GetTitleOffset());
  
  jamesSigma_p->SetMaximum(jamesSigma_p->GetMaximum()*1.5);
  jamesSigma_p->SetMinimum(0.0);
  setColorStyleLabelTitle(jamesSigma_p, 1, 20, 12, 14);
  jamesSigma_p->GetYaxis()->SetTitleOffset(jamesHist_p[0]->GetYaxis()->GetTitleOffset());

  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  TH1F* v2Raw_h[nCentBins];
  TH1F* v2RawCorr_h[nCentBins];
  TH1F* v2ObsCorr_h[nCentBins];

  TH1F* v2Raw_Mean_h = (TH1F*)inFile_p->Get("v2Raw_Mean_h");
  TH1F* v2Raw_Sigma_h = (TH1F*)inFile_p->Get("v2Raw_Sigma_h");

  TH1F* v2RawCorr_Mean_h = (TH1F*)inFile_p->Get("v2RawCorr_Mean_h");
  TH1F* v2RawCorr_Sigma_h = (TH1F*)inFile_p->Get("v2RawCorr_Sigma_h");

  TH1F* v2ObsCorr_Mean_h = (TH1F*)inFile_p->Get("v2ObsCorr_Mean_h");
  TH1F* v2ObsCorr_Sigma_h = (TH1F*)inFile_p->Get("v2ObsCorr_Sigma_h");

  for(Int_t cI = 0; cI < nCentBins; ++cI){
    const std::string centStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);
    v2Raw_h[cI] = (TH1F*)inFile_p->Get(("v2Raw_" + centStr + "_h").c_str());
    v2RawCorr_h[cI] = (TH1F*)inFile_p->Get(("v2RawCorr_" + centStr + "_h").c_str());
    v2ObsCorr_h[cI] = (TH1F*)inFile_p->Get(("v2ObsCorr_" + centStr + "_h").c_str());

    centerTitles({v2Raw_h[cI], v2RawCorr_h[cI], v2ObsCorr_h[cI]});
    setSumW2({v2Raw_h[cI], v2RawCorr_h[cI], v2ObsCorr_h[cI]});

    setColorStyleLabelTitle(v2Raw_h[cI], vg.getColor(0), flowStyleSet[0], 12, 14);
    setColorStyleLabelTitle(v2RawCorr_h[cI], vg.getColor(1), flowStyleSet[1], 12, 14);
    setColorStyleLabelTitle(v2ObsCorr_h[cI], vg.getColor(2), flowStyleSet[2], 12, 14);

    v2Raw_h[cI]->Scale(1./v2Raw_h[cI]->Integral());
    v2RawCorr_h[cI]->Scale(1./v2RawCorr_h[cI]->Integral());
    v2ObsCorr_h[cI]->Scale(1./v2ObsCorr_h[cI]->Integral());
    
    v2Raw_h[cI]->GetYaxis()->SetTitleOffset(jamesHist_p[cI]->GetYaxis()->GetTitleOffset());
    v2Raw_h[cI]->GetXaxis()->SetTitleOffset(jamesHist_p[cI]->GetYaxis()->GetTitleOffset()*1.5);       

    if(cI == 0){
      leg_p->AddEntry(v2Raw_h[cI], "Jet Extraction", "P L");
      leg_p->AddEntry(v2RawCorr_h[cI], "Eff. Corr.", "P L");
      leg_p->AddEntry(v2ObsCorr_h[cI], "Eff.+Det. Corr.", "P L");
    }
  }

  centerTitles({v2Raw_Mean_h, v2Raw_Sigma_h, v2RawCorr_Mean_h, v2RawCorr_Sigma_h, v2ObsCorr_Mean_h, v2ObsCorr_Sigma_h});
  setSumW2({v2Raw_Mean_h, v2Raw_Sigma_h, v2RawCorr_Mean_h, v2RawCorr_Sigma_h, v2ObsCorr_Mean_h, v2ObsCorr_Sigma_h});

  setColorStyleLabelTitle(v2Raw_Mean_h, vg.getColor(0), flowStyleSet[0], 12, 14);
  setColorStyleLabelTitle(v2RawCorr_Mean_h, vg.getColor(1), flowStyleSet[1], 12, 14);
  setColorStyleLabelTitle(v2ObsCorr_Mean_h, vg.getColor(2), flowStyleSet[2], 12, 14);

  setColorStyleLabelTitle(v2Raw_Sigma_h, vg.getColor(0), flowStyleSet[0], 12, 14);
  setColorStyleLabelTitle(v2RawCorr_Sigma_h, vg.getColor(1), flowStyleSet[1], 12, 14);
  setColorStyleLabelTitle(v2ObsCorr_Sigma_h, vg.getColor(2), flowStyleSet[2], 12, 14);

  Float_t yPadBottomFrac = 0.35;

  for(Int_t cI = 0; cI < nCentBins; ++cI){
    TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);
    canv_p->SetTopMargin(0.01);
    canv_p->SetBottomMargin(0.01);
    canv_p->SetRightMargin(0.01);
    canv_p->SetLeftMargin(0.01);

    TPad* pads[2];
    pads[0] = new TPad("pad0", "", 0.0, yPadBottomFrac, 1.0, 1.0);
    canv_p->cd();
    pads[0]->SetRightMargin(0.01);
    pads[0]->SetBottomMargin(0.001);
    pads[0]->SetLeftMargin(pads[0]->GetLeftMargin()*1.3);
    pads[0]->SetTopMargin(pads[0]->GetLeftMargin()/(1.0-yPadBottomFrac)*3./6.);
    pads[0]->Draw("SAME");
    pads[0]->cd();
    
    jamesHist_p[cI]->DrawCopy("HIST E1 P");
    v2Raw_h[cI]->DrawCopy("HIST E1 P SAME");
    v2RawCorr_h[cI]->DrawCopy("HIST E1 P SAME");
    v2ObsCorr_h[cI]->DrawCopy("HIST E1 P SAME");

    const std::string centStr2 = std::to_string(centBinsLow[cI]) + "-" + std::to_string(centBinsHi[cI]) + "%";

    label_p->DrawLatex(0.15, 0.94, "#bf{CMS Preliminary}");
    label_p->DrawLatex(0.25, 0.82, ("#bf{" + centStr2 + "}").c_str());
    leg_p->Draw("SAME");

    pads[1] = new TPad("pad1", "", 0.0, 0.0, 1.0, yPadBottomFrac);
    canv_p->cd();
    pads[1]->SetTopMargin(0.001);
    pads[1]->SetRightMargin(0.01);
    pads[1]->SetBottomMargin(pads[0]->GetLeftMargin()/yPadBottomFrac);
    pads[1]->SetLeftMargin(pads[0]->GetLeftMargin());
    pads[1]->Draw("SAME");
    pads[1]->cd();


    v2Raw_h[cI]->Divide(jamesHist_p[cI]);
    v2RawCorr_h[cI]->Divide(jamesHist_p[cI]);
    v2ObsCorr_h[cI]->Divide(jamesHist_p[cI]);
    
    v2Raw_h[cI]->SetMinimum(0.55);
    v2Raw_h[cI]->SetMaximum(1.45);
    
    v2Raw_h[cI]->GetYaxis()->SetNdivisions(505);

    v2Raw_h[cI]->DrawCopy("HIST E1 P");
    v2RawCorr_h[cI]->DrawCopy("SAME HIST E1 P");
    v2ObsCorr_h[cI]->DrawCopy("SAME HIST E1 P");
    
    gStyle->SetOptStat(0);

    const std::string centStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);
    const std::string saveName = "pdfDir/plotRecreateV2James_" + centStr + "_" + dateStr + ".pdf";

    canv_p->SaveAs(saveName.c_str());

    delete pads[0];
    delete pads[1];
    delete canv_p;
  }

  inFile_p->Close();
  delete inFile_p;

  jamesFile_p->Close();
  delete jamesFile_p;


  delete leg_p;
  delete label_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 3){
    std::cout << "Usage: ./plotRecreateV2V3.exe <inFileName> <inJamesFileName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += plotRecreateV2V3(argv[1], argv[2]);
  return retVal;
}
