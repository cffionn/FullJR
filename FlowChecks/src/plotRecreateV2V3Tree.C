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
#include "Utility/include/returnRootFileContentsList.h"
#include "Utility/include/vectorStringUtility.h"

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


int plotRecreateV2V3(const std::string inFileName, const std::string inJamesFileName, const std::string pfOrTrkStr)
{
  if(pfOrTrkStr.size() == 2 && pfOrTrkStr.find("PF") != std::string::npos) std::cout << "Running in mode: " << pfOrTrkStr << std::endl;
  else if(pfOrTrkStr.size() == 3 && pfOrTrkStr.find("Trk") != std::string::npos) std::cout << "Running in mode: " << pfOrTrkStr << std::endl;
  else{
    std::cout << "Given pfOrTrkStr \'" << pfOrTrkStr << "\' is invalid return 1" << std::endl;
    return 1;
  }

  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  std::vector<std::string> listOfHist = returnRootFileContentsList(inFile_p, "TH1F");
  removeVectorDuplicates(&listOfHist);
  bool hasPFOrTrkStr = false;
  for(unsigned int vI = 0; vI < listOfHist.size(); ++vI){
    if(listOfHist.at(vI).find("_"+pfOrTrkStr+"_") != std::string::npos){
      hasPFOrTrkStr = true;
      break;
    }
  }
  
  inFile_p->Close();
  delete inFile_p;

  if(!hasPFOrTrkStr){
    std::cout << "File \'" << inFileName << "\' does not contain \'" << pfOrTrkStr << "\'. return 1" << std::endl;
    return 1;
  }

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
  TH1F* jamesSigma_p = new TH1F("jamesSigma_p", ";Centrality (%);#sigma(v_{2}^{obs}", nCentBins, centBins);

  Double_t globalMax = 0.0;

  for(Int_t cI = 0; cI < nCentBins; ++cI){
    jamesHist_p[cI] = (TH1F*)jamesFile_p->Get(("qwebye/hVnFull_c" + std::to_string(cI+1)).c_str());

    jamesMean_p->SetBinContent(cI+1, jamesHist_p[cI]->GetMean());
    jamesMean_p->SetBinError(cI+1, jamesHist_p[cI]->GetMeanError());

    jamesSigma_p->SetBinContent(cI+1, jamesHist_p[cI]->GetStdDev());
    jamesSigma_p->SetBinError(cI+1, jamesHist_p[cI]->GetStdDevError());

    centerTitles(jamesHist_p[cI]);
    setSumW2(jamesHist_p[cI]);


    jamesHist_p[cI]->Scale(1./jamesHist_p[cI]->Integral());
    for(Int_t bIX = 0; bIX < jamesHist_p[cI]->GetNbinsX(); ++bIX){      
      jamesHist_p[cI]->SetBinContent(bIX+1, jamesHist_p[cI]->GetBinContent(bIX+1)/jamesHist_p[cI]->GetBinWidth(bIX+1));
      jamesHist_p[cI]->SetBinError(bIX+1, jamesHist_p[cI]->GetBinError(bIX+1)/jamesHist_p[cI]->GetBinWidth(bIX+1));
    }

    if(jamesHist_p[cI]->GetMaximum() > globalMax) globalMax = jamesHist_p[cI]->GetMaximum();

    jamesHist_p[cI]->SetMaximum(jamesHist_p[cI]->GetMaximum()*1.5);

    setColorStyleLabelTitle(jamesHist_p[cI], 1, 20, 12, 14);
    jamesHist_p[cI]->GetYaxis()->SetTitleOffset(jamesHist_p[cI]->GetYaxis()->GetTitleOffset()*1.3);
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

  inFile_p = new TFile(inFileName.c_str(), "READ");
  TH1F* v2Raw_h[nCentBins];
  TH1F* v2RawCorr_h[nCentBins];
  TH1F* v2ObsCorr_h[nCentBins];

  TH1F* v2Raw_Mean_h = (TH1F*)inFile_p->Get(("v2Raw_Mean_" + pfOrTrkStr + "_h").c_str());
  TH1F* v2Raw_Sigma_h = (TH1F*)inFile_p->Get(("v2Raw_Sigma_" + pfOrTrkStr + "_h").c_str());

  TH1F* v2RawCorr_Mean_h = (TH1F*)inFile_p->Get(("v2RawCorr_Mean_" + pfOrTrkStr + "_h").c_str());
  TH1F* v2RawCorr_Sigma_h = (TH1F*)inFile_p->Get(("v2RawCorr_Sigma_" + pfOrTrkStr + "_h").c_str());

  TH1F* v2ObsCorr_Mean_h = (TH1F*)inFile_p->Get(("v2ObsCorr_Mean_" + pfOrTrkStr + "_h").c_str());
  TH1F* v2ObsCorr_Sigma_h = (TH1F*)inFile_p->Get(("v2ObsCorr_Sigma_" + pfOrTrkStr + "_h").c_str());

  for(Int_t cI = 0; cI < nCentBins; ++cI){
    const std::string centStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);
    v2Raw_h[cI] = (TH1F*)inFile_p->Get(("v2Raw_" + centStr + "_" + pfOrTrkStr + "_h").c_str());
    v2RawCorr_h[cI] = (TH1F*)inFile_p->Get(("v2RawCorr_" + centStr + "_" + pfOrTrkStr + "_h").c_str());
    v2ObsCorr_h[cI] = (TH1F*)inFile_p->Get(("v2ObsCorr_" + centStr + "_" + pfOrTrkStr + "_h").c_str());

    centerTitles({v2Raw_h[cI], v2RawCorr_h[cI], v2ObsCorr_h[cI]});
    setSumW2({v2Raw_h[cI], v2RawCorr_h[cI], v2ObsCorr_h[cI]});

    setColorStyleLabelTitle(v2Raw_h[cI], vg.getColor(0), flowStyleSet[0], 12, 14);
    setColorStyleLabelTitle(v2RawCorr_h[cI], vg.getColor(1), flowStyleSet[1], 12, 14);
    setColorStyleLabelTitle(v2ObsCorr_h[cI], vg.getColor(2), flowStyleSet[2], 12, 14);

    v2Raw_h[cI]->Scale(1./v2Raw_h[cI]->Integral());
    v2RawCorr_h[cI]->Scale(1./v2RawCorr_h[cI]->Integral());
    v2ObsCorr_h[cI]->Scale(1./v2ObsCorr_h[cI]->Integral());

    for(Int_t bIX = 0; bIX < v2Raw_h[cI]->GetNbinsX(); ++bIX){      
      v2Raw_h[cI]->SetBinContent(bIX+1, v2Raw_h[cI]->GetBinContent(bIX+1)/v2Raw_h[cI]->GetBinWidth(bIX+1));
      v2Raw_h[cI]->SetBinError(bIX+1, v2Raw_h[cI]->GetBinError(bIX+1)/v2Raw_h[cI]->GetBinWidth(bIX+1));
    }

    for(Int_t bIX = 0; bIX < v2RawCorr_h[cI]->GetNbinsX(); ++bIX){      
      v2RawCorr_h[cI]->SetBinContent(bIX+1, v2RawCorr_h[cI]->GetBinContent(bIX+1)/v2RawCorr_h[cI]->GetBinWidth(bIX+1));
      v2RawCorr_h[cI]->SetBinError(bIX+1, v2RawCorr_h[cI]->GetBinError(bIX+1)/v2RawCorr_h[cI]->GetBinWidth(bIX+1));
    }

    for(Int_t bIX = 0; bIX < v2ObsCorr_h[cI]->GetNbinsX(); ++bIX){      
      v2ObsCorr_h[cI]->SetBinContent(bIX+1, v2ObsCorr_h[cI]->GetBinContent(bIX+1)/v2ObsCorr_h[cI]->GetBinWidth(bIX+1));
      v2ObsCorr_h[cI]->SetBinError(bIX+1, v2ObsCorr_h[cI]->GetBinError(bIX+1)/v2ObsCorr_h[cI]->GetBinWidth(bIX+1));
    }

    if(v2Raw_h[cI]->GetMaximum() > globalMax) globalMax = v2Raw_h[cI]->GetMaximum();
    if(v2RawCorr_h[cI]->GetMaximum() > globalMax) globalMax = v2RawCorr_h[cI]->GetMaximum();
    if(v2ObsCorr_h[cI]->GetMaximum() > globalMax) globalMax = v2ObsCorr_h[cI]->GetMaximum();

    if(cI == 0){
      leg_p->AddEntry(v2Raw_h[cI], "No Corr.", "P L");
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

  std::vector<TH1*> jHist_p;
  std::vector<TH1*> v2RawHist_p;
  std::vector<TH1*> v2RawCorrHist_p;
  std::vector<TH1*> v2ObsCorrHist_p;
  std::vector<std::string> centStrName;
  std::vector<std::string> centStrLabel;
  std::vector<double> maxVal;
  std::vector<double> minVal;
  std::vector<std::string> yTitle;
  std::vector<std::string> yTitleRat;

  for(Int_t cI = 0; cI < nCentBins; ++cI){
    jHist_p.push_back(jamesHist_p[cI]);
    v2RawHist_p.push_back(v2Raw_h[cI]);
    v2RawCorrHist_p.push_back(v2RawCorr_h[cI]);
    v2ObsCorrHist_p.push_back(v2ObsCorr_h[cI]);
    centStrLabel.push_back(std::to_string(centBinsLow[cI]) + "-" + std::to_string(centBinsHi[cI]) + "%");
    centStrName.push_back("Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]));
    maxVal.push_back(globalMax);
    minVal.push_back(0.0);
    yTitle.push_back("#frac{1}{N_{evt}} #frac{dN_{evt}}{dv_{2}}");
    yTitleRat.push_back("#frac{v_{2}}{HIN-16-019}");
  }

  jHist_p.push_back(jamesMean_p);
  v2RawHist_p.push_back(v2Raw_Mean_h);
  v2RawCorrHist_p.push_back(v2RawCorr_Mean_h);
  v2ObsCorrHist_p.push_back(v2ObsCorr_Mean_h);
  centStrLabel.push_back(std::to_string(centBinsLow[0]) + "-" + std::to_string(centBinsHi[nCentBins-1]) + "%");
  centStrName.push_back("Mean_Cent" + std::to_string(centBinsLow[0]) + "to" + std::to_string(centBinsHi[nCentBins-1]));
  maxVal.push_back(0.20);
  minVal.push_back(0.0);
  yTitle.push_back("#LTv_{2}#GT");
  yTitleRat.push_back("#frac{#LTv_{2}#GT}{HIN-16-019}");

  jHist_p.push_back(jamesSigma_p);
  v2RawHist_p.push_back(v2Raw_Sigma_h);
  v2RawCorrHist_p.push_back(v2RawCorr_Sigma_h);
  v2ObsCorrHist_p.push_back(v2ObsCorr_Sigma_h);
  centStrLabel.push_back(std::to_string(centBinsLow[0]) + "-" + std::to_string(centBinsHi[nCentBins-1]) + "%");
  centStrName.push_back("Sigma_Cent" + std::to_string(centBinsLow[0]) + "to" + std::to_string(centBinsHi[nCentBins-1]));
  maxVal.push_back(0.10);
  minVal.push_back(0.0);
  yTitle.push_back("#sigma(v_{2})");
  yTitleRat.push_back("#frac{#sigma(v_{2})}{HIN-16-019}");
  
  Float_t yPadBottomFrac = 0.35;

  for(unsigned int cI = 0; cI < jHist_p.size(); ++cI){
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

    jHist_p.at(cI)->GetYaxis()->SetTitle(yTitle.at(cI).c_str());

    jHist_p.at(cI)->SetMaximum(maxVal.at(cI)*1.1);
    jHist_p.at(cI)->SetMinimum(minVal.at(cI));


    jHist_p.at(cI)->GetYaxis()->SetNdivisions(505);
    jHist_p.at(cI)->DrawCopy("HIST E1 P");
    v2RawHist_p.at(cI)->DrawCopy("HIST E1 P SAME");
    v2RawCorrHist_p.at(cI)->DrawCopy("HIST E1 P SAME");
    v2ObsCorrHist_p.at(cI)->DrawCopy("HIST E1 P SAME");

    gPad->RedrawAxis();

    label_p->DrawLatex(0.15, 0.94, "#bf{CMS Preliminary}");
    label_p->DrawLatex(0.25, 0.82, ("#bf{" + centStrLabel.at(cI) + "}").c_str());
    label_p->DrawLatex(0.25, 0.74, ("#bf{Using " + pfOrTrkStr + "}").c_str());
    leg_p->Draw("SAME");

    pads[1] = new TPad("pad1", "", 0.0, 0.0, 1.0, yPadBottomFrac);
    canv_p->cd();
    pads[1]->SetTopMargin(0.001);
    pads[1]->SetRightMargin(0.01);
    pads[1]->SetBottomMargin(pads[0]->GetLeftMargin()/yPadBottomFrac);
    pads[1]->SetLeftMargin(pads[0]->GetLeftMargin());
    pads[1]->Draw("SAME");
    pads[1]->cd();


    v2RawHist_p.at(cI)->Divide(jHist_p.at(cI));
    v2RawCorrHist_p.at(cI)->Divide(jHist_p.at(cI));
    v2ObsCorrHist_p.at(cI)->Divide(jHist_p.at(cI));
    
    v2RawHist_p.at(cI)->SetMinimum(0.55);
    v2RawHist_p.at(cI)->SetMaximum(1.45);
    
    v2RawHist_p.at(cI)->GetYaxis()->SetNdivisions(505);


    v2RawHist_p.at(cI)->GetYaxis()->SetTitle(yTitleRat.at(cI).c_str());

    v2RawHist_p.at(cI)->GetYaxis()->SetTitleOffset(jHist_p.at(cI)->GetYaxis()->GetTitleOffset());
    v2RawHist_p.at(cI)->GetXaxis()->SetTitleOffset(jHist_p.at(cI)->GetYaxis()->GetTitleOffset()*1.8);       

    v2RawHist_p.at(cI)->DrawCopy("HIST E1 P");
    v2RawCorrHist_p.at(cI)->DrawCopy("SAME HIST E1 P");
    v2ObsCorrHist_p.at(cI)->DrawCopy("SAME HIST E1 P");

    gPad->RedrawAxis();
    
    gStyle->SetOptStat(0);

    const std::string saveName = "pdfDir/plotRecreateV2James_" + centStrName.at(cI) + "_" + pfOrTrkStr + "_" + dateStr + ".pdf";

    canv_p->SaveAs(saveName.c_str());

    delete pads[0];
    delete pads[1];
    delete canv_p;
  }

  inFile_p->Close();
  delete inFile_p;

  delete jamesMean_p;
  delete jamesSigma_p;

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
  retVal += plotRecreateV2V3(argv[1], argv[2], "PF");
  retVal += plotRecreateV2V3(argv[1], argv[2], "Trk");
  return retVal;
}
