//cpp dependencies
#include <iostream>
#include <string>

//ROOT dependencies
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TDatime.h"
#include "TStyle.h"

//RooUnfold dependencies
#include "src/RooUnfoldResponse.h"
#include "src/RooUnfoldBayes.h"

//Non-local FullJR (Utility, etc.) dependencies                                                                                                      
#include "Utility/include/checkMakeDir.h"
#include "Utility/include/returnRootFileContentsList.h"
#include "Utility/include/histDefUtility.h"
#include "Utility/include/vanGoghPalette.h"


int validateJetResponse(const std::string inFileName)
{
  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  std::vector<std::string> jetDirList = returnRootFileContentsList(inFile_p, "TDirectoryFile", "JetAnalyzer");
  std::vector<std::string> cutDirList = returnRootFileContentsList(inFile_p, "TNamed", "");

  std::cout << "Validating " << jetDirList.size() << " jets..." << std::endl;
  for(unsigned int jI = 0; jI < jetDirList.size(); ++jI){
    std::cout << " " << jI << "/" << jetDirList.size() << ": " << jetDirList.at(jI) << std::endl;
  }

  Int_t nCentBins = -1;
  std::vector<Int_t> centBinsLow;
  std::vector<Int_t> centBinsHi;

  std::cout << "Using cuts: " << std::endl;
  for(unsigned int cI = 0; cI < cutDirList.size(); ++cI){
    std::cout << " " << cI << "/" << cutDirList.size() << ": " << cutDirList.at(cI) << std::endl;

    std::string tempStr = cutDirList.at(cI);
    while(tempStr.find("/") != std::string::npos){tempStr.replace(0, tempStr.find("/")+1, "");}

    if(tempStr.find("nCentBins") != std::string::npos && tempStr.size() == std::string("nCentBins").size()) nCentBins = std::stoi(((TNamed*)inFile_p->Get(cutDirList.at(cI).c_str()))->GetTitle());
    else if(tempStr.find("centBinsLow") != std::string::npos && tempStr.size() == std::string("centBinsLow").size()){
      std::string tempStr2 = ((TNamed*)inFile_p->Get(cutDirList.at(cI).c_str()))->GetTitle();
      while(tempStr2.find(",") != std::string::npos){
	centBinsLow.push_back(std::stoi(tempStr2.substr(0, tempStr2.find(","))));
	tempStr2.replace(0, tempStr2.find(",")+1, "");
      }
      if(tempStr2.size() != 0) centBinsLow.push_back(std::stoi(tempStr2));
    }
    else if(tempStr.find("centBinsHi") != std::string::npos && tempStr.size() == std::string("centBinsHi").size()){
      std::string tempStr2 = ((TNamed*)inFile_p->Get(cutDirList.at(cI).c_str()))->GetTitle();
      while(tempStr2.find(",") != std::string::npos){
	centBinsHi.push_back(std::stoi(tempStr2.substr(0, tempStr2.find(","))));
	tempStr2.replace(0, tempStr2.find(",")+1, "");
      }
      if(tempStr2.size() != 0) centBinsHi.push_back(std::stoi(tempStr2));
    }
  }

  std::cout << "nCentBins: " << nCentBins << std::endl;
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    std::cout << " " << cI << "/" << nCentBins << ": " << centBinsLow.at(cI) << "-" << centBinsHi.at(cI) << std::endl;
  }

  vanGoghPalette vg;
  const Int_t nStyles = 3;
  const Int_t styles[nStyles] = {20, 47, 34};
  const Int_t colors[nStyles] = {vg.getColor(0), vg.getColor(1), vg.getColor(2)};

  const Int_t nJets = jetDirList.size(); 

  const Int_t nX = 2;
  const Int_t nY = 1;

  const Int_t nPad = 3;
  const Double_t padXLow[nPad] = {0.0, 0.5, 0.5};
  const Double_t padXHi[nPad] = {0.5, 1.0, 1.0};
  const Double_t padYLow[nPad] = {0.0, 0.35, 0.0};
  const Double_t padYHi[nPad] = {1.0, 1.0, 0.35};
  
  const Double_t globalMargin = 0.06;
  
  const Double_t leftMargin[nPad] = {globalMargin*1./(padXHi[0] - padXLow[0]), globalMargin*1./(padXHi[1] - padXLow[1]), globalMargin*1./(padXHi[2] - padXLow[2])};
  const Double_t rightMargin[nPad] = {0.001, 0.001, 0.001};
  const Double_t topMargin[nPad] = {leftMargin[0], leftMargin[0]*(padYHi[0] - padYLow[0])/(padYHi[1] - padYLow[1]), 0.001};
  const Double_t bottomMargin[nPad] = {leftMargin[0], 0.001, leftMargin[0]*(padYHi[0] - padYLow[0])/(padYHi[2] - padYLow[2])};

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  checkMakeDir("pdfDir");

  for(Int_t jI = 0; jI < nJets; ++jI){
    std::string dirName = jetDirList.at(jI);
    dirName = dirName.substr(0, dirName.find("/"));

    for(Int_t cI = 0; cI < nCentBins; ++cI){
      const std::string centStr = "Cent" + std::to_string(centBinsLow.at(cI)) + "to" + std::to_string(centBinsHi.at(cI));
      RooUnfoldResponse* rooRes_p = (RooUnfoldResponse*)inFile_p->Get((dirName + "/rooResponse_" + dirName + "_" + centStr + "_RecoTrunc_h").c_str());
      TH1D* recoJtPt_RecoTrunc_h = (TH1D*)inFile_p->Get((dirName + "/recoJtPt_" + dirName + "_" + centStr + "_RecoTrunc_h").c_str());
      TH1D* genJtPt_h = (TH1D*)inFile_p->Get((dirName + "/genJtPt_" + dirName + "_" + centStr + "_h").c_str());

      recoJtPt_RecoTrunc_h->GetXaxis()->SetTitleFont(43);
      recoJtPt_RecoTrunc_h->GetYaxis()->SetTitleFont(43);
      recoJtPt_RecoTrunc_h->GetXaxis()->SetLabelFont(43);
      recoJtPt_RecoTrunc_h->GetYaxis()->SetLabelFont(43);

      recoJtPt_RecoTrunc_h->GetXaxis()->SetTitleSize(14);
      recoJtPt_RecoTrunc_h->GetYaxis()->SetTitleSize(14);
      recoJtPt_RecoTrunc_h->GetXaxis()->SetLabelSize(14);
      recoJtPt_RecoTrunc_h->GetYaxis()->SetLabelSize(14);

      genJtPt_h->GetXaxis()->SetTitleFont(43);
      genJtPt_h->GetYaxis()->SetTitleFont(43);
      genJtPt_h->GetXaxis()->SetLabelFont(43);
      genJtPt_h->GetYaxis()->SetLabelFont(43);

      genJtPt_h->GetXaxis()->SetTitleSize(14);
      genJtPt_h->GetYaxis()->SetTitleSize(14);
      genJtPt_h->GetXaxis()->SetLabelSize(14);
      genJtPt_h->GetYaxis()->SetLabelSize(14);
     
      TCanvas* canv_p = new TCanvas("canv_p", "", 450*nX, 450*nY);
      canv_p->SetTopMargin(0.01);
      canv_p->SetBottomMargin(0.01);
      canv_p->SetLeftMargin(0.01);
      canv_p->SetRightMargin(0.01);

      gStyle->SetOptStat(0);

      TPad* pads_p[nPad];
      for(Int_t pI = 0; pI < nPad; ++pI){
	canv_p->cd();
	pads_p[pI] = new TPad(("pad" + std::to_string(pI)).c_str(), "", padXLow[pI], padYLow[pI], padXHi[pI], padYHi[pI]);
	canv_p->cd();
	pads_p[pI]->Draw("SAME");
	gStyle->SetOptStat(0);
	pads_p[pI]->SetLeftMargin(leftMargin[pI]);
	pads_p[pI]->SetRightMargin(rightMargin[pI]);
	pads_p[pI]->SetTopMargin(topMargin[pI]);
	pads_p[pI]->SetBottomMargin(bottomMargin[pI]);

	std::cout << "Creating pad " << pI << "..." << std::endl;
	std::cout << " Dimensions: x " << padXLow[pI] << "-" << padXHi[pI] << ", y " << padYLow[pI] << "-" << padYHi[pI] << std::endl;
	std::cout << " Margins: left " << leftMargin[pI] << ", right " << rightMargin[pI] << ", top " << topMargin[pI] << ", bottom " << bottomMargin[pI] << std::endl;
      }

      canv_p->cd();
      pads_p[0]->cd();

      recoJtPt_RecoTrunc_h->SetMarkerColor(colors[0]);
      recoJtPt_RecoTrunc_h->SetLineColor(colors[0]);
      recoJtPt_RecoTrunc_h->SetMarkerStyle(styles[0]);
      recoJtPt_RecoTrunc_h->SetMarkerSize(0.8);

      genJtPt_h->SetMarkerColor(colors[1]);
      genJtPt_h->SetLineColor(colors[1]);
      genJtPt_h->SetMarkerStyle(styles[1]);
      genJtPt_h->SetMarkerSize(0.8);

      recoJtPt_RecoTrunc_h->GetXaxis()->SetNdivisions(505);
      recoJtPt_RecoTrunc_h->GetYaxis()->SetNdivisions(505);

      
      Double_t maxVal = 0.0;
      Double_t minVal = 999999.;

      for(Int_t bI = 0; bI < recoJtPt_RecoTrunc_h->GetNbinsX(); ++bI){
	if(recoJtPt_RecoTrunc_h->GetBinContent(bI+1) > maxVal) maxVal = recoJtPt_RecoTrunc_h->GetBinContent(bI+1);
	if(genJtPt_h->GetBinContent(bI+1) > maxVal) maxVal = genJtPt_h->GetBinContent(bI+1);

	if(recoJtPt_RecoTrunc_h->GetBinContent(bI+1) < minVal && recoJtPt_RecoTrunc_h->GetBinContent(bI+1) > 0) minVal = recoJtPt_RecoTrunc_h->GetBinContent(bI+1);
	if(genJtPt_h->GetBinContent(bI+1) < minVal && genJtPt_h->GetBinContent(bI+1) > 0) minVal = genJtPt_h->GetBinContent(bI+1);
      }

      maxVal *= 5.;
      minVal /= 5.;

      recoJtPt_RecoTrunc_h->SetMaximum(maxVal);
      recoJtPt_RecoTrunc_h->SetMinimum(minVal);

      recoJtPt_RecoTrunc_h->DrawCopy("HIST E1 P");
      genJtPt_h->DrawCopy("HIST E1 P SAME");

      gPad->SetLogy();

      canv_p->cd();
      pads_p[1]->cd();

      RooUnfoldBayes bayes(rooRes_p, recoJtPt_RecoTrunc_h, 1, false, "name");
      bayes.SetVerbose(0);

      TH1D* unfold_h = (TH1D*)bayes.Hreco(RooUnfold::kCovToy);

      unfold_h->SetMarkerColor(colors[2]);
      unfold_h->SetLineColor(colors[2]);
      unfold_h->SetMarkerStyle(styles[2]);
      unfold_h->SetMarkerSize(0.8);

      unfold_h->GetXaxis()->SetNdivisions(505);
      unfold_h->GetYaxis()->SetNdivisions(505);
      centerTitles(unfold_h);
      setSumW2(unfold_h);
      unfold_h->SetTitle("");

      unfold_h->GetXaxis()->SetTitleFont(43);
      unfold_h->GetYaxis()->SetTitleFont(43);
      unfold_h->GetXaxis()->SetLabelFont(43);
      unfold_h->GetYaxis()->SetLabelFont(43);

      unfold_h->GetXaxis()->SetTitleSize(14);
      unfold_h->GetYaxis()->SetTitleSize(14);
      unfold_h->GetXaxis()->SetLabelSize(14);
      unfold_h->GetYaxis()->SetLabelSize(14);

      unfold_h->SetMaximum(maxVal);
      unfold_h->SetMinimum(minVal);

      unfold_h->DrawCopy("HIST E1 P");
      genJtPt_h->DrawCopy("HIST E1 P SAME");

      gPad->SetLogy();

      canv_p->cd();
      pads_p[2]->cd();

      unfold_h->Divide(genJtPt_h);
      unfold_h->SetMaximum(1.25);
      unfold_h->SetMinimum(0.75);
      unfold_h->DrawCopy("HIST E1 P");

      canv_p->SaveAs(("pdfDir/" + dirName + "_" + centStr + "_" + dateStr + ".pdf").c_str());

      for(Int_t pI = 0; pI < nPad; ++pI){
	delete pads_p[pI];
      }

      delete canv_p;
    }
  }
  

  inFile_p->Close();
  delete inFile_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/validateJetResponse.exe <inFileName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += validateJetResponse(argv[1]);
  return retVal;
}
