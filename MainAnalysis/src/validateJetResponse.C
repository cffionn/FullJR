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

//Local dependencies                                                                                                                             
#include "MainAnalysis/include/cutPropagator.h"

//Non-local FullJR (Utility, etc.) dependencies                                                                                                      
#include "Utility/include/checkMakeDir.h"
#include "Utility/include/returnRootFileContentsList.h"
#include "Utility/include/histDefUtility.h"
#include "Utility/include/plotUtilities.h"
#include "Utility/include/vanGoghPalette.h"


int validateJetResponse(const std::string inResponseName, const bool doParaFills)
{
  if(!checkFile(inResponseName)){
    std::cout << "Given inResponseName \'" << inResponseName << "\' is invalid. return 1" << std::endl;
    return 1;
  }

  TFile* responseFile_p = new TFile(inResponseName.c_str(), "READ");
  std::vector<std::string> jetDirList = returnRootFileContentsList(responseFile_p, "TDirectoryFile", "JetAnalyzer");

  std::cout << "Validating " << jetDirList.size() << " jets..." << std::endl;
  for(unsigned int jI = 0; jI < jetDirList.size(); ++jI){
    std::cout << " " << jI << "/" << jetDirList.size() << ": " << jetDirList.at(jI) << std::endl;
  }

  std::string paraString = "";
  Int_t nBayes = 1;
  if(doParaFills){
    paraString = "_ParaFills";
    nBayes = 6;
  }

  

  cutPropagator cutProp;
  cutProp.Clean();
  cutProp.GetAllVarFromFile(responseFile_p);

  Int_t nCentBinsTemp = cutProp.GetNCentBins();
  std::vector<Int_t> centBinsLow = cutProp.GetCentBinsLow();
  std::vector<Int_t> centBinsHi = cutProp.GetCentBinsHi();

  Int_t nJtPtBinsTemp = cutProp.GetNJtPtBins();
  std::vector<Double_t> jtPtBinsTemp = cutProp.GetJtPtBins();

  Int_t nJtAbsEtaBinsTemp = cutProp.GetNJtAbsEtaBins();
  std::vector<Double_t> jtAbsEtaBinsLowTemp = cutProp.GetJtAbsEtaBinsLow();
  std::vector<Double_t> jtAbsEtaBinsHiTemp = cutProp.GetJtAbsEtaBinsHi();

  bool isResponsePP = cutProp.GetIsPP();

  Int_t nIDTemp = cutProp.GetNID();
  std::vector<std::string> idStr = cutProp.GetIdStr();
  std::vector<double> jtPfCHMFCutLow = cutProp.GetJtPfCHMFCutLow();
  std::vector<double> jtPfCHMFCutHi = cutProp.GetJtPfCHMFCutHi();
  std::vector<double> jtPfMUMFCutLow = cutProp.GetJtPfMUMFCutLow();
  std::vector<double> jtPfMUMFCutHi = cutProp.GetJtPfMUMFCutHi();


  const Int_t nCentBins = nCentBinsTemp;                                                                                                         
                                                                                                                                                 
  const Int_t nJtPtBins = nJtPtBinsTemp;                                                                                                         
  Double_t jtPtBins[nJtPtBins+1];                                                                                                                
  std::cout << "nJtPtBins: ";                                                                                                                    
  for(Int_t jI = 0; jI < nJtPtBins+1; ++jI){                                                                                                     
    jtPtBins[jI] = jtPtBinsTemp.at(jI);                                                                                                          
    std::cout << " " << jtPtBins[jI] << ",";                                                                                                     
  }                                                                                                                                              
  std::cout << std::endl;                                                                                                                        
                                                                                                                                                 
  const Int_t nJtAbsEtaBins = nJtAbsEtaBinsTemp;                                                                                                 
  Double_t jtAbsEtaBinsLow[nJtAbsEtaBins];                                                                                                       
  Double_t jtAbsEtaBinsHi[nJtAbsEtaBins];                                                                                                        
  std::cout << "nJtAbsEtaBins: ";                                                                                                                
  for(Int_t jI = 0; jI < nJtAbsEtaBins; ++jI){                                                                                                   
    jtAbsEtaBinsLow[jI] = jtAbsEtaBinsLowTemp.at(jI);                                                                                            
    jtAbsEtaBinsHi[jI] = jtAbsEtaBinsHiTemp.at(jI);                                                                                              
    std::cout << " " << jtAbsEtaBinsLow[jI] << "-" << jtAbsEtaBinsHi[jI] << ",";                                                                 
  }                                                                                                                                              
  std::cout << std::endl;                                                                                                                        
                                                                                                                                                 
  const Int_t nID = nIDTemp;     

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
      std::string centStr = "PP";
      if(!isResponsePP) centStr = "Cent" + std::to_string(centBinsLow.at(cI)) + "to" + std::to_string(centBinsHi.at(cI));

      for(Int_t idI = 0; idI < nID; ++idI){
	for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
	  const std::string jtAbsEtaStr = "AbsEta" + prettyString(jtAbsEtaBinsLow[aI], 1, true) + "to" + prettyString(jtAbsEtaBinsHi[aI], 1, true);

	  RooUnfoldResponse* rooRes_p = (RooUnfoldResponse*)responseFile_p->Get((dirName + "/rooResponse_" + dirName + "_" + centStr + "_" + idStr.at(idI) + "_" + jtAbsEtaStr + "_RecoTrunc_h").c_str());
	  TH1D* recoJtPt_RecoTrunc_h = (TH1D*)responseFile_p->Get((dirName + "/recoJtPt_" + dirName + "_" + centStr + "_" + idStr.at(idI) + "_" + jtAbsEtaStr + "_RecoTrunc" + paraString + "_h").c_str());
	  TH1D* genJtPt_h = (TH1D*)responseFile_p->Get((dirName + "/genJtPt_" + dirName + "_" + centStr + "_" + idStr.at(idI) + "_" + jtAbsEtaStr + paraString + "_h").c_str());

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

	  for(Int_t bI = 0; bI < nBayes; ++bI){
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
	    
	    bool doLogX = false;
	    if(recoJtPt_RecoTrunc_h->GetBinWidth(1)*3 < recoJtPt_RecoTrunc_h->GetBinWidth(recoJtPt_RecoTrunc_h->GetNbinsX()-1)) doLogX = true;

	    recoJtPt_RecoTrunc_h->DrawCopy("HIST E1 P");
	    genJtPt_h->DrawCopy("HIST E1 P SAME");
	    
	    gPad->SetLogy();
	    if(doLogX) gPad->SetLogx();

	    canv_p->cd();
	    pads_p[1]->cd();
	    
	    RooUnfoldBayes bayes(rooRes_p, recoJtPt_RecoTrunc_h, 1+bI, false, "name");
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
	    if(doLogX) gPad->SetLogx();

	    canv_p->cd();
	    pads_p[2]->cd();
	    
	    unfold_h->Divide(genJtPt_h);
	    unfold_h->SetMaximum(1.25);
	    unfold_h->SetMinimum(0.75);
	    unfold_h->DrawCopy("HIST E1 P");
	    if(doLogX) gPad->SetLogx();

	    canv_p->SaveAs(("pdfDir/" + dirName + "_" + centStr + "_" + idStr.at(idI) + "_" + jtAbsEtaStr + "_Bayes" + std::to_string(bI+1) + paraString + "_" + dateStr + ".pdf").c_str());
	  
	    for(Int_t pI = 0; pI < nPad; ++pI){
	      delete pads_p[pI];
	    }
	    
	    delete canv_p;
	  }
	}
      }
    }
  }

  responseFile_p->Close();
  delete responseFile_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/validateJetResponse.exe <inResponseName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += validateJetResponse(argv[1], false);
  retVal += validateJetResponse(argv[1], true);

  std::cout << "Job complete. Return " << retVal << "." << std::endl;
  return retVal;
}
