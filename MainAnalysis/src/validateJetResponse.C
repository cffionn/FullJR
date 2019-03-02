//cpp dependenciGGGSjtptes
#include <iostream>
#include <string>
#include <vector>

//ROOT dependencies
#include "TCanvas.h"
#include "TDatime.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TPad.h"
#include "TStyle.h"

//RooUnfold dependencies
#include "src/RooUnfoldBayes.h"
#include "src/RooUnfoldResponse.h"

//Local dependencies
#include "MainAnalysis/include/cutPropagator.h"
#include "MainAnalysis/include/texSlideCreator.h"

//Non-local FullJR (Utility, etc.) dependencies
#include "Utility/include/checkMakeDir.h"
#include "Utility/include/getLinBins.h"
#include "Utility/include/getLogBins.h"
#include "Utility/include/histDefUtility.h"
#include "Utility/include/plotUtilities.h"
#include "Utility/include/returnRootFileContentsList.h"
#include "Utility/include/vanGoghPalette.h"

int validateJetResponse(const std::string inResponseName, const bool doParaFills, std::vector<std::string>* slideTitles, std::vector<std::vector<std::string> >* pdfPerSlide)
{
  if(!checkFile(inResponseName)){
    std::cout << "Given inResponseName \'" << inResponseName << "\' is invalid. return 1" << std::endl;
    return 1;
  }

  TFile* responseFile_p = new TFile(inResponseName.c_str(), "READ");
  std::vector<std::string> jetDirList = returnRootFileContentsList(responseFile_p, "TDirectoryFile", "JetAnalyzer");

  std::cout << "Validating " << jetDirList.size() << " jets..." << std::endl;
  for(unsigned int jI = 0; jI < jetDirList.size(); ++jI){
    std::cout << " " << jI << "/" << jetDirList.size() << ": " << jetDirList[jI] << std::endl;
  }

  std::string paraString = "";
  Int_t nBayes = 1;
  if(doParaFills){
    paraString = "_ParaFills";
    nBayes = 4;
  }

 
  cutPropagator cutProp;
  cutProp.Clean();
  cutProp.GetAllVarFromFile(responseFile_p);

  Int_t nCentBinsTemp = cutProp.GetNCentBins();
  std::vector<Int_t> centBinsLow = cutProp.GetCentBinsLow();
  std::vector<Int_t> centBinsHi = cutProp.GetCentBinsHi();

  /*
  Int_t nJtAbsEtaBinsTemp = cutProp.GetNJtAbsEtaBins();
  std::vector<Double_t> jtAbsEtaBinsLowTemp = cutProp.GetJtAbsEtaBinsLow();
  std::vector<Double_t> jtAbsEtaBinsHiTemp = cutProp.GetJtAbsEtaBinsHi();
  */

  Int_t nJtAbsEtaBinsTemp = 1;
  std::vector<Double_t> jtAbsEtaBinsLowTemp = {0.0};
  std::vector<Double_t> jtAbsEtaBinsHiTemp = {2.0};

  bool isResponsePP = cutProp.GetIsPP();

  /*
  Int_t nIDTemp = cutProp.GetNID();
  std::vector<std::string> idStr = cutProp.GetIdStr();
  */

  Int_t nIDTemp = 1;
  std::vector<std::string> idStr = {"LightMUAndCHID"};

  const Int_t nMaxCentBins = 4;
  const Int_t nCentBins = nCentBinsTemp;

  if(nCentBins > nMaxCentBins){
    std::cout << "nCentBins \'" << nCentBins << "\' is greater than nMaxCentBins \'" << nMaxCentBins << "\'. return 1" << std::endl;
    return 1;
  }

  std::cout << std::endl;
  const Int_t nMaxJtAbsEtaBins = 6;
  const Int_t nJtAbsEtaBins = nJtAbsEtaBinsTemp;
  if(nJtAbsEtaBins > nMaxJtAbsEtaBins){
    std::cout << "nJtAbsEtaBins \'" << nJtAbsEtaBins << "\' is greater than nMaxJtAbsEtaBins \'" << nMaxJtAbsEtaBins << "\'. return 1" << std::endl;
    return 1;
  }

  Double_t jtAbsEtaBinsLow[nMaxJtAbsEtaBins];
  Double_t jtAbsEtaBinsHi[nMaxJtAbsEtaBins];
  std::cout << "nJtAbsEtaBins: ";

  for(Int_t jI = 0; jI < nJtAbsEtaBins; ++jI){
    jtAbsEtaBinsLow[jI] = jtAbsEtaBinsLowTemp[jI];
    jtAbsEtaBinsHi[jI] = jtAbsEtaBinsHiTemp[jI];
    std::cout << " " << jtAbsEtaBinsLow[jI] << "-" << jtAbsEtaBinsHi[jI] << ",";
  }
  std::cout << std::endl;
                               
  const Int_t nID = nIDTemp;     

  std::cout << "nCentBins: " << nCentBins << std::endl;
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    std::cout << " " << cI << "/" << nCentBins << ": " << centBinsLow[cI] << "-" << centBinsHi[cI] << std::endl;
  }

  const Int_t nResponseMod = cutProp.GetNResponseMod();
  std::vector<double> responseMod = cutProp.GetResponseMod();
  std::vector<double> jerVarData = cutProp.GetJERVarData();

  const Int_t nSyst = cutProp.GetNSyst();
  std::vector<std::string> systStr = cutProp.GetSystStr();

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
  const Double_t rightMargin[nPad] = {0.005, 0.001, 0.001};
  const Double_t topMargin[nPad] = {leftMargin[0], leftMargin[0]*(padYHi[0] - padYLow[0])/(padYHi[1] - padYLow[1]), 0.001};
  const Double_t bottomMargin[nPad] = {leftMargin[0], 0.001, leftMargin[0]*(padYHi[0] - padYLow[0])/(padYHi[2] - padYLow[2])};

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  checkMakeDir("pdfDir");
  checkMakeDir("pdfDir/" + dateStr);
  checkMakeDir("pdfDir/" + dateStr + "/Validate");

  for(Int_t jI = 0; jI < nJets; ++jI){
    std::string dirName = jetDirList[jI];
    dirName = dirName.substr(0, dirName.find("/"));

    for(Int_t cI = 0; cI < nCentBins; ++cI){
      std::string centStr = "PP";
      std::string centStr2 = "PP";
      if(!isResponsePP){
	centStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);
	centStr2 = std::to_string(centBinsLow[cI]) + "-" + std::to_string(centBinsHi[cI]) + "%";
      }

      for(Int_t idI = 0; idI < nID; ++idI){
	for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
	  const std::string jtAbsEtaStr = "AbsEta" + prettyString(jtAbsEtaBinsLow[aI], 1, true) + "to" + prettyString(jtAbsEtaBinsHi[aI], 1, true);
	  const std::string jtAbsEtaStr2 = prettyString(jtAbsEtaBinsLow[aI], 1, false) + "<|#eta|<" + prettyString(jtAbsEtaBinsHi[aI], 1, false);

	  for(Int_t mI = 0; mI < nResponseMod; ++mI){
	    const std::string resStr = "ResponseMod" + prettyString(responseMod[mI], 2, true);
	    const std::string resStr2 = "#sigma_{#frac{Data}{MC}}=" + prettyString(responseMod[mI], 2, false);
	   
	    for(Int_t sI = 0; sI < nSyst; ++sI){
	      std::string tempSystStr = "_" + systStr[sI] + "_";
	      while(tempSystStr.find("__") != std::string::npos){tempSystStr.replace(tempSystStr.find("__"), 2, "_");}

	      RooUnfoldResponse* rooRes_p = (RooUnfoldResponse*)responseFile_p->Get((dirName + "/rooResponse_" + dirName + "_" + centStr + "_" + idStr[idI] + "_" + resStr + "_" + jtAbsEtaStr + tempSystStr + "RecoGenAsymm_h").c_str());

	      TH1D* recoJtPt_RecoGenAsymm_h = (TH1D*)responseFile_p->Get((dirName + "/recoJtPt_" + dirName + "_" + centStr + "_" + idStr[idI] + "_" + resStr + "_" + jtAbsEtaStr + tempSystStr + "GoodGen" + paraString + "_h").c_str());

	      TH1D* genJtPt_h = (TH1D*)responseFile_p->Get((dirName + "/genJtPt_" + dirName + "_" + centStr + "_" + idStr[idI] + "_" + resStr + "_" + jtAbsEtaStr + tempSystStr+ "GoodReco" + paraString + "_h").c_str());
	     
	      slideTitles->push_back(rooRes_p->GetName());
	      /*
	      std::cout << "Unfold list (" << tempSystStr << "): " << std::endl;
	      std::cout << " " << rooRes_p->GetName() << std::endl;	
	      std::cout << " " << recoJtPt_RecoGenAsymm_h->GetName() << std::endl;
	      std::cout << " " << genJtPt_h->GetName() << std::endl;
      
	      std::cout << " Reco Print: " << std::endl;
	      std::cout << "  " << recoJtPt_RecoGenAsymm_h->GetName() << std::endl;
	      recoJtPt_RecoGenAsymm_h->Print("ALL");
	      std::cout << "  " << rooRes_p->Hmeasured()->GetName() << std::endl;
	      rooRes_p->Hmeasured()->Print("ALL");
	      std::cout << "  " << rooRes_p->Hfakes()->GetName() << std::endl;
	      rooRes_p->Hfakes()->Print("ALL");

	      std::cout << " Gen Print: " << std::endl;
	      std::cout << "  " << genJtPt_h->GetName() << std::endl;
	      genJtPt_h->Print("ALL");
	      std::cout << "  " << rooRes_p->Htruth()->GetName() << std::endl;
	      rooRes_p->Htruth()->Print("ALL");
	      */

	      recoJtPt_RecoGenAsymm_h->GetXaxis()->SetTitleFont(43);
	      recoJtPt_RecoGenAsymm_h->GetYaxis()->SetTitleFont(43);
	      recoJtPt_RecoGenAsymm_h->GetXaxis()->SetLabelFont(43);
	      recoJtPt_RecoGenAsymm_h->GetYaxis()->SetLabelFont(43);

	      recoJtPt_RecoGenAsymm_h->GetXaxis()->SetTitleSize(14);
	      recoJtPt_RecoGenAsymm_h->GetYaxis()->SetTitleSize(14);
	      recoJtPt_RecoGenAsymm_h->GetXaxis()->SetLabelSize(14);
	      recoJtPt_RecoGenAsymm_h->GetYaxis()->SetLabelSize(14);

	      genJtPt_h->GetXaxis()->SetTitleFont(43);
	      genJtPt_h->GetYaxis()->SetTitleFont(43);
	      genJtPt_h->GetXaxis()->SetLabelFont(43);
	      genJtPt_h->GetYaxis()->SetLabelFont(43);

	      genJtPt_h->GetXaxis()->SetTitleSize(14);
	      genJtPt_h->GetYaxis()->SetTitleSize(14);
	      genJtPt_h->GetXaxis()->SetLabelSize(14);
	      genJtPt_h->GetYaxis()->SetLabelSize(14);
	    
	      pdfPerSlide->push_back({});

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
		}

		canv_p->cd();
		pads_p[0]->cd();
		
		recoJtPt_RecoGenAsymm_h->SetMarkerColor(colors[0]);
		recoJtPt_RecoGenAsymm_h->SetLineColor(colors[0]);
		recoJtPt_RecoGenAsymm_h->SetMarkerStyle(styles[0]);
		recoJtPt_RecoGenAsymm_h->SetMarkerSize(0.8);
		
		genJtPt_h->SetMarkerColor(colors[1]);
		genJtPt_h->SetLineColor(colors[1]);
		genJtPt_h->SetMarkerStyle(styles[1]);
		genJtPt_h->SetMarkerSize(0.8);
		
		recoJtPt_RecoGenAsymm_h->GetXaxis()->SetNdivisions(505);
		recoJtPt_RecoGenAsymm_h->GetYaxis()->SetNdivisions(505);
		
		Double_t maxVal = 0.0;
		Double_t minVal = 999999.;
		
		for(Int_t bI = 0; bI < recoJtPt_RecoGenAsymm_h->GetNbinsX(); ++bI){
		  if(recoJtPt_RecoGenAsymm_h->GetBinContent(bI+1) > maxVal) maxVal = recoJtPt_RecoGenAsymm_h->GetBinContent(bI+1);
		  if(genJtPt_h->GetBinContent(bI+1) > maxVal) maxVal = genJtPt_h->GetBinContent(bI+1);
	      
		  if(recoJtPt_RecoGenAsymm_h->GetBinContent(bI+1) < minVal && recoJtPt_RecoGenAsymm_h->GetBinContent(bI+1) > 0) minVal = recoJtPt_RecoGenAsymm_h->GetBinContent(bI+1);
		  if(genJtPt_h->GetBinContent(bI+1) < minVal && genJtPt_h->GetBinContent(bI+1) > 0) minVal = genJtPt_h->GetBinContent(bI+1);
		}
	  
		maxVal *= 5.;
		minVal /= 5.;

		Int_t nBinsForDraw = 100;
		Double_t binsForDraw[nBinsForDraw+1];
		getLogBins(minVal, maxVal, nBinsForDraw, binsForDraw);
		
		Double_t xLowReco = recoJtPt_RecoGenAsymm_h->GetBinLowEdge(1);
		//		Double_t xLowGen = genJtPt_h->GetBinLowEdge(1);
		Double_t drawHigh = maxVal*maxVal/binsForDraw[90];
		Double_t drawHigh2 = maxVal*maxVal/binsForDraw[95];
		Double_t drawLow = minVal*minVal/binsForDraw[4];
		Double_t drawLowTight = minVal*minVal/binsForDraw[1];
		Double_t drawLowLoose = minVal*minVal/binsForDraw[6];

		recoJtPt_RecoGenAsymm_h->SetMaximum(maxVal);
		recoJtPt_RecoGenAsymm_h->SetMinimum(minVal);

		bool doLogX = false;
		if(recoJtPt_RecoGenAsymm_h->GetBinWidth(1)*3 < recoJtPt_RecoGenAsymm_h->GetBinWidth(recoJtPt_RecoGenAsymm_h->GetNbinsX()-1)) doLogX = true;
		
		recoJtPt_RecoGenAsymm_h->DrawCopy("HIST E1 P");
		genJtPt_h->DrawCopy("HIST E1 P SAME");

		TLatex* label_p = new TLatex();
		label_p->SetTextFont(43);
		label_p->SetTextSize(14);
	      
		gPad->SetLogy();
		if(doLogX){
		  gPad->SetLogx();

		  drawWhiteBox(900, 1200, drawLowLoose, drawLowTight);

		  label_p->DrawLatex(300, drawLow, "300");
		  label_p->DrawLatex(400, drawLow, "400");
		  label_p->DrawLatex(600, drawLow, "600");
		  label_p->DrawLatex(800, drawLow, "800");
		  label_p->DrawLatex(1000, drawLow, "1000");
		}


		std::string labelStr = dirName + ", " + centStr2 + ", " + idStr[idI] + ", " + jtAbsEtaStr2;
		std::string labelStr2 = resStr2 + ", " + tempSystStr + ", Bayes=" + std::to_string(bI);
		label_p->DrawLatex(xLowReco, drawHigh, labelStr.c_str());
		label_p->DrawLatex(xLowReco, drawHigh2, labelStr2.c_str());

		TLegend* leg_p = new TLegend(0.2, 0.2, 0.5, 0.5);
		leg_p->SetTextFont(43);
		leg_p->SetTextSize(14);
		leg_p->SetBorderSize(0);
		leg_p->SetFillColor(0);
		leg_p->SetFillStyle(0);

		leg_p->AddEntry(recoJtPt_RecoGenAsymm_h, "Reco.", "P L");
		leg_p->AddEntry(genJtPt_h, "Gen.", "P L");

		leg_p->Draw("SAME");
		
		canv_p->cd();
		pads_p[1]->cd();
		
		RooUnfoldBayes bayes(rooRes_p, recoJtPt_RecoGenAsymm_h, 1+bI, false, "name");
		bayes.SetVerbose(0);
		
		TH1D* unfold_h = (TH1D*)bayes.Hreco(RooUnfold::kCovToy);
		
		unfold_h->SetMarkerColor(colors[2]);
		unfold_h->SetLineColor(colors[2]);
		unfold_h->SetMarkerStyle(styles[2]);
		unfold_h->SetMarkerSize(0.8);

		unfold_h->GetXaxis()->SetNdivisions(505);
		unfold_h->GetYaxis()->SetNdivisions(505);
		centerTitles(unfold_h);
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

		TLegend* leg2_p = new TLegend(0.2, 0.2, 0.5, 0.5);
		leg2_p->SetTextFont(43);
		leg2_p->SetTextSize(14);
		leg2_p->SetBorderSize(0);
		leg2_p->SetFillColor(0);
		leg2_p->SetFillStyle(0);

		leg2_p->AddEntry(unfold_h, "Unfolded", "P L");
		leg2_p->AddEntry(genJtPt_h, "Gen.", "P L");

		leg2_p->Draw("SAME");
		
		gPad->SetLogy();
		if(doLogX) gPad->SetLogx();

		canv_p->cd();
		pads_p[2]->cd();

		minVal = 0.75;
		maxVal = 1.25;
		
		if(!doParaFills){
		  minVal = 0.98;
		  maxVal = 1.02;
		}
		
		unfold_h->Divide(genJtPt_h);
		unfold_h->SetMaximum(maxVal);
		unfold_h->SetMinimum(minVal);
		unfold_h->DrawCopy("HIST E1 P");

		getLinBins(minVal, maxVal, nBinsForDraw, binsForDraw);
		drawLow = minVal*minVal/binsForDraw[10];
		drawLowTight = minVal*minVal/binsForDraw[1];
		drawLowLoose = minVal*minVal/binsForDraw[20];

		if(doLogX){
		  gPad->SetLogx();

		  drawWhiteBox(900, 1200, drawLowLoose, drawLowTight);

		  label_p->DrawLatex(300, drawLow, "300");
		  label_p->DrawLatex(400, drawLow, "400");
		  label_p->DrawLatex(600, drawLow, "600");
		  label_p->DrawLatex(800, drawLow, "800");
		  label_p->DrawLatex(1000, drawLow, "1000");
		}
		
		const std::string saveName = dirName + "_" + centStr + "_" + idStr[idI] + "_" + resStr + "_" + jtAbsEtaStr + tempSystStr + "Bayes" + std::to_string(bI+1) + paraString + "_" + dateStr + ".pdf";

		pdfPerSlide->at(pdfPerSlide->size()-1).push_back(saveName);
		quietSaveAs(canv_p, "pdfDir/" + dateStr + "/Validate/" + saveName);

		delete leg_p;
		delete leg2_p;
		delete label_p;
		for(Int_t pI = 0; pI < nPad; ++pI){
		  delete pads_p[pI];
		}	    
		delete canv_p;
	      }
	    }
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

  texSlideCreator tex;
  tex.Clean();
  tex.Init(argv[1]);
  tex.SetAuthor("Chris McGinn");
  std::vector<std::string> slideTitles;
  std::vector<std::vector<std::string> > pdfPerSlide;

  int retVal = 0;
  retVal += validateJetResponse(argv[1], false, &slideTitles, &pdfPerSlide);
  retVal += validateJetResponse(argv[1], true, &slideTitles, &pdfPerSlide);

  tex.SetSlideTitles(slideTitles);
  tex.SetSlidePdfs(pdfPerSlide);
  if(!(tex.CreateTexSlides())){
    std::cout << "Warning: .tex slide creation failed" << std::endl;
  }

  std::cout << "Job complete. Return " << retVal << "." << std::endl;
  return retVal;
}
