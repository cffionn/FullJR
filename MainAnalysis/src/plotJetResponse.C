//cpp dependencies
#include <iostream>
#include <string>
#include <vector>

//ROOT dependencies
#include "TCanvas.h"
#include "TDatime.h"
#include "TFile.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TPad.h"
#include "TStyle.h"

//Local dependencies
#include "MainAnalysis/include/cutPropagator.h"
#include "MainAnalysis/include/plotToSuperPlotDim.h"
#include "MainAnalysis/include/smallOrLargeR.h"
#include "MainAnalysis/include/texSlideCreator.h"

//Non-local FullJR (Utility, etc.) dependencies
#include "Utility/include/checkMakeDir.h"
#include "Utility/include/getLinBins.h"
#include "Utility/include/getLogBins.h"
#include "Utility/include/histDefUtility.h"
#include "Utility/include/plotUtilities.h"
#include "Utility/include/returnRootFileContentsList.h"
#include "Utility/include/vanGoghPalette.h"

int plotJetResponse(const std::string inResponseName)
{
  std::vector<std::string> slideTitles;
  std::vector<std::vector<std::string > > pdfPerSlide;

  TFile* responseFile_p = new TFile(inResponseName.c_str(), "READ");
  std::vector<std::string> jetDirList = returnRootFileContentsList(responseFile_p, "TDirectoryFile", "JetAnalyzer");

  std::cout << "Validating " << jetDirList.size() << " jets..." << std::endl;
  for(unsigned int jI = 0; jI < jetDirList.size(); ++jI){
    std::cout << " " << jI << "/" << jetDirList.size() << ": " << jetDirList.at(jI) << std::endl;
  }

  cutPropagator cutProp;
  cutProp.Clean();
  cutProp.GetAllVarFromFile(responseFile_p); 

  Bool_t isPP = cutProp.GetIsPP();

  Int_t nCentBins = cutProp.GetNCentBins();
  std::vector<Int_t> centBinsLow = cutProp.GetCentBinsLow();
  std::vector<Int_t> centBinsHi = cutProp.GetCentBinsHi();

  Int_t nGenJtPtBinsSmallRTemp = cutProp.GetNGenJtPtBinsSmallR();
  std::vector<Double_t> genJtPtBinsSmallRTemp = cutProp.GetGenJtPtBinsSmallR();

  Int_t nGenJtPtBinsLargeRTemp = cutProp.GetNGenJtPtBinsLargeR();
  std::vector<Double_t> genJtPtBinsLargeRTemp = cutProp.GetGenJtPtBinsLargeR();

  Int_t nJtAbsEtaBinsTemp = cutProp.GetNJtAbsEtaBins();
  std::vector<Double_t> jtAbsEtaBinsLowTemp = cutProp.GetJtAbsEtaBinsLow();
  std::vector<Double_t> jtAbsEtaBinsHiTemp = cutProp.GetJtAbsEtaBinsHi();

  Int_t nIDTemp = cutProp.GetNID();
  std::vector<std::string> idStrTemp = cutProp.GetIdStr();

  const Int_t nMaxResponseMod = 4;
  const Int_t nResponseMod = cutProp.GetNResponseMod();

  if(nResponseMod > nMaxResponseMod){
    std::cout << "nResponseMod \'" << nResponseMod << "\' is greater than nMaxResponseMod \'" << nMaxResponseMod << "\'. return 1" << std::endl;
    return 1;
  }

  std::vector<double> responseMod = cutProp.GetResponseMod();

  if(nCentBins < 0) std::cout << "nCentBins less than 0. please check input file. return 1" << std::endl;
  if(nGenJtPtBinsSmallRTemp < 0) std::cout << "nGenJtPtBinsSmallRTemp less than 0. please check input file. return 1" << std::endl;
  if(nGenJtPtBinsLargeRTemp < 0) std::cout << "nGenJtPtBinsLargeRTemp less than 0. please check input file. return 1" << std::endl;
  if(nJtAbsEtaBinsTemp < 0) std::cout << "nJtAbsEtaBinsTemp less than 0. please check input file. return 1" << std::endl;
  if(nIDTemp < 0) std::cout << "nIDTemp less than 0. please check input file. return 1" << std::endl;
  
  if(nCentBins < 0 || nGenJtPtBinsSmallRTemp < 0 || nGenJtPtBinsLargeRTemp < 0 || nJtAbsEtaBinsTemp < 0 || nIDTemp < 0){
    responseFile_p->Close();
    delete responseFile_p;
    return 1;
  }

  smallOrLargeR rReader;
  if(!rReader.CheckNGenJtPtBinsSmallR(nGenJtPtBinsSmallRTemp)) return 1;
  if(!rReader.CheckNGenJtPtBinsLargeR(nGenJtPtBinsLargeRTemp)) return 1;
  if(!rReader.CheckGenJtPtBinsSmallR(genJtPtBinsSmallRTemp)) return 1;
  if(!rReader.CheckGenJtPtBinsLargeR(genJtPtBinsLargeRTemp)) return 1;

  const Int_t nMaxJtAbsEtaBins = 6;
  const Int_t nJtAbsEtaBins = nJtAbsEtaBinsTemp;

  if(nJtAbsEtaBins > nMaxJtAbsEtaBins){
    std::cout << "nJtAbsEtaBins \'" << nJtAbsEtaBins << "\' is greater than nMaxJtAbsEtaBins \'" << nMaxJtAbsEtaBins << "\'. return 1" << std::endl;
    return 1;
  }

  Double_t jtAbsEtaBinsLow[nJtAbsEtaBins];
  Double_t jtAbsEtaBinsHi[nJtAbsEtaBins];
  for(Int_t jI = 0; jI < nJtAbsEtaBins; ++jI){
    jtAbsEtaBinsLow[jI] = jtAbsEtaBinsLowTemp.at(jI);
    jtAbsEtaBinsHi[jI] = jtAbsEtaBinsHiTemp.at(jI);
  }

  std::cout << "nCentBins: " << nCentBins << std::endl;
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    std::cout << " " << cI << "/" << nCentBins << ": " << centBinsLow.at(cI) << "-" << centBinsHi.at(cI) << std::endl;
  }

  const Int_t nMaxID = 6;
  const Int_t nID = nIDTemp;
  std::string idStr[nID];
  for(Int_t i = 0; i < nID; ++i){
    idStr[i] = idStrTemp.at(i);
  }

  if(nID > nMaxID){
    std::cout << "nID \'" << nID << "\' is greater than nMaxID \'" << nMaxID << "\'. return 1" << std::endl;
    return 1;
  }

  
  const Int_t nManip = 6;
  const std::string recoTruncStr[nManip] = {"_RecoGenSymm", "_RecoGenAsymm", "_RecoGenSymm", "_RecoGenAsymm", "_RecoGenSymm", "_RecoGenAsymm"};
  Bool_t renormX[nManip] = {false, false, true, true, false, false};
  Bool_t renormY[nManip] = {false, false, false, false, true, true};

  const Int_t nMaxJets = 10;
  const Int_t nJets = jetDirList.size(); 

  if(nJets > nMaxJets){
    std::cout << "nJets \'" << nJets << "\' is greater than nMaxJets \'" << nMaxJets << "\'. return 1" << std::endl;
    return 1;
  }

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  checkMakeDir("pdfDir");
  checkMakeDir("pdfDir/" + dateStr);


  for(Int_t jI = 0; jI < nJets; ++jI){
    std::string dirName = jetDirList.at(jI);
    dirName = dirName.substr(0, dirName.find("/"));

    const std::string jetName = dirName.substr(0, dirName.find("JetAna"));

    const Int_t rVal = getRVal(dirName);
    const bool isSmallR = rReader.GetIsSmallR(rVal);
    const Int_t nMaxPtBins = 50;
    const Int_t nGenJtPtBins = rReader.GetSmallOrLargeRNBins(isSmallR, true);
    Double_t genJtPtBins[nMaxPtBins+1];
    rReader.GetSmallOrLargeRBins(isSmallR, true, nGenJtPtBins+1, genJtPtBins);

    std::cout << "Bins for rVal = " << rVal << std::endl;
    std::cout << " ";
    for(Int_t gI = 0; gI < nGenJtPtBins+1; ++gI){
      std::cout << genJtPtBins[gI] << ", ";
    }
    std::cout << std::endl;
  
    for(Int_t cI = 0; cI < nCentBins; ++cI){
      std::string centStr = "PP";
      std::string centStr2 = "PP";
      if(!isPP){
	centStr = "Cent" + std::to_string(centBinsLow.at(cI)) + "to" + std::to_string(centBinsHi.at(cI));
	centStr2 = std::to_string(centBinsLow.at(cI)) + "-" + std::to_string(centBinsHi.at(cI)) + "%";
      }

      for(Int_t iI = 0; iI < nID; ++iI){	

	for(Int_t modI = 0; modI < nResponseMod; ++modI){
	  std::string resStr = "ResponseMod" + prettyString(responseMod[modI], 2, true);
	  std::string resStr2 = "Res.x" + prettyString(responseMod[modI], 2, false);

	  for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
	    const std::string jtAbsEtaStr = "AbsEta" + prettyString(jtAbsEtaBinsLow[aI], 1, true) + "to" + prettyString(jtAbsEtaBinsHi[aI], 1, true);
	    const std::string jtAbsEtaStr2 = prettyString(jtAbsEtaBinsLow[aI], 1, false) + "<|\\eta|<" + prettyString(jtAbsEtaBinsHi[aI], 1, false);

	    bool doSlides = isStrSame(jtAbsEtaStr, "AbsEta0p0to2p0");
	    doSlides = doSlides && isStrSame(idStr[iI], "LightMUAndCHID");

	    if(doSlides){
	      const std::string slideTitle = jetName + ", " + centStr2 + ", " + idStr[iI] + ", " + resStr2 + ", $" + jtAbsEtaStr2 + "$";
	      slideTitles.push_back(slideTitle);
	      pdfPerSlide.push_back({});
	    }

	    for(Int_t mI = 0; mI < nManip; ++mI){
	      TH2D* response_h = (TH2D*)(responseFile_p->Get((dirName + "/response_" + dirName + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_" + jtAbsEtaStr + recoTruncStr[mI] + "_h").c_str())->Clone("response_h"));
	      
	      bool doPrint = true;
	      if(jtAbsEtaStr.find("AbsEta0p0to2p0") == std::string::npos) doPrint = false;
	      else if(resStr.find("ResponseMod0p00") == std::string::npos) doPrint = false;
	      else if(idStr[iI].find("LightMUAndCHID") == std::string::npos) doPrint = false;
	      else if(centStr.find("Cent0to10") == std::string::npos) doPrint = false;
	     
	      if(doPrint){
		std::cout << "PRINTING response \'" << response_h->GetName() << "\'." << std::endl;
		response_h->Print("ALL");
		std::cout << "END PRINT" << std::endl;
	      }
	      
	      bool doLogX = false;
	      if(response_h->GetXaxis()->GetBinWidth(1)*3 < response_h->GetXaxis()->GetBinWidth(response_h->GetNbinsX()-1)) doLogX = true;
	      
	      TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);
	      canv_p->SetTopMargin(0.01);
	      canv_p->SetBottomMargin(0.14);
	      canv_p->SetLeftMargin(0.14);
	      canv_p->SetRightMargin(0.14);
	      
	      gStyle->SetOptStat(0);
	      
	      if(renormX[mI]){
		for(Int_t bIX = 0; bIX < response_h->GetNbinsX(); ++bIX){
		  Double_t xSum = 0.;
		  
		  for(Int_t bIY = 0; bIY < response_h->GetNbinsY(); ++bIY){
		    xSum += response_h->GetBinContent(bIX+1, bIY+1);
		  }
		  
		  if(xSum <= 0) continue;
		  
		  for(Int_t bIY = 0; bIY < response_h->GetNbinsY(); ++bIY){
		    response_h->SetBinContent(bIX+1, bIY+1, response_h->GetBinContent(bIX+1, bIY+1)/xSum);
		    response_h->SetBinError(bIX+1, bIY+1, response_h->GetBinError(bIX+1, bIY+1)/xSum);
		  }
		}	
	      }
	      else if(renormY[mI]){
		for(Int_t bIY = 0; bIY < response_h->GetNbinsY(); ++bIY){
		  Double_t ySum = 0.;
		  
		  for(Int_t bIX = 0; bIX < response_h->GetNbinsX(); ++bIX){
		    ySum += response_h->GetBinContent(bIX+1, bIY+1);
		  }
		  
		  if(ySum <= 0) continue;
		  
		  for(Int_t bIX = 0; bIX < response_h->GetNbinsX(); ++bIX){
		    response_h->SetBinContent(bIX+1, bIY+1, response_h->GetBinContent(bIX+1, bIY+1)/ySum);
		    response_h->SetBinError(bIX+1, bIY+1, response_h->GetBinError(bIX+1, bIY+1)/ySum);
		  }
		}	
	      }	
	      
	      response_h->DrawCopy("COLZ TEXT");
	      gStyle->SetPaintTextFormat("1.3f");
	      gPad->SetLogz();
	      if(doLogX){
		gPad->SetLogx();
		gPad->SetLogy();
	      }
	      
	      std::string renormStr = "nonorm";
	      if(renormX[mI]) renormStr = "renormX";
	      else if(renormY[mI]) renormStr = "renormY";

	      Double_t histMinX = response_h->GetXaxis()->GetBinLowEdge(1);
	      Double_t histMaxX = response_h->GetXaxis()->GetBinLowEdge(response_h->GetNbinsX());
	      Double_t histMinY = response_h->GetYaxis()->GetBinLowEdge(1);
	      Double_t histMaxY = response_h->GetYaxis()->GetBinLowEdge(response_h->GetNbinsY());

	      delete response_h;

	      const Int_t nBins = 40;
	      Double_t binsLinX[nBins+1];
	      Double_t binsLogX[nBins+1];
	      Double_t binsLinY[nBins+1];
	      Double_t binsLogY[nBins+1];

	      getLinBins(histMinX, histMaxX, nBins, binsLinX);
	      getLogBins(histMinX, histMaxX, nBins, binsLogX);
	      getLinBins(histMinY, histMaxY, nBins, binsLinY);
	      getLogBins(histMinY, histMaxY, nBins, binsLogY);	     

	      TLatex* label_p = new TLatex();
	      label_p->SetTextFont(43);
	      label_p->SetTextSize(12);

	      if(doLogX){
		label_p->DrawLatex(binsLogX[1], binsLogY[nBins-1], renormStr.c_str());
		label_p->DrawLatex(binsLogX[1], binsLogY[nBins-2], recoTruncStr[mI].substr(1, recoTruncStr[mI].size()).c_str());

		label_p->DrawLatex(binsLogX[0]*9./10., 300, "300");
		label_p->DrawLatex(binsLogX[0]*9./10., 400, "400");
		label_p->DrawLatex(binsLogX[0]*9./10., 800, "800");

		label_p->DrawLatex(300, binsLogY[0]*9./10., "300");
		label_p->DrawLatex(400, binsLogY[0]*9./10., "400");
		label_p->DrawLatex(800, binsLogY[0]*9./10., "800");
	      }
	      else{
		label_p->DrawLatex(binsLinX[1], binsLinY[19], renormStr.c_str());
		label_p->DrawLatex(binsLinX[1], binsLinY[18], recoTruncStr[mI].c_str());
	      }

	     
	      delete label_p;


	      const std::string saveName = "response_" + dirName + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_" + jtAbsEtaStr + recoTruncStr[mI] + "_" + renormStr + "_" + dateStr  + ".pdf";
	      if(doSlides) pdfPerSlide.at(pdfPerSlide.size()-1).push_back(saveName);
	      const std::string finalSaveName = "pdfDir/" + dateStr + "/" + saveName;
	      quietSaveAs(canv_p, finalSaveName);
	      delete canv_p;
	    }

	    
	    plotToSuperPlotDim superDim;
	    Int_t nPadX = superDim.GetNPlotsX(nGenJtPtBins);
	    Int_t nPadY = superDim.GetNPlotsY(nGenJtPtBins);
	  	    
	    TCanvas* recoJtPtPerGenPtBin_p = new TCanvas("recoJtPtPerGenPtBin_p", "", 450*nPadX, 450*nPadY);
	    recoJtPtPerGenPtBin_p->SetTopMargin(0.01);
	    recoJtPtPerGenPtBin_p->SetBottomMargin(0.01);
	    recoJtPtPerGenPtBin_p->SetLeftMargin(0.01);
	    recoJtPtPerGenPtBin_p->SetRightMargin(0.01);
	    
	    recoJtPtPerGenPtBin_p->Divide(nPadX, nPadY);
	    
	    for(Int_t jI = 0; jI < nGenJtPtBins; ++jI){
	      const std::string jtPtStr = "Pt" + prettyString(genJtPtBins[jI], 1, true) + "to" + prettyString(genJtPtBins[jI+1], 1, true);
	      
	      recoJtPtPerGenPtBin_p->cd();
	      recoJtPtPerGenPtBin_p->cd(jI+1);
	      
	      TH1D* temp_p = (TH1D*)responseFile_p->Get((dirName + "/recoJtPtPerGenPtBin_" + dirName + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_Gen" + jtPtStr + "_" + jtAbsEtaStr + "_h").c_str());
	      
	      temp_p->DrawCopy("HIST E1 P");
	    }
	    
	    const std::string recoJtSaveName = "pdfDir/" + dateStr + "/recoJtPtPerGenPtBin_" + dirName + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_" + jtAbsEtaStr + "_" + dateStr + ".pdf";
	    quietSaveAs(recoJtPtPerGenPtBin_p, recoJtSaveName);
	    delete recoJtPtPerGenPtBin_p;
	    
	    TCanvas* recoJtPtPerGenPtBinWeighted_p = new TCanvas("recoJtPtPerGenPtBinWeighted_p", "", 450*nPadX, 450*nPadY);
	    recoJtPtPerGenPtBinWeighted_p->SetTopMargin(0.01);
	    recoJtPtPerGenPtBinWeighted_p->SetBottomMargin(0.01);
	    recoJtPtPerGenPtBinWeighted_p->SetLeftMargin(0.01);
	    recoJtPtPerGenPtBinWeighted_p->SetRightMargin(0.01);
	    
	    recoJtPtPerGenPtBinWeighted_p->Divide(nPadX, nPadY);
	    
	    for(Int_t jI = 0; jI < nGenJtPtBins; ++jI){
	      const std::string jtPtStr = "Pt" + prettyString(genJtPtBins[jI], 1, true) + "to" + prettyString(genJtPtBins[jI+1], 1, true);
	    
	      recoJtPtPerGenPtBinWeighted_p->cd();
	      recoJtPtPerGenPtBinWeighted_p->cd(jI+1);
	      
	      
	      TH1D* temp_p = (TH1D*)responseFile_p->Get((dirName + "/recoJtPtPerGenPtBin_" + dirName + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_Gen" + jtPtStr + "_" + jtAbsEtaStr + "_Weighted_h").c_str());
	      
	      temp_p->DrawCopy("HIST E1 P");
	    }
	    
	    const std::string tempSaveName = "pdfDir/" + dateStr + "/recoJtPtPerGenPtBin_Weighted_" + dirName + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_" + jtAbsEtaStr + "_" + dateStr + ".pdf";
	    quietSaveAs(recoJtPtPerGenPtBinWeighted_p, tempSaveName);
	    delete recoJtPtPerGenPtBinWeighted_p;
	  }
	}
      }
    }
  }
  
  responseFile_p->Close();
  delete responseFile_p;

  //produce some latex slides
  texSlideCreator tex;
  tex.Clean();
  tex.Init(inResponseName);

  if(nJets == 1){
    std::string jetTag = jetDirList.at(0);
    jetTag = jetTag.substr(0, jetTag.find("/"));

    tex.InitTag("AllPlots_" + jetTag);
  }
  else tex.InitTag("AllPlots");

  tex.SetAuthor("Christopher McGinn");
  tex.SetSlideTitles(slideTitles);
  tex.SetSlidePdfs(pdfPerSlide);
  if(!(tex.CreateTexSlides())){
    std::cout << "Warning: .tex slide creation failed" << std::endl;
  }

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/plotJetResponse.exe <inResponseName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += plotJetResponse(argv[1]);

  std::cout << "Job complete. Return " << retVal << "." << std::endl;
  return retVal;
}
