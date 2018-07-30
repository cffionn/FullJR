//cpp dependencies
#include <iostream>
#include <string>

//ROOT dependencies
#include "TFile.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TDatime.h"
#include "TStyle.h"

//Local dependencies                                                                                                                             
#include "MainAnalysis/include/cutPropagator.h"

//Non-local FullJR (Utility, etc.) dependencies                                                                                                      
#include "Utility/include/checkMakeDir.h"
#include "Utility/include/returnRootFileContentsList.h"
#include "Utility/include/histDefUtility.h"
#include "Utility/include/plotUtilities.h"
#include "Utility/include/vanGoghPalette.h"


void getPadsXY(const Int_t nPlots, Int_t *nPadX, Int_t *nPadY)
{
  if(nPlots == 1){
    (*nPadX) = 1;
    (*nPadY) = 1;
  }
  else if(nPlots == 2){
    (*nPadX) = 2;
    (*nPadY) = 1;
  }
  else if(nPlots == 3){
    (*nPadX) = 2;
    (*nPadY) = 2;
  }
  else if(nPlots == 4){
    (*nPadX) = 2;
    (*nPadY) = 2;
  }
  else if(nPlots == 5){
    (*nPadX) = 3;
    (*nPadY) = 2;
  }
  else if(nPlots == 6){
    (*nPadX) = 3;
    (*nPadY) = 2;
  }
  else if(nPlots == 7){
    (*nPadX) = 4;
    (*nPadY) = 2;
  }
  else if(nPlots == 8){
    (*nPadX) = 4;
    (*nPadY) = 2;
  }
  else if(nPlots == 9){
    (*nPadX) = 3;
    (*nPadY) = 3;
  }
  else if(nPlots == 10){
    (*nPadX) = 4;
    (*nPadY) = 3;
  }
  else if(nPlots == 11){
    (*nPadX) = 4;
    (*nPadY) = 3;
  }
  else if(nPlots == 12){
    (*nPadX) = 4;
    (*nPadY) = 3;
  }
  else if(nPlots == 13){
    (*nPadX) = 4;
    (*nPadY) = 4;
  }
  else if(nPlots == 14){
    (*nPadX) = 4;
    (*nPadY) = 4;
  }
  else if(nPlots == 15){
    (*nPadX) = 4;
    (*nPadY) = 4;
  }
  else if(nPlots == 16){
    (*nPadX) = 4;
    (*nPadY) = 4;
  }
  else{
    std::cout << "WARNING: nPlots \'" << nPlots << "\' has no specified value. return" << std::endl;
  }

  return;
}


int plotJetResponse(const std::string inResponseName)
{
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

  Int_t nJtPtBinsTemp = cutProp.GetNJtPtBins();
  std::vector<Double_t> jtPtBinsTemp = cutProp.GetJtPtBins();

  Int_t nJtAbsEtaBinsTemp = cutProp.GetNJtAbsEtaBins();
  std::vector<Double_t> jtAbsEtaBinsLowTemp = cutProp.GetJtAbsEtaBinsLow();
  std::vector<Double_t> jtAbsEtaBinsHiTemp = cutProp.GetJtAbsEtaBinsHi();

  Int_t nIDTemp = cutProp.GetNID();
  std::vector<std::string> idStrTemp = cutProp.GetIdStr();

  const Int_t nResponseMod = cutProp.GetNResponseMod();
  std::vector<double> responseMod = cutProp.GetResponseMod();

  if(nCentBins < 0) std::cout << "nCentBins less than 0. please check input file. return 1" << std::endl;
  if(nJtPtBinsTemp < 0) std::cout << "nJtPtBinsTemp less than 0. please check input file. return 1" << std::endl;
  if(nJtAbsEtaBinsTemp < 0) std::cout << "nJtAbsEtaBinsTemp less than 0. please check input file. return 1" << std::endl;
  if(nIDTemp < 0) std::cout << "nIDTemp less than 0. please check input file. return 1" << std::endl;
  
  if(nCentBins < 0 || nJtPtBinsTemp < 0 || nJtAbsEtaBinsTemp < 0 || nIDTemp < 0){
    responseFile_p->Close();
    delete responseFile_p;
    return 1;
  }

  const Int_t nJtPtBins = nJtPtBinsTemp;
  Double_t jtPtBins[nJtPtBins+1];
  for(Int_t jI = 0; jI < nJtPtBins+1; ++jI){
    jtPtBins[jI] = jtPtBinsTemp.at(jI);
  }

  const Int_t nJtAbsEtaBins = nJtAbsEtaBinsTemp;
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

  const Int_t nID = nIDTemp;
  std::string idStr[nID];
  for(Int_t i = 0; i < nID; ++i){
    idStr[i] = idStrTemp.at(i);
  }
  

  const Int_t nManip = 4;
  const std::string recoTruncStr[nManip] = {"", "", "_RecoTrunc", "_RecoTrunc"};
  Bool_t renormX[nManip] = {true, false, true, false};

  const Int_t nJets = jetDirList.size(); 

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  checkMakeDir("pdfDir");

  for(Int_t jI = 0; jI < nJets; ++jI){
    std::string dirName = jetDirList.at(jI);
    dirName = dirName.substr(0, dirName.find("/"));

    for(Int_t cI = 0; cI < nCentBins; ++cI){
      std::string centStr = "PP";
      if(!isPP) centStr = "Cent" + std::to_string(centBinsLow.at(cI)) + "to" + std::to_string(centBinsHi.at(cI));

      for(Int_t iI = 0; iI < nID; ++iI){	

	for(Int_t modI = 0; modI < nResponseMod; ++modI){
	  std::string resStr = "ResponseMod" + prettyString(responseMod[modI], 2, true);

	  for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
	    const std::string jtAbsEtaStr = "AbsEta" + prettyString(jtAbsEtaBinsLow[aI], 1, true) + "to" + prettyString(jtAbsEtaBinsHi[aI], 1, true);	
	  
	    for(Int_t mI = 0; mI < nManip; ++mI){
	      TH2D* response_h = (TH2D*)responseFile_p->Get((dirName + "/response_" + dirName + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_" + jtAbsEtaStr + recoTruncStr[mI] + "_h").c_str());
	      
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
	      else{
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
	      
	      std::string renormStr = "renormX";
	      if(!renormX[mI]) renormStr = "renormY";
	      canv_p->SaveAs(("pdfDir/response_" + dirName + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_" + jtAbsEtaStr + recoTruncStr[mI] + "_" + renormStr + "_" + dateStr  + ".pdf").c_str());
	      delete canv_p;
	    }
	    
	    Int_t nPadX = -1;
	    Int_t nPadY = -1;
	    getPadsXY(nJtPtBins, &nPadX, &nPadY);
	    
	    TCanvas* recoJtPtPerGenPtBin_p = new TCanvas("recoJtPtPerGenPtBin_p", "", 450*nPadX, 450*nPadY);
	    recoJtPtPerGenPtBin_p->SetTopMargin(0.01);
	    recoJtPtPerGenPtBin_p->SetBottomMargin(0.01);
	    recoJtPtPerGenPtBin_p->SetLeftMargin(0.01);
	    recoJtPtPerGenPtBin_p->SetRightMargin(0.01);
	    
	    recoJtPtPerGenPtBin_p->Divide(nPadX, nPadY);
	    
	    for(Int_t jI = 0; jI < nJtPtBins; ++jI){
	      const std::string jtPtStr = "Pt" + prettyString(jtPtBins[jI], 1, true) + "to" + prettyString(jtPtBins[jI+1], 1, true);
	      
	      recoJtPtPerGenPtBin_p->cd();
	      recoJtPtPerGenPtBin_p->cd(jI+1);
	      
	      TH1D* temp_p = (TH1D*)responseFile_p->Get((dirName + "/recoJtPtPerGenPtBin_" + dirName + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_Gen" + jtPtStr + "_" + jtAbsEtaStr + "_h").c_str());
	      
	      temp_p->DrawCopy("HIST E1 P");
	    }
	    
	    recoJtPtPerGenPtBin_p->SaveAs(("pdfDir/recoJtPtPerGenPtBin_" + dirName + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_" + jtAbsEtaStr + "_" + dateStr + ".pdf").c_str());
	    delete recoJtPtPerGenPtBin_p;
	    
	    TCanvas* recoJtPtPerGenPtBinWeighted_p = new TCanvas("recoJtPtPerGenPtBinWeighted_p", "", 450*nPadX, 450*nPadY);
	    recoJtPtPerGenPtBinWeighted_p->SetTopMargin(0.01);
	    recoJtPtPerGenPtBinWeighted_p->SetBottomMargin(0.01);
	    recoJtPtPerGenPtBinWeighted_p->SetLeftMargin(0.01);
	    recoJtPtPerGenPtBinWeighted_p->SetRightMargin(0.01);
	    
	    recoJtPtPerGenPtBinWeighted_p->Divide(nPadX, nPadY);
	    
	    for(Int_t jI = 0; jI < nJtPtBins; ++jI){
	      const std::string jtPtStr = "Pt" + prettyString(jtPtBins[jI], 1, true) + "to" + prettyString(jtPtBins[jI+1], 1, true);
	    
	      recoJtPtPerGenPtBinWeighted_p->cd();
	      recoJtPtPerGenPtBinWeighted_p->cd(jI+1);
	      
	      
	      TH1D* temp_p = (TH1D*)responseFile_p->Get((dirName + "/recoJtPtPerGenPtBin_" + dirName + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_Gen" + jtPtStr + "_" + jtAbsEtaStr + "_Weighted_h").c_str());
	      
	      temp_p->DrawCopy("HIST E1 P");
	    }
	    
	    recoJtPtPerGenPtBinWeighted_p->SaveAs(("pdfDir/recoJtPtPerGenPtBin_Weighted_" + dirName + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_" + jtAbsEtaStr + "_" + dateStr + ".pdf").c_str());

	    delete recoJtPtPerGenPtBinWeighted_p;
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
    std::cout << "Usage: ./bin/plotJetResponse.exe <inResponseName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += plotJetResponse(argv[1]);

  std::cout << "Job complete. Return " << retVal << "." << std::endl;
  return retVal;
}
