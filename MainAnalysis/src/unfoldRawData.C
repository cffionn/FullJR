//cpp dependencies
#include <iostream>
#include <string>
#include <fstream>
#include <vector>

//ROOT dependencies
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TDatime.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"

//RooUnfold dependencies
#include "src/RooUnfoldResponse.h"
#include "src/RooUnfoldBayes.h"

//Local dependencies
#include "MainAnalysis/include/cutPropagator.h"

//Non-local FullJR dependencies (Utility, etc.)
#include "Utility/include/checkMakeDir.h"
#include "Utility/include/getLinBins.h"
#include "Utility/include/histDefUtility.h"
#include "Utility/include/mntToXRootdFileString.h"
#include "Utility/include/plotUtilities.h"
#include "Utility/include/returnRootFileContentsList.h"
#include "Utility/include/vanGoghPalette.h"

int unfoldRawData(const std::string inDataFileName, const std::string inResponseName, const std::string selectJtAlgo = "")
{
  vanGoghPalette vg;

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  const std::string dateStr2 = std::to_string(date->GetYear()) + "." + std::to_string(date->GetMonth()) + "." + std::to_string(date->GetDay());
  delete date;

  const std::string preDirPerma= "/data/cmcginn/FullJR/";
  std::string preDirTemp = "";
  if(checkDir(preDirPerma)) preDirTemp = preDirPerma;
  const std::string preDir = preDirTemp;

  checkMakeDir(preDir + "pdfDir");
  checkMakeDir(preDir + "pdfDir/" + dateStr);

  TFile* responseFile_p = new TFile(inResponseName.c_str(), "READ");
  std::vector<std::string> responseJetDirList = returnRootFileContentsList(responseFile_p, "TDirectoryFile", "JetAnalyzer");
  std::cout << "Printing " << responseJetDirList.size() << " response jets..." << std::endl;
  for(unsigned int jI = 0; jI < responseJetDirList.size(); ++jI){
    std::cout << " " << jI << "/" << responseJetDirList.size() << ": " << responseJetDirList.at(jI) << std::endl;

    checkMakeDir(preDir + "pdfDir/" + dateStr + "/" + responseJetDirList.at(jI));
  }

  cutPropagator cutPropResponse;
  cutPropResponse.Clean();
  cutPropResponse.GetAllVarFromFile(responseFile_p);

  responseFile_p->Close();
  delete responseFile_p;

  TFile* dataFile_p = new TFile(inDataFileName.c_str(), "READ");
  std::vector<std::string> dataJetDirList = returnRootFileContentsList(dataFile_p, "TDirectoryFile", "JetAnalyzer");
  std::cout << "Printing " << dataJetDirList.size() << " data jets..." << std::endl;
  for(unsigned int jI = 0; jI < dataJetDirList.size(); ++jI){
    std::cout << " " << jI << "/" << dataJetDirList.size() << ": " << dataJetDirList.at(jI) << std::endl;
  }

  cutPropagator cutPropData;
  cutPropData.Clean();
  cutPropData.GetAllVarFromFile(dataFile_p);

  dataFile_p->Close();
  delete dataFile_p;

  std::cout << cutPropResponse.GetPthats().size() << std::endl;

  if(!cutPropResponse.CheckPropagatorsMatch(cutPropData, false, true)){
    std::cout << "unfoldRawData - Cuts listed in data file \'" << inDataFileName << "\' and response file \'" << inResponseName << "\' do not match. return 1" << std::endl;
    return 1;
  }

  Int_t valForForLoops = 100000000;
  if(doLocalDebug || doGlobalDebug){
    std::cout << "DOLOCALDEBUG or DOGLOBALDEBUG in unfoldRawData: Setting all histogram array sizes to 1 for faster processing" << std::endl;
    valForForLoops = 1;
  }

  const Int_t isDataPP = cutPropData.GetIsPP();
  const Int_t isResponsePP = cutPropResponse.GetIsPP();

  const Int_t nSyst = TMath::Min(valForForLoops, cutPropResponse.GetNSyst());
  std::vector<std::string> systStr = cutPropResponse.GetSystStr();
  
  const Int_t nResponseMod = TMath::Min(valForForLoops, cutPropResponse.GetNResponseMod());
  std::vector<double> responseMod = cutPropResponse.GetResponseMod();
  std::vector<double> jerVarData = cutPropResponse.GetJERVarData();

  const Int_t nCentBins = TMath::Min(valForForLoops, cutPropData.GetNCentBins());
  std::vector<Int_t> centBinsLow = cutPropData.GetCentBinsLow();
  std::vector<Int_t> centBinsHi = cutPropData.GetCentBinsHi();

  const Int_t nJtPtBins = cutPropData.GetNJtPtBins();
  std::vector<Double_t> jtPtBinsTemp = cutPropData.GetJtPtBins();

  const Int_t nJtAbsEtaBins = TMath::Min(valForForLoops, cutPropData.GetNJtAbsEtaBins());
  std::vector<Double_t> jtAbsEtaBinsLowTemp = cutPropData.GetJtAbsEtaBinsLow();
  std::vector<Double_t> jtAbsEtaBinsHiTemp = cutPropData.GetJtAbsEtaBinsHi();

  const Int_t nID = TMath::Min(valForForLoops, cutPropData.GetNID());
  std::vector<std::string> idStr = cutPropData.GetIdStr();
  std::vector<double> jtPfCHMFCutLow = cutPropData.GetJtPfCHMFCutLow();
  std::vector<double> jtPfCHMFCutHi = cutPropData.GetJtPfCHMFCutHi();
  std::vector<double> jtPfMUMFCutLow = cutPropData.GetJtPfMUMFCutLow();
  std::vector<double> jtPfMUMFCutHi = cutPropData.GetJtPfMUMFCutHi();

  Double_t jtPtBins[nJtPtBins+1];
  std::cout << "nJtPtBins: ";
  for(Int_t jI = 0; jI < nJtPtBins+1; ++jI){
    jtPtBins[jI] = jtPtBinsTemp.at(jI);
    std::cout << " " << jtPtBins[jI] << ",";
  }
  std::cout << std::endl;

  Double_t jtAbsEtaBinsLow[nJtAbsEtaBins];
  Double_t jtAbsEtaBinsHi[nJtAbsEtaBins];
  std::cout << "nJtAbsEtaBins: ";
  for(Int_t jI = 0; jI < nJtAbsEtaBins; ++jI){
    jtAbsEtaBinsLow[jI] = jtAbsEtaBinsLowTemp.at(jI);
    jtAbsEtaBinsHi[jI] = jtAbsEtaBinsHiTemp.at(jI);
    std::cout << " " << jtAbsEtaBinsLow[jI] << "-" << jtAbsEtaBinsHi[jI] << ",";
  }
  std::cout << std::endl;

  std::cout << "nCentBins: " << nCentBins << std::endl;
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    std::cout << " " << cI << "/" << nCentBins << ": " << centBinsLow.at(cI) << "-" << centBinsHi.at(cI) << std::endl;
  }

  std::cout << "Raw data from file: \'" << inDataFileName << "\'" << std::endl;
  std::cout << "Response data from file: \'" << inResponseName << "\'" << std::endl;

  //Reduce to match jet dirs
  unsigned int pos = 0;
  while(pos < dataJetDirList.size()){
    bool isFound = false;
    for(unsigned int i = 0; i < responseJetDirList.size(); ++i){
      if(dataJetDirList.at(pos).size() != responseJetDirList.at(i).size()) continue;
      if(dataJetDirList.at(pos).find(responseJetDirList.at(i)) == std::string::npos) continue;
      
      isFound = true;
      break;
    }

    if(isFound) ++pos;
    else dataJetDirList.erase(dataJetDirList.begin()+pos);
  }

  pos = 0;
  while(pos < responseJetDirList.size()){
    bool isFound = false;
    for(unsigned int i = 0; i < dataJetDirList.size(); ++i){
      if(responseJetDirList.at(pos).size() != dataJetDirList.at(i).size()) continue;
      if(responseJetDirList.at(pos).find(dataJetDirList.at(i)) == std::string::npos) continue;
      
      isFound = true;
      break;
    }

    if(isFound) ++pos;
    else responseJetDirList.erase(responseJetDirList.begin()+pos);
  }

  std::cout << "Shared jets to process: " << std::endl;
  for(unsigned int i = 0; i < responseJetDirList.size(); ++i){
    std::cout << " " << i << "/" << responseJetDirList.size() << ": " << responseJetDirList.at(i) << std::endl;
  }

  if(selectJtAlgo.size() != 0){
    std::cout << "Restricting to \'" << selectJtAlgo << "\'." << std::endl;
    unsigned int pos = 0;
    while(pos < responseJetDirList.size()){
      if(responseJetDirList.at(pos).find(selectJtAlgo) == std::string::npos) responseJetDirList.erase(responseJetDirList.begin()+pos);
      else ++pos;
    }

    if(responseJetDirList.size() == 0){
      std::cout << "No jet dirs after selection for " << selectJtAlgo << std::endl;
      return 1;
    }
  }


  const Int_t nDataJet = responseJetDirList.size();

  std::string outFileName = inDataFileName;
  while(outFileName.find("/") != std::string::npos){outFileName.replace(0, outFileName.find("/")+1, "");}
  if(outFileName.find(".txt") != std::string::npos) outFileName.replace(outFileName.find(".txt"), std::string(".txt").size(), "");
  else if(outFileName.find(".root") != std::string::npos) outFileName.replace(outFileName.find(".root"), std::string(".root").size(), "");

  std::string outFileName2 = inResponseName;
  while(outFileName2.find("/") != std::string::npos){outFileName2.replace(0, outFileName2.find("/")+1, "");}
  if(outFileName2.find(".txt") != std::string::npos) outFileName2.replace(outFileName2.find(".txt"), std::string(".txt").size(), "");
  else if(outFileName2.find(".root") != std::string::npos) outFileName2.replace(outFileName2.find(".root"), std::string(".root").size(), "");

  std::string debugStr = "";
  if(doLocalDebug || doGlobalDebug) debugStr = "DEBUG_";

  std::string selectJtAlgoStr = selectJtAlgo;
  if(selectJtAlgoStr.size() != 0) selectJtAlgoStr = selectJtAlgoStr + "_";

  const Int_t sizeToTruncName = 40;
  while(outFileName.size() > sizeToTruncName){outFileName = outFileName.substr(0,outFileName.size()-1);}
  while(outFileName2.size() > sizeToTruncName){outFileName2 = outFileName2.substr(0,outFileName2.size()-1);}

  const Int_t nSuperBayes = 1;
  outFileName = "output/" + outFileName + "_" + outFileName2 + "_UnfoldRawData_NSuperBayes" + std::to_string(nSuperBayes) + "_" + selectJtAlgoStr + debugStr + dateStr + ".root";

  while(outFileName.find("__") != std::string::npos){outFileName.replace(outFileName.find("__"), 2, "_");}

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  //Following two lines necessary else you get absolutely clobbered on deletion timing. See threads here:
  //https://root-forum.cern.ch/t/closing-root-files-is-dead-slow-is-this-a-normal-thing/5273/16
  //https://root-forum.cern.ch/t/tfile-speed/17549/25
  //Bizarre
  outFile_p->SetBit(TFile::kDevNull);
  TH1::AddDirectory(kFALSE);

  const Int_t nBayes = 37;
  const Int_t nBigBayesSymm = 3;
  const Int_t nBayesBig = nBayes - (2*nBigBayesSymm + 1);
  Int_t temp100Pos = -1;
  Int_t bayesVal[nBayes];
  for(Int_t bI = 0; bI < nBayes; ++bI){
    if(bI >= nBayesBig) bayesVal[bI] = 100 + 1 + nBigBayesSymm - (nBayes - bI);
    else bayesVal[bI] = bI+1;

    if(bayesVal[bI] == 100) temp100Pos = bI;
  }


  const Int_t bayes100Pos = temp100Pos;

  std::cout << "Bayes: " << std::endl;
  for(Int_t bI = 0; bI < nBayes; ++bI){
    std::cout << " " << bI << "/" << nBayes << ": " << bayesVal[bI] << std::endl;
  }

  std::cout << "N Super Bayes: " << nSuperBayes << std::endl;

  //  return 0;

  const Int_t nHistDim = nDataJet*nCentBins*nID*nResponseMod*nJtAbsEtaBins*nSyst;
  std::vector<std::string> histTag;
  std::vector<int> histBestBayes;
  
  cutPropData.SetNBayes(nBayes);
  cutPropData.SetNBigBayesSymm(nBigBayesSymm);
  cutPropData.SetBayesVal(nBayes, bayesVal);
  cutPropData.SetNSuperBayes(nSuperBayes);

  cutPropData.SetNResponseMod(nResponseMod);
  cutPropData.SetResponseMod(responseMod);
  cutPropData.SetJERVarData(jerVarData);

  cutPropData.SetNSyst(nSyst);
  cutPropData.SetSystStr(systStr);

  const Int_t nBayesDraw = TMath::Min(valForForLoops, 6);
  TDirectory* dir_p[nDataJet];
  TH1D* jtPtUnfolded_h[nDataJet][nCentBins][nID][nResponseMod][nJtAbsEtaBins][nSyst][nBayes];
  TH1D* jtPtUnfolded_RecoTrunc_h[nDataJet][nCentBins][nID][nResponseMod][nJtAbsEtaBins][nSyst][nBayes];

  for(Int_t jI = 0; jI < nDataJet; ++jI){
    std::string tempStr = responseJetDirList.at(jI);
    if(tempStr.find("/") != std::string::npos) tempStr.replace(tempStr.find("/"), tempStr.size() - tempStr.find("/"), "");
    dir_p[jI] = (TDirectory*)outFile_p->mkdir(tempStr.c_str());

    for(Int_t cI = 0; cI < nCentBins; ++cI){
      std::string centStr = "PP";
      if(!isDataPP) centStr = "Cent" + std::to_string(centBinsLow.at(cI)) + "to" + std::to_string(centBinsHi.at(cI));
      
      for(Int_t idI = 0; idI < nID; ++idI){

	for(Int_t mI = 0; mI < nResponseMod; ++mI){
          const std::string resStr = "ResponseMod" + prettyString(responseMod[mI], 2, true);

	  for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
	    const std::string jtAbsEtaStr = "AbsEta" + prettyString(jtAbsEtaBinsLow[aI], 1, true) + "to" + prettyString(jtAbsEtaBinsHi[aI], 1, true);

	    for(Int_t sI = 0; sI < nSyst; ++sI){
	      const std::string tempSystStr = systStr.at(sI) + "_";

	      for(Int_t bI = 0; bI < nBayes; ++bI){
		std::string bayesStr = "Bayes" + std::to_string(bayesVal[bI]);
		
		jtPtUnfolded_h[jI][cI][idI][mI][aI][sI][bI] = new TH1D(("jtPtUnfolded_" + tempStr + "_" + centStr + "_" + idStr.at(idI) + "_" + resStr + "_" + jtAbsEtaStr + "_" + tempSystStr + bayesStr + "_h").c_str(), ";Unfolded Jet p_{T};Counts", nJtPtBins, jtPtBins);
		jtPtUnfolded_RecoTrunc_h[jI][cI][idI][mI][aI][sI][bI] = new TH1D(("jtPtUnfolded_RecoTrunc_" + tempStr + "_" + centStr + "_" + idStr.at(idI) + "_" + resStr + "_" + jtAbsEtaStr + "_" + tempSystStr + bayesStr + "_h").c_str(), ";Unfolded Jet p_{T};Counts", nJtPtBins, jtPtBins);
		centerTitles({jtPtUnfolded_h[jI][cI][idI][mI][aI][sI][bI], jtPtUnfolded_RecoTrunc_h[jI][cI][idI][mI][aI][sI][bI]});
		setSumW2({jtPtUnfolded_h[jI][cI][idI][mI][aI][sI][bI], jtPtUnfolded_RecoTrunc_h[jI][cI][idI][mI][aI][sI][bI]});
	      }
	    }
	  }
	}
      }
    }
  }

  responseFile_p = new TFile(inResponseName.c_str(), "READ");
  RooUnfoldResponse* rooResponse_RecoTrunc_h[nDataJet][nCentBins][nID][nResponseMod][nJtAbsEtaBins][nSyst];

  for(Int_t jI = 0; jI < nDataJet; ++jI){
    std::string tempStr = responseJetDirList.at(jI);
    if(tempStr.find("/") != std::string::npos) tempStr.replace(tempStr.find("/"), tempStr.size() - tempStr.find("/"), "");

    for(Int_t cI = 0; cI < nCentBins; ++cI){
      std::string centStr = "PP";
      if(!isResponsePP) centStr = "Cent" + std::to_string(centBinsLow.at(cI)) + "to" + std::to_string(centBinsHi.at(cI));

      for(Int_t idI = 0; idI < nID; ++idI){
	for(Int_t mI = 0; mI < nResponseMod; ++mI){
	  const std::string resStr = "ResponseMod" + prettyString(responseMod[mI], 2, true);

	  for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
	    const std::string jtAbsEtaStr = "AbsEta" + prettyString(jtAbsEtaBinsLow[aI], 1, true) + "to" + prettyString(jtAbsEtaBinsHi[aI], 1, true);
	    
	    for(Int_t sI = 0; sI < nSyst; ++sI){
	      const std::string tempSystStr = systStr.at(sI) + "_";
	     
	      rooResponse_RecoTrunc_h[jI][cI][idI][mI][aI][sI] = (RooUnfoldResponse*)responseFile_p->Get((tempStr + "/rooResponse_" + tempStr + "_" + centStr + "_" + idStr.at(idI) + "_" + resStr + "_" + jtAbsEtaStr + "_" + tempSystStr + "RecoTrunc_h").c_str());
	    }
	  }
	}
      }
    }
  }


  dataFile_p = new TFile(inDataFileName.c_str(), "READ");
  TH1D* jtPtRaw_RecoTrunc_h[nDataJet][nCentBins][nID][nJtAbsEtaBins];

  for(Int_t jI = 0; jI < nDataJet; ++jI){
    std::string tempStr = responseJetDirList.at(jI);
    if(tempStr.find("/") != std::string::npos) tempStr.replace(tempStr.find("/"), tempStr.size() - tempStr.find("/"), "");

    for(Int_t cI = 0; cI < nCentBins; ++cI){
      std::string centStr = "PP";
      if(!isDataPP) centStr = "Cent" + std::to_string(centBinsLow.at(cI)) + "to" + std::to_string(centBinsHi.at(cI));

      for(Int_t idI = 0; idI < nID; ++idI){
        for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
          const std::string jtAbsEtaStr = "AbsEta" + prettyString(jtAbsEtaBinsLow[aI], 1, true) + "to" + prettyString(jtAbsEtaBinsHi[aI], 1, true);

	  jtPtRaw_RecoTrunc_h[jI][cI][idI][aI] = (TH1D*)dataFile_p->Get((tempStr + "/jtPtRaw_RecoTrunc_" + tempStr + "_" + centStr + "_" + idStr.at(idI) + "_" + jtAbsEtaStr + "_h").c_str());
        }
      }
    }
  }

  std::cout << "Start Unfolding..." << std::endl;
  for(Int_t jI = 0; jI < nDataJet; ++jI){
    std::string tempStr = responseJetDirList.at(jI);
    if(tempStr.find("/") != std::string::npos) tempStr.replace(tempStr.find("/"), tempStr.size() - tempStr.find("/"), "");

    std::cout << " Unfolding " << jI << "/" << nDataJet << ": " << tempStr << std::endl;

    for(Int_t cI = 0; cI < nCentBins; ++cI){
      if(!isDataPP) std::cout << "  " << centBinsLow.at(cI) << "-" << centBinsHi.at(cI) << "%..." << std::endl;

      for(Int_t idI = 0; idI < nID; ++idI){

	for(Int_t mI = 0; mI < nResponseMod; ++mI){	 
	  for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
	    for(Int_t sI = 0; sI < nSyst; ++sI){

	      std::string histName = jtPtUnfolded_RecoTrunc_h[jI][cI][idI][mI][aI][sI][0]->GetName();
	      bool highLight = true;
	      if(histName.find("ak3PFJetAnalyzer") == std::string::npos) highLight = false;
	      else if(histName.find("NoID") == std::string::npos) highLight = false;
	      else if(histName.find("ResponseMod0p00") == std::string::npos) highLight = false;
	      else if(histName.find("AbsEta0p5to1p0") == std::string::npos) highLight = false;
	      else if(histName.find("JECUpMC") == std::string::npos) highLight = false;

	      RooUnfoldResponse* rooResSuperClone_p = (RooUnfoldResponse*)rooResponse_RecoTrunc_h[jI][cI][idI][mI][aI][sI]->Clone("rooResSuperClone_p");
	      TH2D* initRes_p = (TH2D*)rooResSuperClone_p->Hresponse()->Clone("initRes_p");
	      TH1D* initMeas_p = (TH1D*)rooResSuperClone_p->Hmeasured()->Clone("initMeas_p");
	      TH1D* initTrue_p = (TH1D*)rooResSuperClone_p->Htruth()->Clone("initTrue_p");

	      if(highLight){
		std::cout << "DOING PRE-HIGHLIGHT" << std::endl;
		std::cout << " Print initTrue: " << std::endl;
		initTrue_p->Print("ALL");

		std::cout << " Print initMeas: " << std::endl;
		initMeas_p->Print("ALL");

		std::cout << " Print from rooUnfold: " << std::endl;
		rooResSuperClone_p->Print("ALL");

		std::cout << " Print from TH2" << std::endl;
		for(Int_t bIY = 0; bIY < initRes_p->GetYaxis()->GetNbins(); ++bIY){
		  std::string firstString = std::to_string(bIY+1);
		  if(firstString.size() == 0) firstString = " " + firstString;
		  std::cout << firstString;

		  for(Int_t bIX = 0; bIX < initRes_p->GetXaxis()->GetNbins(); ++bIX){
		    std::string value = prettyString(initRes_p->GetBinContent(bIX+1, bIY+1), 6, false);
		    
		    while(value.size() < 10){value = " " + value;}
		    while(value.size() > 10){value.replace(value.size()-2, 1, "");}
		    value = "'" + value + "'";
		    std::cout << value;
		  }
		  std::cout << std::endl;
		}
	      }
	      
	      for(Int_t bsI = 0; bsI < nSuperBayes; ++bsI){
		RooUnfoldBayes superBayes(rooResSuperClone_p, jtPtRaw_RecoTrunc_h[jI][cI][idI][aI], 3, false, "name");
		superBayes.SetVerbose(-1);
		TH1D* unfold_h = (TH1D*)superBayes.Hreco(RooUnfold::kCovToy);
		Double_t tot = unfold_h->Integral();
		unfold_h->Scale(1./tot);
		
		for(Int_t bIY = 0; bIY < initRes_p->GetNbinsY(); ++bIY){
		  Double_t sumOverX = 0.0;
		  for(Int_t bIX = 0; bIX < initRes_p->GetNbinsX(); ++bIX){
		    sumOverX += initRes_p->GetBinContent(bIX+1, bIY+1);
		  }

		  Double_t scaleFactor = 1;
		  if(sumOverX > 0) scaleFactor = unfold_h->GetBinContent(bIY+1)/sumOverX;

		  //		  std::cout << "Check bin width match: " << unfold_h->GetBinLowEdge(bIY+1) << "-" << unfold_h->GetBinLowEdge(bIY+2) << ", " << initRes_p->GetYaxis()->GetBinLowEdge(bIY+1) << "-" << initRes_p->GetYaxis()->GetBinLowEdge(bIY+2) << "." << std::endl;

		  for(Int_t bIX = 0; bIX < initRes_p->GetNbinsX(); ++bIX){
		    initRes_p->SetBinContent(bIX+1, bIY+1, initRes_p->GetBinContent(bIX+1, bIY+1)*scaleFactor);
		    initRes_p->SetBinError(bIX+1, bIY+1, initRes_p->GetBinError(bIX+1, bIY+1)*scaleFactor);
		  }
		}
		
		delete rooResSuperClone_p;
		rooResSuperClone_p = NULL;
		rooResSuperClone_p = new RooUnfoldResponse(initMeas_p, unfold_h, initRes_p, "rooResSuperClone_p");
	      }

	      if(highLight){
		std::cout << "DOING POST-HIGHLIGHT" << std::endl;
		std::cout << " Print from roo" << std::endl;
		rooResSuperClone_p->Print("ALL");

		std::cout << " Print from TH2" << std::endl;
		for(Int_t bIY = 0; bIY < initRes_p->GetYaxis()->GetNbins(); ++bIY){
		  std::string firstString = std::to_string(bIY+1);
		  if(firstString.size() == 0) firstString = " " + firstString;
		  std::cout << firstString;

		  for(Int_t bIX = 0; bIX < initRes_p->GetXaxis()->GetNbins(); ++bIX){
		    std::string value = prettyString(initRes_p->GetBinContent(bIX+1, bIY+1), 6, false);
		    
		    while(value.size() < 10){value = " " + value;}
		    while(value.size() > 10){value.replace(value.size()-2, 1, "");}
		    value = "'" + value + "'";
		    std::cout << value;
		  }
		  std::cout << std::endl;
		}
	      }
	    
	      delete initRes_p;
	      delete initMeas_p;
	      delete initTrue_p;

	      for(Int_t bI = 0; bI < nBayes; ++bI){	    		
		RooUnfoldResponse* rooResClone_p = (RooUnfoldResponse*)rooResSuperClone_p->Clone("rooResClone_p");

		RooUnfoldBayes bayes(rooResClone_p, jtPtRaw_RecoTrunc_h[jI][cI][idI][aI], bayesVal[bI], false, ("name_" + std::to_string(bI)).c_str());
		bayes.SetVerbose(-1);	    
		TH1D* unfold_h = (TH1D*)bayes.Hreco(RooUnfold::kCovToy);	  
		
		for(Int_t bIX = 0; bIX < unfold_h->GetNbinsX(); ++bIX){
		  jtPtUnfolded_RecoTrunc_h[jI][cI][idI][mI][aI][sI][bI]->SetBinContent(bIX+1, unfold_h->GetBinContent(bIX+1));
		  jtPtUnfolded_RecoTrunc_h[jI][cI][idI][mI][aI][sI][bI]->SetBinError(bIX+1, unfold_h->GetBinError(bIX+1));
		}
		delete rooResClone_p;
	      }
	      delete rooResSuperClone_p;
	    }
	  }
	}
      }
    }
  }

  std::cout << "Unfolding complete." << std::endl;

  dataFile_p->Close();
  delete dataFile_p;

  responseFile_p->Close();
  delete responseFile_p;

  outFile_p->cd();

  const Double_t yPadFrac = 0.45;
  const Double_t marg = 0.12;

  const Int_t nStyles = 6;
  const Int_t styles[nStyles] = {24, 25, 27, 28, 46, 44};

  const Int_t nColors = 5;
  const Int_t colors[nColors] = {1, vg.getColor(0), vg.getColor(1), vg.getColor(2), vg.getColor(4)};


  std::vector<std::vector<std::string> > pdfNames;

  for(Int_t jI = 0; jI < nDataJet; ++jI){
    outFile_p->cd();
    dir_p[jI]->cd();

    std::string tempStr = responseJetDirList.at(jI);
    if(tempStr.find("/") != std::string::npos) tempStr.replace(tempStr.find("/"), tempStr.size() - tempStr.find("/"), "");

    Double_t lowPtTruncVal = 200.;
    Bool_t isBigJet = tempStr.find("ak8") != std::string::npos || tempStr.find("ak10") != std::string::npos || tempStr.find("akCs8") != std::string::npos || tempStr.find("akCs10") != std::string::npos;
    if(isBigJet) lowPtTruncVal = 300.;
    
    for(Int_t cI = 0; cI < nCentBins; ++cI){
      std::string centStr = "PP";
      std::string centStr2 = "PP";
      if(!isDataPP){
	centStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);
	centStr2 = std::to_string(centBinsLow[cI]) + "-" + std::to_string(centBinsHi[cI]) + "%";
      }

      for(Int_t idI = 0; idI < nID; ++idI){
	for(Int_t mI = 0; mI < nResponseMod; ++mI){	
	  const std::string resStr = "ResponseMod" + prettyString(responseMod[mI], 2, true);

	  for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
	    const std::string jtAbsEtaStr = "AbsEta" + prettyString(jtAbsEtaBinsLow[aI], 1, true) + "to" + prettyString(jtAbsEtaBinsHi[aI], 1, true);
	    //	    if(jtAbsEtaStr.find("AbsEta0p0to2p0") == std::string::npos) continue;
	    const std::string jtAbsEtaStr2 = prettyString(jtAbsEtaBinsLow[aI], 1, false) + "<|#eta|<" +  prettyString(jtAbsEtaBinsHi[aI], 1, false);

	    
	    pdfNames.push_back({});

	    for(Int_t sI = 0; sI < nSyst; ++sI){	    	   
              const std::string tempSystStr = systStr.at(sI) + "_";

	      TLegend* leg_p = new TLegend(0.15, 0.05, 0.5, 0.6);
	      leg_p->SetBorderSize(0.0);
	      leg_p->SetFillStyle(0);
	      leg_p->SetFillColor(0);
	      leg_p->SetTextFont(43);
	      leg_p->SetTextSize(14);

	      TLatex* label_p = new TLatex();
	      label_p->SetTextFont(43);
	      label_p->SetTextSize(14);
	      label_p->SetNDC();                  

	      TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);
	      TPad* pads[3];
	      canv_p->cd();
	      pads[0] = new TPad("pad0", "", 0.0, yPadFrac, 1.0, 1.0);
	      pads[0]->SetLeftMargin(marg);
	      pads[0]->SetTopMargin(0.01);
	      pads[0]->SetBottomMargin(0.001);
	      pads[0]->SetRightMargin(0.002);
	      pads[0]->Draw();

	      canv_p->cd();
	      pads[1] = new TPad("pad1", "", 0.0, yPadFrac - (yPadFrac - marg)/2., 1.0, yPadFrac);
	      pads[1]->SetLeftMargin(marg);
	      pads[1]->SetTopMargin(0.001);
	      pads[1]->SetBottomMargin(0.001);
	      pads[1]->SetRightMargin(0.002);
	      pads[1]->Draw();

	      canv_p->cd();
	      pads[2] = new TPad("pad2", "", 0.0, 0.0, 1.0, yPadFrac - (yPadFrac - marg)/2.);
	      pads[2]->SetLeftMargin(marg);
	      pads[2]->SetTopMargin(0.001);
	      pads[2]->SetBottomMargin(marg/(yPadFrac - (yPadFrac - marg)/2.));
	      pads[2]->SetRightMargin(0.002);
	      pads[2]->Draw();

	      canv_p->cd();	       
	      pads[0]->cd();

	      Double_t min = 1000000000;
	      Double_t max = -1;

	      for(Int_t bI = 0; bI < nBayesDraw; ++bI){
		for(Int_t bIX = 0; bIX < jtPtUnfolded_RecoTrunc_h[jI][cI][idI][mI][aI][sI][bI]->GetNbinsX(); ++bIX){
		  if(max < jtPtUnfolded_RecoTrunc_h[jI][cI][idI][mI][aI][sI][bI]->GetBinContent(bIX+1)) max = jtPtUnfolded_RecoTrunc_h[jI][cI][idI][mI][aI][sI][bI]->GetBinContent(bIX+1);
		  if(min > jtPtUnfolded_RecoTrunc_h[jI][cI][idI][mI][aI][sI][bI]->GetBinContent(bIX+1) && jtPtUnfolded_RecoTrunc_h[jI][cI][idI][mI][aI][sI][bI]->GetBinContent(bIX+1) > 0) min = jtPtUnfolded_RecoTrunc_h[jI][cI][idI][mI][aI][sI][bI]->GetBinContent(bIX+1);
		}
	      }

	      

	      for(Int_t bI = 0; bI < nBayes; ++bI){
		jtPtUnfolded_RecoTrunc_h[jI][cI][idI][mI][aI][sI][bI]->SetMaximum(20.*max);
		jtPtUnfolded_RecoTrunc_h[jI][cI][idI][mI][aI][sI][bI]->SetMinimum(min/20.);

		jtPtUnfolded_RecoTrunc_h[jI][cI][idI][mI][aI][sI][bI]->SetMarkerStyle(styles[bI%nStyles]);
		jtPtUnfolded_RecoTrunc_h[jI][cI][idI][mI][aI][sI][bI]->SetMarkerColor(colors[bI%nColors]);
		jtPtUnfolded_RecoTrunc_h[jI][cI][idI][mI][aI][sI][bI]->SetLineColor(colors[bI%nColors]);
		jtPtUnfolded_RecoTrunc_h[jI][cI][idI][mI][aI][sI][bI]->SetMarkerSize(0.8);

		jtPtUnfolded_RecoTrunc_h[jI][cI][idI][mI][aI][sI][bI]->GetXaxis()->SetTitleFont(43);
		jtPtUnfolded_RecoTrunc_h[jI][cI][idI][mI][aI][sI][bI]->GetYaxis()->SetTitleFont(43);
		jtPtUnfolded_RecoTrunc_h[jI][cI][idI][mI][aI][sI][bI]->GetXaxis()->SetLabelFont(43);
		jtPtUnfolded_RecoTrunc_h[jI][cI][idI][mI][aI][sI][bI]->GetYaxis()->SetLabelFont(43);

		jtPtUnfolded_RecoTrunc_h[jI][cI][idI][mI][aI][sI][bI]->GetXaxis()->SetTitleSize(14);
		jtPtUnfolded_RecoTrunc_h[jI][cI][idI][mI][aI][sI][bI]->GetYaxis()->SetTitleSize(14);
		jtPtUnfolded_RecoTrunc_h[jI][cI][idI][mI][aI][sI][bI]->GetXaxis()->SetLabelSize(14);
		jtPtUnfolded_RecoTrunc_h[jI][cI][idI][mI][aI][sI][bI]->GetYaxis()->SetLabelSize(14);

		jtPtUnfolded_RecoTrunc_h[jI][cI][idI][mI][aI][sI][bI]->GetXaxis()->SetTitleOffset(5.);
		jtPtUnfolded_RecoTrunc_h[jI][cI][idI][mI][aI][sI][bI]->GetYaxis()->SetTitleOffset(1.5);

		if(bI < nBayesDraw){
		  if(bI == 0){
		    jtPtUnfolded_RecoTrunc_h[jI][cI][idI][mI][aI][sI][bI]->DrawCopy("HIST E1 P");

		    jtPtRaw_RecoTrunc_h[jI][cI][idI][aI]->SetMarkerStyle(1);
		    jtPtRaw_RecoTrunc_h[jI][cI][idI][aI]->SetMarkerSize(0.001);
		    jtPtRaw_RecoTrunc_h[jI][cI][idI][aI]->SetMarkerColor(0);
		    jtPtRaw_RecoTrunc_h[jI][cI][idI][aI]->SetLineColor(1);
		    jtPtRaw_RecoTrunc_h[jI][cI][idI][aI]->SetLineWidth(2);
		    jtPtRaw_RecoTrunc_h[jI][cI][idI][aI]->DrawCopy("HIST E1 SAME");
		    
		    leg_p->AddEntry(jtPtRaw_RecoTrunc_h[jI][cI][idI][aI], "Folded", "L");
		  }
		  else jtPtUnfolded_RecoTrunc_h[jI][cI][idI][mI][aI][sI][bI]->DrawCopy("HIST E1 P SAME"); 		
		  leg_p->AddEntry(jtPtUnfolded_RecoTrunc_h[jI][cI][idI][mI][aI][sI][bI], ("Bayes=" + std::to_string(bayesVal[bI])).c_str(), "P L");

		}
	      }

	      canv_p->cd();
	      pads[0]->cd();
	      gStyle->SetOptStat(0);
	      gPad->SetLogy();
	      gPad->SetTicks(1,2);
	      bool doLogX = false;
	      if(jtPtUnfolded_RecoTrunc_h[jI][cI][idI][mI][aI][sI][0]->GetBinWidth(1)*3 < jtPtUnfolded_RecoTrunc_h[jI][cI][idI][mI][aI][sI][0]->GetBinWidth(jtPtUnfolded_RecoTrunc_h[jI][cI][idI][mI][aI][sI][0]->GetNbinsX()-1)) doLogX = true;

	      if(doLogX) gPad->SetLogx();

	      leg_p->Draw("SAME");

	      gPad->RedrawAxis();

	      
	      label_p->DrawLatex(0.55, 0.94, tempStr.c_str());
	      label_p->DrawLatex(0.55, 0.88, centStr2.c_str());
	      label_p->DrawLatex(0.55, 0.82, jtAbsEtaStr2.c_str());
	      label_p->DrawLatex(0.55, 0.76, resStr.c_str());
	      label_p->DrawLatex(0.55, 0.70, idStr.at(idI).c_str());
	      label_p->DrawLatex(0.55, 0.64, systStr.at(sI).c_str());

	      canv_p->cd();
              pads[1]->cd();

	      TH1D* clones_p[nBayes];
	      
	      max = -1;
	      min = 100;

	      for(Int_t bI = 0; bI < nBayes; ++bI){
		clones_p[bI] = (TH1D*)jtPtUnfolded_RecoTrunc_h[jI][cI][idI][mI][aI][sI][bI]->Clone(("clone_" + std::to_string(bI)).c_str());
		clones_p[bI]->Divide(jtPtUnfolded_RecoTrunc_h[jI][cI][idI][mI][aI][sI][0]);
		for(Int_t bIX = 0; bIX < clones_p[bI]->GetNbinsX(); ++bIX){
		  double binContent = clones_p[bI]->GetBinContent(bIX+1);
		  if(binContent > max) max = binContent;
		  if(binContent < min && binContent > 0) min = binContent;
		}
	      }

	      double interval = (max - min)/10.;
	      max += interval;
	      min -= interval;

	      for(Int_t bI = 0; bI < nBayes; ++bI){
		clones_p[bI]->SetMaximum(max);
		clones_p[bI]->SetMinimum(min);

		clones_p[bI]->GetYaxis()->SetTitle("Ratio w/ Bayes0");
		clones_p[bI]->GetYaxis()->SetTitleSize(9);
		clones_p[bI]->GetYaxis()->SetLabelSize(9);
		clones_p[bI]->GetYaxis()->SetTitleOffset(clones_p[bI]->GetYaxis()->GetTitleOffset()*1.5);
		clones_p[bI]->GetYaxis()->SetNdivisions(505);	       		

		if(bI < nBayesDraw){
		  if(bI == 0) clones_p[bI]->DrawCopy("HIST E1 P");
		  else clones_p[bI]->DrawCopy("HIST E1 P SAME");
		}
	      }

	      gStyle->SetOptStat(0);
	      if(doLogX) gPad->SetLogx();

	      for(Int_t bI = 0; bI < nBayes; ++bI){
		delete clones_p[bI];
		clones_p[bI] = NULL;
	      }
	      gPad->SetTicks(1,2);
	      gPad->RedrawAxis();

	      canv_p->cd();
              pads[2]->cd();
	      
	      int terminalPos = -1;

	      if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
	    
	      TH1D* bandValLow_p = (TH1D*)jtPtUnfolded_RecoTrunc_h[jI][cI][idI][mI][aI][sI][bayes100Pos]->Clone("bandValLow");
	      TH1D* bandValHi_p = (TH1D*)jtPtUnfolded_RecoTrunc_h[jI][cI][idI][mI][aI][sI][bayes100Pos]->Clone("bandValHi");

	      for(Int_t bI = bayes100Pos-nBigBayesSymm; bI <= bayes100Pos+nBigBayesSymm; ++bI){
		for(Int_t bIX = 0; bIX < bandValLow_p->GetNbinsX(); ++bIX){
		  if(jtPtUnfolded_RecoTrunc_h[jI][cI][idI][mI][aI][sI][bI]->GetBinContent(bIX+1) < bandValLow_p->GetBinContent(bIX+1)){
		    bandValLow_p->SetBinContent(bIX+1, jtPtUnfolded_RecoTrunc_h[jI][cI][idI][mI][aI][sI][bI]->GetBinContent(bIX+1));
		  }
		  if(jtPtUnfolded_RecoTrunc_h[jI][cI][idI][mI][aI][sI][bI]->GetBinContent(bIX+1) > bandValHi_p->GetBinContent(bIX+1)){
		    bandValHi_p->SetBinContent(bIX+1, jtPtUnfolded_RecoTrunc_h[jI][cI][idI][mI][aI][sI][bI]->GetBinContent(bIX+1));
		  }
		}
	      }

	      if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

	      max = -1;
	      min = 100;
	      for(Int_t bI = 0; bI < nBayes; ++bI){
		clones_p[bI] = (TH1D*)jtPtUnfolded_RecoTrunc_h[jI][cI][idI][mI][aI][sI][bI]->Clone(("clone_" + std::to_string(bI)).c_str());
		
		bool sub1Perc = true;
		for(Int_t bIX = 0; bIX < clones_p[bI]->GetNbinsX(); ++bIX){
		  double binCenter = (clones_p[bI]->GetBinLowEdge(bIX+1) + clones_p[bI]->GetBinLowEdge(bIX+2))/2.;
		  double binContent = clones_p[bI]->GetBinContent(bIX+1);
		  double bandLowContent = bandValLow_p->GetBinContent(bIX+1);
		  double bandHiContent = bandValLow_p->GetBinContent(bIX+1);

		  if(binCenter > lowPtTruncVal && binCenter < 1000. && sub1Perc){
		    if(binContent/bandHiContent > 1.01) sub1Perc = false;
		    else if(binContent/bandLowContent < .99) sub1Perc = false;

		    //From when dividing just by 100 pos rather than considering covergence band
		    //		    if(binContent > 1.01 || binContent < .99) sub1Perc = false;
		  }

		  if(binContent > max) max = binContent;
		  if(binContent < min && binContent > 0) min = binContent;
		}

		if(sub1Perc && terminalPos < 0 && bI < nBayesBig) terminalPos = bI;
	      }

	      if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

	      delete bandValLow_p;
	      delete bandValHi_p;

	      if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

	      interval = (max - min)/10.;
	      max += interval;
	      min -= interval;

	      const Int_t nTempBins = 10;
	      Double_t tempBins[nTempBins+1];
	      getLinBins(min, max, nTempBins, tempBins);

	      if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

	      for(Int_t bI = 1; bI < nBayes; ++bI){
 		if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << ", " << bI << std::endl;

		clones_p[bI]->SetMaximum(max);
		clones_p[bI]->SetMinimum(min);

 		if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << ", " << bI << std::endl;
		
		clones_p[bI]->GetYaxis()->SetTitle("Ratio w/ Bayes 100");
		clones_p[bI]->GetYaxis()->SetTitleSize(9);
		clones_p[bI]->GetYaxis()->SetLabelSize(9);
		clones_p[bI]->GetYaxis()->SetTitleOffset(clones_p[bI]->GetYaxis()->GetTitleOffset()*1.5);
		clones_p[bI]->GetYaxis()->SetNdivisions(505);

		if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << ", " << bI << std::endl;

		if(bI < nBayesDraw){
		  if(bI == 1) clones_p[bI]->DrawCopy("HIST E1 P");
		  else clones_p[bI]->DrawCopy("HIST E1 P SAME");
		}

		if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << ", " << bI << std::endl;
	      }

	      if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

	      gStyle->SetOptStat(0);
	      if(doLogX) gPad->SetLogx();

	      for(Int_t bI = 0; bI < nBayes; ++bI){
		delete clones_p[bI];
		clones_p[bI] = NULL;
	      }

	      drawWhiteBox(900, 1100, .00, min-0.0001);

	      gPad->SetTicks(1,2);
	      gPad->RedrawAxis();

	      //	      label_p->SetNDC(0);
	      canv_p->cd();
	      pads[0]->cd();

	      if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
	  
	      const std::string tempHistTag = tempStr + "_" + centStr + "_" + idStr.at(idI) + "_ResponseMod" + prettyString(responseMod[mI], 2, true) + "_" + jtAbsEtaStr + "_" + systStr.at(sI);
	      histTag.push_back(tempHistTag);

	      if(terminalPos >= 0) histBestBayes.push_back(bayesVal[terminalPos]);
	      else histBestBayes.push_back(-1);
	      
	      if(terminalPos > 0){
		//		label_p->DrawLatex(110, tempBins[8], ("Terminate at Bayes=" + std::to_string(terminalPos+1)).c_str());
		label_p->DrawLatex(0.3, 0.4, ("Terminate at Bayes=" + std::to_string(bayesVal[terminalPos]+1)).c_str());
		
		Double_t maxDeltaCenter = -1;
		Double_t maxDelta = 0;
		Double_t lastDeltaCenter = -1;
		Double_t lastDelta = 0;
		for(Int_t bIX = 0; bIX < jtPtUnfolded_RecoTrunc_h[jI][cI][idI][mI][aI][sI][terminalPos]->GetNbinsX(); ++bIX){
		  double center = (jtPtUnfolded_RecoTrunc_h[jI][cI][idI][mI][aI][sI][terminalPos]->GetBinLowEdge(bIX+1) + jtPtUnfolded_RecoTrunc_h[jI][cI][idI][mI][aI][sI][terminalPos]->GetBinLowEdge(bIX+2))/2.;
		  
		  if(center < lowPtTruncVal) continue;
		  if(center > 1000.) continue;
		  
		  double content1 = jtPtUnfolded_RecoTrunc_h[jI][cI][idI][mI][aI][sI][terminalPos]->GetBinContent(bIX+1);
		  double content1Max = jtPtUnfolded_RecoTrunc_h[jI][cI][idI][mI][aI][sI][nBayes-1]->GetBinContent(bIX+1);
		  double content1MaxMin1 = jtPtUnfolded_RecoTrunc_h[jI][cI][idI][mI][aI][sI][nBayes-2]->GetBinContent(bIX+1);

		  if(content1 == 0 && content1Max != 0){
		    maxDelta = 100;
		    maxDeltaCenter = center;
		  }
		  else if(TMath::Abs(content1 - content1Max)/content1 > maxDelta){
		    maxDelta = TMath::Abs(content1 - content1Max)/content1;
		    maxDeltaCenter = center;
		  }

		  if(content1MaxMin1 == 0 && content1Max != 0){
		    lastDelta = 100;
		    lastDeltaCenter = center;
		  }
		  else if(TMath::Abs(content1MaxMin1 - content1Max)/content1MaxMin1 > lastDelta){
		    lastDelta = TMath::Abs(content1MaxMin1 - content1Max)/content1MaxMin1;
		    maxDeltaCenter = center;
		  }
		}

	      if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
	      
		std::cout << "MaxDelta: " << maxDelta << ", " << maxDeltaCenter << std::endl;
		std::cout << "LastDelta: " << lastDelta << ", " << lastDeltaCenter << std::endl;

		//		label_p->DrawLatex(110, tempBins[6], ("MaxDelta: " + prettyString(maxDelta*100, 2, false) + "%").c_str());
		//		label_p->DrawLatex(110, tempBins[4], ("LastDelta: " + prettyString(lastDelta*100, 2, false) + "%").c_str());
		label_p->DrawLatex(0.3, 0.3, ("MaxDelta: " + prettyString(maxDelta*100, 2, false) + "%").c_str());
		label_p->DrawLatex(0.3, 0.2, ("LastDelta: " + prettyString(lastDelta*100, 2, false) + "%").c_str());
	      }
	      else label_p->DrawLatex(0.3, 0.4, ("Doesn't terminate for nBayes=" + std::to_string(nBayes+1)).c_str());
	      //	      else label_p->DrawLatex(110, (max + min)/2., ("Doesn't terminate for nBayes=" + std::to_string(nBayes+1)).c_str());

	      canv_p->cd();
	      pads[2]->cd();
	      
	      label_p->SetNDC(0);

	      label_p->DrawLatex(100, min - interval*2, "100");
	      label_p->DrawLatex(200, min - interval*2, "200");
	      label_p->DrawLatex(400, min - interval*2, "400");
	      label_p->DrawLatex(600, min - interval*2, "600");
	      label_p->DrawLatex(1000, min - interval*2, "1000");

	      const std::string saveName = "jtPtUnfolded_" + tempStr + "_" + centStr + "_" + idStr.at(idI) + "_" + resStr + "_" + jtAbsEtaStr + "_" + tempSystStr +  "AllBayes_RecoTrunc_" + debugStr + dateStr + ".pdf";
	      pdfNames.at(pdfNames.size()-1).push_back(saveName);
	      const std::string finalSaveName = preDir + "pdfDir/" + dateStr + "/" + responseJetDirList.at(jI) + "/" + saveName;
	      quietSaveAs(canv_p, finalSaveName);
	    
	      delete pads[0];
	      delete pads[1];
	      delete pads[2];
 	      delete canv_p;
	      delete leg_p;
	      delete label_p;
	    }
	  }

	  
	  for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
	    for(Int_t sI = 0; sI < nSyst; ++sI){	    	   
	      for(Int_t bI = 0; bI < nBayes; ++bI){
		jtPtUnfolded_h[jI][cI][idI][mI][aI][sI][bI]->Write("", TObject::kOverwrite);
		delete jtPtUnfolded_h[jI][cI][idI][mI][aI][sI][bI];
		
		jtPtUnfolded_RecoTrunc_h[jI][cI][idI][mI][aI][sI][bI]->Write("", TObject::kOverwrite);
		delete jtPtUnfolded_RecoTrunc_h[jI][cI][idI][mI][aI][sI][bI];		
	      }
	    }
	  }
	}
      }
    }
  }

	      if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  //TEX FILE
  const std::string textWidth = "0.24";
  std::string texFileName = outFileName;
  texFileName.replace(texFileName.find(".root"), std::string(".root").size(), ".tex");
  texFileName.replace(0, std::string("output").size(), preDir + "pdfDir/" + dateStr);

  std::ofstream texFile(texFileName.c_str());

  texFile << "\\RequirePackage{xspace}" << std::endl;
  texFile << "\\RequirePackage{amsmath}" << std::endl;
  texFile << std::endl;

  texFile << "\\documentclass[xcolor=dvipsnames]{beamer}" << std::endl;
  texFile << "\\usetheme{Warsaw}" << std::endl;
  texFile << "\\setbeamercolor{structure}{fg=NavyBlue!90!NavyBlue}" << std::endl;
  texFile << "\\setbeamercolor{footlinecolor}{fg=white,bg=lightgray}" << std::endl;
  texFile << std::endl;

  texFile << "\\newcommand{\\pt}{\\ensuremath{p_{\\mathrm{T}}}\\xspace}" << std::endl;
  texFile << std::endl;

  texFile << "\\setbeamersize{text margin left=3pt,text margin right=3pt}" << std::endl;
  texFile << std::endl;

  texFile << "\\setbeamerfont{frametitle}{size=\\tiny}" << std::endl;
  texFile << "\\setbeamertemplate{frametitle}" << std::endl;
  texFile << "{" << std::endl;
  texFile << "    \\nointerlineskip" << std::endl;
  texFile << "    \\begin{beamercolorbox}[sep=0.05cm, ht=1.0em, wd=\\paperwidth]{frametitle}" << std::endl;
  texFile << "        \\vbox{}\\vskip-2ex%" << std::endl;
  texFile << "        \\strut\\insertframetitle\\strut" << std::endl;
  texFile << "        \\vskip-0.8ex%" << std::endl;
  texFile << "    \\end{beamercolorbox}" << std::endl;
  texFile << "}" << std::endl;
  texFile << std::endl;

  texFile << "\\setbeamertemplate{footline}{%" << std::endl;
  texFile << "  \\begin{beamercolorbox}[sep=.4em,wd=\\paperwidth,leftskip=0.5cm,rightskip=0.5cm]{footlinecolor}" << std::endl;
  texFile << "    \\hspace{0.075cm}%" << std::endl;
  texFile << "    \\hfill\\insertauthor \\hfill\\insertpagenumber" << std::endl;
  texFile << "  \\end{beamercolorbox}%" << std::endl;
  texFile << "}" << std::endl;
  texFile << "\\setbeamertemplate{navigation symbols}{}" << std::endl;
  texFile << std::endl;
  
  texFile << "\\setbeamertemplate{itemize item}[circle]" << std::endl;
  texFile << "\\setbeamertemplate{itemize subitem}[circle]" << std::endl;
  texFile << "\\setbeamertemplate{itemize subsubitem}[circle]" << std::endl;
  texFile << "\\setbeamercolor{itemize item}{fg=black}" << std::endl;
  texFile << "\\setbeamercolor{itemize subitem}{fg=black}" << std::endl;
  texFile << "\\setbeamercolor{itemize subsubitem}{fg=black}" << std::endl;
  texFile << std::endl;

  texFile << "\\definecolor{links}{HTML}{00BFFF}" << std::endl;
  texFile << "\\hypersetup{colorlinks,linkcolor=,urlcolor=links}" << std::endl;
  texFile << std::endl;

  texFile << "\\author[CM]{Chris McGinn}" << std::endl;
  texFile << std::endl;

  texFile << "\\begin{document}" << std::endl;
  texFile << "\\begin{frame}" << std::endl;
  texFile << "\\frametitle{\\centerline{Unfolding termination (" << dateStr2 << ")}}" << std::endl;
  texFile << " \\begin{itemize}" << std::endl;
  texFile << "  \\fontsize{10}{10}\\selectfont" << std::endl;
  texFile << "  \\item{Placeholder}" << std::endl;
  texFile << "  \\begin{itemize}" << std::endl;
  texFile << "   \\fontsize{10}{10}\\selectfont" << std::endl;
  texFile << "   \\item{Placeholder}" << std::endl;
  texFile << "  \\end{itemize}" << std::endl;
  texFile << " \\end{itemize}" << std::endl;
  texFile << "\\end{frame}" << std::endl;


  if(doLocalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  for(unsigned int spI = 0; spI < pdfNames.size(); ++spI){
    if(pdfNames.at(spI).at(0).find("AbsEta0p0to2p0") == std::string::npos) continue;

    std::string centStr = "PP";

    if(!isDataPP){
      centStr = pdfNames.at(spI).at(0).substr(pdfNames.at(spI).at(0).find("Cent"), pdfNames.at(spI).at(0).size());
      centStr.replace(centStr.find("_"), centStr.size(), "");
    }

    std::string jtStr = pdfNames.at(spI).at(0).substr(pdfNames.at(spI).at(0).find("_ak")+1, pdfNames.at(spI).at(0).size());
    jtStr.replace(jtStr.find("_"), jtStr.size(), "");

    std::string resStr = pdfNames.at(spI).at(0).substr(pdfNames.at(spI).at(0).find("_Response")+1, pdfNames.at(spI).at(0).size());
    resStr.replace(resStr.find("_"), resStr.size(), "");

    std::string idStr = "";
    if(pdfNames.at(spI).at(0).find("_NoID") != std::string::npos) idStr = "NoID";
    else if(pdfNames.at(spI).at(0).find("_LightMUID") != std::string::npos) idStr = "LightMUID";
    else if(pdfNames.at(spI).at(0).find("_LightMUAndCHID") != std::string::npos) idStr = "LightMUAndCHID";
    else if(pdfNames.at(spI).at(0).find("_FullLight") != std::string::npos) idStr = "FullLight";
    else if(pdfNames.at(spI).at(0).find("_FullTight") != std::string::npos) idStr = "FullTight";

    if(idStr.find("LightMUAndCHID") == std::string::npos) continue;
    if(resStr.find("0p10") == std::string::npos) continue;
    
    texFile << "\\begin{frame}" << std::endl;
    texFile << "\\frametitle{\\centerline{" << jtStr << ", " << centStr << ", " << resStr << ", " << idStr << "}}" << std::endl;

    for(unsigned int mI = 0; mI < pdfNames.at(spI).size(); ++mI){
      texFile << "\\includegraphics[width=" << textWidth << "\\textwidth]{" << pdfNames.at(spI).at(mI) << "}";
      if(mI == 3 || mI == 7) texFile << "\\\\";
      texFile << std::endl;
    }

    texFile << "\\begin{itemize}" << std::endl;
    texFile << "\\fontsize{8}{8}\\selectfont" << std::endl;
    texFile << "\\item{test}" << std::endl;
    texFile << "\\end{itemize}" << std::endl;
    texFile << "\\end{frame}" << std::endl;
  }


  texFile << "\\end{document}" << std::endl;
  texFile << std::endl;

  texFile.close();


  outFile_p->cd();
  TDirectory* cutDir_p = (TDirectory*)outFile_p->mkdir("cutDir");
  TDirectory* subDir_p = (TDirectory*)cutDir_p->mkdir("subDir");
  TDirectory* unfoldDir_p = (TDirectory*)cutDir_p->mkdir("unfoldDir");

  if(nHistDim != (int)histTag.size()) std::cout << "WARNING: nHistDim (" << nHistDim << ") != histTag.size() (" << histTag.size() << ")" << std::endl;
  if(nHistDim != (int)histBestBayes.size()) std::cout << "WARNING: nHistDim (" << nHistDim << ") != histBestBayes.size() (" << histBestBayes.size() << ")" << std::endl;

  cutPropData.SetNHistDim(nHistDim);
  cutPropData.SetHistTag(histTag);
  cutPropData.SetHistBestBayes(histBestBayes);

  if(!cutPropData.WriteAllVarToFile(outFile_p, cutDir_p, subDir_p, unfoldDir_p)) std::cout << "Warning: Cut writing has failed" << std::endl;

  outFile_p->Close();
  delete outFile_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 3 && argc != 4){
    std::cout << "Usage: ./bin/unfoldRawData.exe <inDataFileName> <inResponseName> <selectJtAlgo-Opt>" << std::endl;
    return 1;
  }

  int retVal = 0;

  if(argc == 3) retVal += unfoldRawData(argv[1], argv[2]);
  else if(argc == 4) retVal += unfoldRawData(argv[1], argv[2], argv[3]);

  std::cout << "Job complete. Return " << retVal << "." << std::endl;
  return retVal;
}
