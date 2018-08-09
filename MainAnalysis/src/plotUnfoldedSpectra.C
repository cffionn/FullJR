//cpp dependencies
#include <iostream>
#include <string>

//ROOT dependencies
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TDatime.h"

//Local dependencies
#include "MainAnalysis/include/cutPropagator.h"

//Non-local dependencies
#include "Utility/include/checkMakeDir.h"
#include "Utility/include/lumiAndTAAUtil.h"
#include "Utility/include/returnRootFileContentsList.h"


int plotUnfoldedSpectra(const std::string inFileNamePP, const std::string inFileNamePbPb)
{
  if(!checkFile(inFileNamePP) || !checkFile(inFileNamePbPb)){
    if(!checkFile(inFileNamePP)) std::cout << "inFileNamePP, \'" << inFileNamePP << "\', is invalid. return 1" << std::endl;
    if(!checkFile(inFileNamePbPb)) std::cout << "inFileNamePbPb, \'" << inFileNamePbPb << "\', is invalid. return 1" << std::endl;
    return 1;
  }

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  const std::string dirStr = "pdfDir/" + dateStr;
  checkMakeDir(dirStr);

  TFile* inFilePP_p = new TFile(inFileNamePP.c_str(), "READ");
  cutPropagator cutPropPP;
  cutPropPP.Clean();
  cutPropPP.GetAllVarFromFile(inFilePP_p);
  std::vector<std::string> jetPPList = returnRootFileContentsList(inFilePP_p, "TDirectoryFile", "JetAnalyzer");

  TFile* inFilePbPb_p = new TFile(inFileNamePbPb.c_str(), "READ");
  cutPropagator cutPropPbPb;
  cutPropPbPb.Clean();
  cutPropPbPb.GetAllVarFromFile(inFilePbPb_p);
  std::vector<std::string> jetPbPbList = returnRootFileContentsList(inFilePbPb_p, "TDirectoryFile", "JetAnalyzer");

  if(!cutPropPP.GetIsPP() || cutPropPbPb.GetIsPP() || !cutPropPP.CheckPropagatorsMatch(cutPropPbPb, true, false)){
    if(!cutPropPP.GetIsPP()) std::cout << "inFileNamePP \'" << inFileNamePP << "\' is not pp, return 1" << std::endl;
    if(cutPropPbPb.GetIsPP()) std::cout << "inFileNamePbPb \'" << inFileNamePbPb << "\' is not PbPb, return 1" << std::endl;
    if(!cutPropPP.CheckPropagatorsMatch(cutPropPbPb, true, false)) std::cout << "Cut propagators of pp and pbpb do not match. inspect. return 1" << std::endl;
  
    inFilePP_p->Close();
    delete inFilePP_p;
    
    inFilePbPb_p->Close();
    delete inFilePbPb_p;
    
    return 1;
  }

  const Int_t nAdditionalSyst = 4;
  const std::string additionalSyst[nAdditionalSyst] = {"LumiUp", "LumiDown", "TAAUp", "TAADown"};

  const Int_t nJtPP = jetPPList.size();
  const Int_t nJtPbPb = jetPbPbList.size();

  const Int_t nSystOrig = cutPropPbPb.GetNSyst();
  const Int_t nSyst = nSystOrig + nAdditionalSyst;
  std::vector<std::string> systStr = cutPropPbPb.GetSystStr();
  for(Int_t sI = 0; sI < nAdditionalSyst; ++sI){
    systStr.push_back(additionalSyst[sI]);
  }

  const Int_t nBayes = cutPropPbPb.GetNBayes();
  std::vector<int> bayesVal = cutPropPbPb.GetBayesVal();

  const Int_t nResponseMod = cutPropPbPb.GetNResponseMod();
  std::vector<double> responseMod = cutPropPbPb.GetResponseMod();

  const Int_t nCentBins = cutPropPbPb.GetNCentBins();
  std::vector<Int_t> centBinsLow = cutPropPbPb.GetCentBinsLow();
  std::vector<Int_t> centBinsHi = cutPropPbPb.GetCentBinsHi();

  //  const Int_t nJtPtBins = cutPropPbPb.GetNJtPtBins();
  //  std::vector<Double_t> jtPtBinsTemp = cutPropPbPb.GetJtPtBins();

  const Int_t nJtAbsEtaBins = cutPropPbPb.GetNJtAbsEtaBins();
  std::vector<Double_t> jtAbsEtaBinsLow = cutPropPbPb.GetJtAbsEtaBinsLow();
  std::vector<Double_t> jtAbsEtaBinsHi = cutPropPbPb.GetJtAbsEtaBinsHi();

  const Int_t nID = cutPropPbPb.GetNID();
  std::vector<std::string> idStr = cutPropPbPb.GetIdStr();
  std::vector<double> jtPfCHMFCutLow = cutPropPbPb.GetJtPfCHMFCutLow();
  std::vector<double> jtPfCHMFCutHi = cutPropPbPb.GetJtPfCHMFCutHi();
  std::vector<double> jtPfMUMFCutLow = cutPropPbPb.GetJtPfMUMFCutLow();
  std::vector<double> jtPfMUMFCutHi = cutPropPbPb.GetJtPfMUMFCutHi();

  std::vector<int> jetPPMatchedToPbPb;
  for(Int_t jI = 0; jI < nJtPbPb; ++jI){
    const Int_t rValPbPb = getRVal(jetPbPbList.at(jI));

    bool isMatched = false;
    for(Int_t jI2 = 0; jI2 < nJtPP; ++jI2){
      const Int_t rValPP = getRVal(jetPPList.at(jI2));
      
      if(rValPbPb == rValPP){
	jetPPMatchedToPbPb.push_back(jI2);
	isMatched = true;
	break;
      }
    }

    if(!isMatched){
      std::cout << "Warning: " << jetPbPbList.at(jI) << " has no match." << std::endl;
      std::cout << " Options: ";
      for(Int_t jI2 = 0; jI2 < nJtPP; ++jI2){
	std::cout << jetPPList.at(jI2) << ", ";
      }
      std::cout << std::endl;
    }

  }

  if(jetPPMatchedToPbPb.size() != jetPbPbList.size()){
    std::cout << "Size of pp algos, " << jetPPMatchedToPbPb.size() << ",  matched to jetPbPbList of different size, " << jetPbPbList.size() << ". return 1" << std::endl;
    
    inFilePP_p->Close();
    delete inFilePP_p;
    
    inFilePbPb_p->Close();
    delete inFilePbPb_p;    

    return 1;
  }

  std::cout << "Processing the following pairs of jets: " << std::endl;
  for(Int_t jI = 0; jI < nJtPbPb; ++jI){
    std::cout << " " << jI << "/" << nJtPbPb << ": " << jetPbPbList.at(jI) << ", " << jetPPList.at(jetPPMatchedToPbPb.at(jI)) << std::endl;
  }

  inFilePbPb_p->cd();
  TH1D* jtPtUnfolded_RecoTrunc_PbPb_h[nJtPbPb][nCentBins][nID][nResponseMod][nJtAbsEtaBins][nSyst][nBayes];

  inFilePP_p->cd();
  TH1D* jtPtUnfolded_RecoTrunc_PP_h[nJtPbPb][nID][nResponseMod][nJtAbsEtaBins][nSyst][nBayes];


  for(Int_t idI = 0; idI < nID; ++idI){
    for(Int_t mI = 0; mI < nResponseMod; ++mI){
      const std::string resStr = "ResponseMod" + prettyString(responseMod.at(mI), 2, true);

      for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
	const std::string jtAbsEtaStr = "AbsEta" + prettyString(jtAbsEtaBinsLow.at(aI), 1, true) + "to" + prettyString(jtAbsEtaBinsHi.at(aI), 1, true);	

	for(Int_t sI = 0; sI < nSyst; ++sI){
	  const std::string tempSystStr = systStr.at(sI) + "_";

	  for(Int_t bI = 0; bI < nBayes; ++bI){
	    const std::string bayesStr = "Bayes" + std::to_string(bayesVal.at(bI));
	    
	    for(Int_t tI = 0; tI < nJtPbPb; ++tI){
	      for(Int_t cI = 0; cI < nCentBins; ++cI){		
		const std::string centStr = "Cent" + std::to_string(centBinsLow.at(cI)) + "to" + std::to_string(centBinsHi.at(cI));

		const std::string name = jetPbPbList.at(tI) + "/jtPtUnfolded_RecoTrunc_" + jetPbPbList.at(tI) + "_" + centStr + "_" + idStr.at(idI) + "_" + resStr + "_" + jtAbsEtaStr + "_" + tempSystStr + bayesStr + "_h";

		if(sI < nSystOrig){
		  jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idI][mI][aI][sI][bI] = (TH1D*)inFilePbPb_p->Get(name.c_str());	    
		  jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idI][mI][aI][sI][bI]->GetName();
		}
		else{
		  jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idI][mI][aI][sI][bI] = (TH1D*)jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idI][mI][aI][0][bI]->Clone(name.c_str());	    
		}
	      }
	    }

	    for(Int_t tI = 0; tI < nJtPP; ++tI){
	      const std::string name = jetPPList.at(tI) + "/jtPtUnfolded_RecoTrunc_" + jetPPList.at(tI) + "_" + idStr.at(idI) + "_" + resStr + "_" + jtAbsEtaStr + "_" + tempSystStr + bayesStr + "_h";

	      if(sI < nSystOrig){
		jtPtUnfolded_RecoTrunc_PP_h[tI][idI][mI][aI][sI][bI] = (TH1D*)inFilePP_p->Get(name.c_str());
		jtPtUnfolded_RecoTrunc_PP_h[tI][idI][mI][aI][sI][bI]->GetName();
	      }
	      else{
		jtPtUnfolded_RecoTrunc_PP_h[tI][idI][mI][aI][sI][bI] = (TH1D*)jtPtUnfolded_RecoTrunc_PP_h[tI][idI][mI][aI][0][bI]->Clone(name.c_str());
	      }
	    }
	  }
	}
      }
    }
  }

  const Double_t lumiFactor = getLumiFactor();
  const Double_t nMBEvents = getNMBEvents();

  for(Int_t idI = 0; idI < nID; ++idI){
    for(Int_t mI = 0; mI < nResponseMod; ++mI){
      for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
	const Double_t etaBinWidth = TMath::Abs(jtAbsEtaBinsHi[aI] - jtAbsEtaBinsLow[aI]);

	for(Int_t sI = 0; sI < nSyst; ++sI){
	  for(Int_t bI = 0; bI < nBayes; ++bI){
	    for(Int_t tI = 0; tI < nJtPbPb; ++tI){
	      for(Int_t cI = 0; cI < nCentBins; ++cI){		
		const std::string centStr = "Cent" + std::to_string(centBinsLow.at(cI)) + "to" + std::to_string(centBinsHi.at(cI));
		const Double_t centBinWidth = TMath::Abs(centBinsHi[cI] - centBinsLow[cI])/100.;

		Double_t tempTAAFactor = getTAAScaleFactor(centStr);
		if(isStrSame(systStr[sI], "TAAUp")) tempTAAFactor += tempTAAFactor*getTAAScaleFactorUp(centStr);
		else if(isStrSame(systStr[sI], "TAADown")) tempTAAFactor -= tempTAAFactor*getTAAScaleFactorDown(centStr);

		Double_t totalPbPbFactor = tempTAAFactor*nMBEvents*2.*etaBinWidth*centBinWidth;

		jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idI][mI][aI][sI][bI]->Scale(1./totalPbPbFactor);
		divBinWidth(jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idI][mI][aI][sI][bI]);
	      }
	    }

	    for(Int_t tI = 0; tI < nJtPP; ++tI){
	      Double_t tempLumiFactor = lumiFactor;
	      if(isStrSame(systStr[sI], "LumiUp")) tempLumiFactor += tempLumiFactor*getLumiPercentError();
	      else if(isStrSame(systStr[sI], "LumiDown")) tempLumiFactor -= tempLumiFactor*getLumiPercentError();	      

	      Double_t totalPPFactor = tempLumiFactor*2.*etaBinWidth;
	      jtPtUnfolded_RecoTrunc_PP_h[tI][idI][mI][aI][sI][bI]->Scale(1./totalPPFactor);
	      divBinWidth(jtPtUnfolded_RecoTrunc_PP_h[tI][idI][mI][aI][sI][bI]);
	    }
	  }
	}
      }
    }
  }

  

  

  inFilePP_p->Close();
  delete inFilePP_p;

  inFilePbPb_p->Close();
  delete inFilePbPb_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 3){
    std::cout << "Usage: ./bin/plotUnfoldedSpectra.exe <inFileNamePP> <inFileNamePbPb>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += plotUnfoldedSpectra(argv[1], argv[2]);
  return retVal;
}
