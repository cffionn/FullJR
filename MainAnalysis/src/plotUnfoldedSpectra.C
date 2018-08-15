//cpp dependencies
#include <iostream>
#include <string>
#include <vector>

//ROOT dependencies
#include "TFile.h"
#include "TH1D.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TDatime.h"
#include "TMath.h"

//Local dependencies
#include "MainAnalysis/include/cutPropagator.h"
#include "MainAnalysis/include/doLocalDebug.h"

//Non-local dependencies
#include "Utility/include/checkMakeDir.h"
#include "Utility/include/cppWatch.h"
#include "Utility/include/doGlobalDebug.h"
#include "Utility/include/histDefUtility.h"
#include "Utility/include/kirchnerPalette.h"
#include "Utility/include/lumiAndTAAUtil.h"
#include "Utility/include/returnRootFileContentsList.h"

int plotUnfoldedSpectra(const std::string inFileNamePP, const std::string inFileNamePbPb)
{
  double total = 0;
  cppWatch toPPFile;
  cppWatch ppFile;
  cppWatch pbpbFile;
  cppWatch comparisonOfFiles;

  toPPFile.start();

  std::cout << "Initializing files..." << std::endl;

  if(!checkFile(inFileNamePP) || !checkFile(inFileNamePbPb)){
    if(!checkFile(inFileNamePP)) std::cout << "inFileNamePP, \'" << inFileNamePP << "\', is invalid. return 1" << std::endl;
    if(!checkFile(inFileNamePbPb)) std::cout << "inFileNamePbPb, \'" << inFileNamePbPb << "\', is invalid. return 1" << std::endl;
    return 1;
  }

  kirchnerPalette kPalette;

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  const std::string dirStr = "pdfDir/" + dateStr;
  checkMakeDir("pdfDir");
  checkMakeDir(dirStr);

  toPPFile.stop();
  total += toPPFile.total();
  std::cout << "To PP File time: " << toPPFile.total() << "/" << total << std::endl;
  ppFile.start();

  TFile* inFilePP_p = new TFile(inFileNamePP.c_str(), "READ");
  cutPropagator cutPropPP;
  cutPropPP.Clean();
  cutPropPP.GetAllVarFromFile(inFilePP_p);

  std::vector<std::string> histTagPP = cutPropPP.GetHistTag();
  std::vector<int> histBestBayesPP = cutPropPP.GetHistBestBayes();

  std::map<std::string, int> histTagMapPP;
  for(unsigned int i = 0; i < histTagPP.size(); ++i){
    histTagMapPP[histTagPP.at(i)] = histBestBayesPP.at(i);
  }


  std::vector<std::string> jetPPList = returnRootFileContentsList(inFilePP_p, "TDirectoryFile", "JetAnalyzer", 1);

  std::cout << "JetPPList: " << std::endl;
  for(unsigned int jI = 0; jI < jetPPList.size(); ++jI){
    std::cout << " " << jI << "/" << jetPPList.size() << ": " << jetPPList.at(jI) << std::endl;
  }

  ppFile.stop();
  total += ppFile.total();
  std::cout << "PP File time: " << ppFile.total() << "/" << total << std::endl;
  pbpbFile.start();

  TFile* inFilePbPb_p = new TFile(inFileNamePbPb.c_str(), "READ");
  cutPropagator cutPropPbPb;
  cutPropPbPb.Clean();
  cutPropPbPb.GetAllVarFromFile(inFilePbPb_p);
  std::vector<std::string> jetPbPbList = returnRootFileContentsList(inFilePbPb_p, "TDirectoryFile", "JetAnalyzer", 1);

  std::vector<std::string> histTagPbPb = cutPropPbPb.GetHistTag();
  std::vector<int> histBestBayesPbPb = cutPropPbPb.GetHistBestBayes();
  std::map<std::string, int> histTagMapPbPb;
  for(unsigned int i = 0; i < histTagPbPb.size(); ++i){
    histTagMapPbPb[histTagPbPb.at(i)] = histBestBayesPbPb.at(i);
  }

  std::cout << "JetPbPbList: " << std::endl;
  for(unsigned int jI = 0; jI < jetPbPbList.size(); ++jI){
    std::cout << " " << jI << "/" << jetPbPbList.size() << ": " << jetPbPbList.at(jI) << std::endl;
  }

  pbpbFile.stop();
  total += pbpbFile.total();
  std::cout << "PBPB File time: " << pbpbFile.total() << "/" << total << std::endl;
  comparisonOfFiles.start();

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

  comparisonOfFiles.stop();
  total += comparisonOfFiles.total();
  std::cout << "Comp File time: " << comparisonOfFiles.total() << "/" << total << std::endl;

  Int_t minForForLoop = 10000000;
  if(doLocalDebug || doGlobalDebug) minForForLoop = 1;

  std::cout << "Extracting cuts..." << std::endl;

  const Int_t nAdditionalSyst = 4;
  const std::string additionalSyst[nAdditionalSyst] = {"LumiUp", "LumiDown", "TAAUp", "TAADown"};

  const Int_t nJtPP = jetPPList.size();
  const Int_t nJtPbPb = jetPbPbList.size();

  const Int_t nSystOrig = cutPropPbPb.GetNSyst();
  //  const Int_t nSyst = TMath::Min(minForForLoop, nSystOrig + nAdditionalSyst);
  const Int_t nSyst = nSystOrig + nAdditionalSyst;
  std::vector<std::string> systStr = cutPropPbPb.GetSystStr();
  for(Int_t sI = 0; sI < nAdditionalSyst; ++sI){
    systStr.push_back(additionalSyst[sI]);
  }

  const Int_t nBayesCap = 10;
  const Int_t nBayes = TMath::Min(nBayesCap, cutPropPbPb.GetNBayes());
  std::vector<int> bayesVal = cutPropPbPb.GetBayesVal();

  const Int_t nResponseMod = TMath::Min(minForForLoop, cutPropPbPb.GetNResponseMod());
  //  const Int_t nResponseMod = cutPropPbPb.GetNResponseMod();
  std::vector<double> responseMod = cutPropPbPb.GetResponseMod();

  const Int_t nCentBins = cutPropPbPb.GetNCentBins();
  std::vector<Int_t> centBinsLow = cutPropPbPb.GetCentBinsLow();
  std::vector<Int_t> centBinsHi = cutPropPbPb.GetCentBinsHi();
  std::vector<Double_t> centBinsScalingFact;
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    centBinsScalingFact.push_back(TMath::Power(10., cI+1.));
  }

  //  const Int_t nJtPtBins = cutPropPbPb.GetNJtPtBins();
  //  std::vector<Double_t> jtPtBinsTemp = cutPropPbPb.GetJtPtBins();

  //  const Int_t nJtAbsEtaBins = TMath::Min(minForForLoop, cutPropPbPb.GetNJtAbsEtaBins());
  const Int_t nJtAbsEtaBins = cutPropPbPb.GetNJtAbsEtaBins();
  std::vector<Double_t> jtAbsEtaBinsLow = cutPropPbPb.GetJtAbsEtaBinsLow();
  std::vector<Double_t> jtAbsEtaBinsHi = cutPropPbPb.GetJtAbsEtaBinsHi();

  //  const Int_t nID = cutPropPbPb.GetNID();
  const Int_t nID = TMath::Min(minForForLoop, cutPropPbPb.GetNID());
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

  if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  inFilePbPb_p->cd();
  
  if(doLocalDebug || doGlobalDebug){
    std::cout << "nJtPbPb*nCentBins*nID*nResponseMod*nJtAbsEtaBins*nSyst*nBayes=Total" << std::endl;
    Int_t Total = nJtPbPb*nCentBins*nID*nResponseMod*nJtAbsEtaBins*nSyst*nBayes;
    std::cout << nJtPbPb << "*" << nCentBins << "*" << nID << "*" << nResponseMod << "*" << nJtAbsEtaBins << "*" << nSyst << "*" << nBayes << "=" << Total << std::endl;
  }

  TH1D* jtPtUnfolded_RecoTrunc_PbPb_h[nJtPbPb][nCentBins][nID][nResponseMod][nJtAbsEtaBins][nSyst][nBayes];

  if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  inFilePP_p->cd();
  TH1D* jtPtUnfolded_RecoTrunc_PP_h[nJtPP][nID][nResponseMod][nJtAbsEtaBins][nSyst][nBayes];

  if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  std::cout << "Grabbing histograms..." << std::endl;

  if(doLocalDebug || doGlobalDebug) std::cout << "idI,mI,aI,sI,bI,tI,cI" << std::endl;

  cppWatch histGrab;
  histGrab.start();

  inFilePbPb_p->cd();

  for(Int_t idI = 0; idI < nID; ++idI){
    std::cout << "Grabbing id " << idI << "/" << nID << ": " << idStr.at(idI) << std::endl;
    for(Int_t mI = 0; mI < nResponseMod; ++mI){
      const std::string resStr = "ResponseMod" + prettyString(responseMod.at(mI), 2, true);

      std::cout << " Grabbing mI " << mI << "/" << nResponseMod << ": " << resStr << std::endl;

      for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
	const std::string jtAbsEtaStr = "AbsEta" + prettyString(jtAbsEtaBinsLow.at(aI), 1, true) + "to" + prettyString(jtAbsEtaBinsHi.at(aI), 1, true);	

	std::cout << "  Grabbing aI " << aI << "/" << nJtAbsEtaBins << ": " << jtAbsEtaStr << std::endl;

	for(Int_t sI = 0; sI < nSyst; ++sI){
	  const std::string tempSystStr = systStr.at(sI) + "_";
	  
	  std::cout << "   Grabbing sI " << sI << "/" << nSyst << ": " << systStr.at(sI) << std::endl;

	  for(Int_t bI = 0; bI < nBayes; ++bI){
	    const std::string bayesStr = "Bayes" + std::to_string(bayesVal.at(bI));
	    
	    for(Int_t tI = 0; tI < nJtPbPb; ++tI){
	      for(Int_t cI = 0; cI < nCentBins; ++cI){		
		//		if(doLocalDebug || doGlobalDebug) std::cout << idI << "," << mI << "," << aI << "," << sI << "," << bI << "," << tI << "," << cI << std::endl;

		const std::string centStr = "Cent" + std::to_string(centBinsLow.at(cI)) + "to" + std::to_string(centBinsHi.at(cI));

		const std::string name = jetPbPbList.at(tI) + "/jtPtUnfolded_RecoTrunc_" + jetPbPbList.at(tI) + "_" + centStr + "_" + idStr.at(idI) + "_" + resStr + "_" + jtAbsEtaStr + "_" + tempSystStr + bayesStr + "_h";

		if(sI < nSystOrig){
		  jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idI][mI][aI][sI][bI] = (TH1D*)inFilePbPb_p->Get(name.c_str());	    
		  centerTitles(jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idI][mI][aI][sI][bI]);
		  if(jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idI][mI][aI][sI][bI]->GetSumw2()->GetSize() == 0) setSumW2(jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idI][mI][aI][sI][bI]);
		}
		else{
		  jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idI][mI][aI][sI][bI] = (TH1D*)jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idI][mI][aI][0][bI]->Clone(name.c_str());	    
		}
	      }
	    }
	  }
	}
      }
    }
  }


  inFilePP_p->cd();

  for(Int_t idI = 0; idI < nID; ++idI){
    std::cout << "Grabbing id " << idI << "/" << nID << ": " << idStr.at(idI) << std::endl;
    for(Int_t mI = 0; mI < nResponseMod; ++mI){
      const std::string resStr = "ResponseMod" + prettyString(responseMod.at(mI), 2, true);

      std::cout << " Grabbing mI " << mI << "/" << nResponseMod << ": " << resStr << std::endl;

      for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
	const std::string jtAbsEtaStr = "AbsEta" + prettyString(jtAbsEtaBinsLow.at(aI), 1, true) + "to" + prettyString(jtAbsEtaBinsHi.at(aI), 1, true);	

	std::cout << "  Grabbing aI " << aI << "/" << nJtAbsEtaBins << ": " << jtAbsEtaStr << std::endl;

	for(Int_t sI = 0; sI < nSyst; ++sI){
	  const std::string tempSystStr = systStr.at(sI) + "_";
	  
	  std::cout << "   Grabbing sI " << sI << "/" << nSyst << ": " << systStr.at(sI) << std::endl;

	  for(Int_t bI = 0; bI < nBayes; ++bI){
	    const std::string bayesStr = "Bayes" + std::to_string(bayesVal.at(bI));
	    
	    for(Int_t tI = 0; tI < nJtPP; ++tI){
	      const std::string name = jetPPList.at(tI) + "/jtPtUnfolded_RecoTrunc_" + jetPPList.at(tI) + "_PP_" + idStr.at(idI) + "_" + resStr + "_" + jtAbsEtaStr + "_" + tempSystStr + bayesStr + "_h";

	      if(sI < nSystOrig){
		jtPtUnfolded_RecoTrunc_PP_h[tI][idI][mI][aI][sI][bI] = (TH1D*)inFilePP_p->Get(name.c_str());
		centerTitles(jtPtUnfolded_RecoTrunc_PP_h[tI][idI][mI][aI][sI][bI]);
		if(jtPtUnfolded_RecoTrunc_PP_h[tI][idI][mI][aI][sI][bI]->GetSumw2()->GetSize() == 0) setSumW2(jtPtUnfolded_RecoTrunc_PP_h[tI][idI][mI][aI][sI][bI]);
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

  histGrab.stop();
  std::cout << "Total histGrab: " << histGrab.total() << std::endl;

  if(doLocalDebug || doGlobalDebug){
    std::cout << "Started return" << std::endl;
    return 1;
  }


  const Double_t lumiFactor = getLumiFactor();
  const Double_t nMBEvents = getNMBEvents();

  std::cout << "Scaling histograms..." << std::endl;

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

  std::cout << "Writing histograms..." << std::endl;

  std::string outFileName = inFileNamePP;
  while(outFileName.find("/") != std::string::npos){outFileName.replace(0, outFileName.find("/")+1, "");}
  if(outFileName.find(".txt") != std::string::npos) outFileName.replace(outFileName.find(".txt"), std::string(".txt").size(), "");
  else if(outFileName.find(".root") != std::string::npos) outFileName.replace(outFileName.find(".root"), std::string(".root").size(), "");

  std::string outFileName2 = inFileNamePbPb;
  while(outFileName2.find("/") != std::string::npos){outFileName2.replace(0, outFileName2.find("/")+1, "");}
  if(outFileName2.find(".txt") != std::string::npos) outFileName2.replace(outFileName2.find(".txt"), std::string(".txt").size(), "");
  else if(outFileName2.find(".root") != std::string::npos) outFileName2.replace(outFileName2.find(".root"), std::string(".root").size(), "");

  const Int_t sizeToTruncName = 40;
  while(outFileName.size() > sizeToTruncName){outFileName = outFileName.substr(0,outFileName.size()-1);}
  while(outFileName2.size() > sizeToTruncName){outFileName2 = outFileName2.substr(0,outFileName2.size()-1);}
  outFileName = "output/" + outFileName + "_" + outFileName2 + "_PlotUnfoldedSpectra_" + dateStr + ".root";

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  //Following two lines necessary else you get absolutely clobbered on deletion timing. See threads here:
  //https://root-forum.cern.ch/t/closing-root-files-is-dead-slow-is-this-a-normal-thing/5273/16
  //https://root-forum.cern.ch/t/tfile-speed/17549/25
  //Bizarre

  outFile_p->SetBit(TFile::kDevNull);
  TH1::AddDirectory(kFALSE);
 
  TDirectory* dirPP_p[nJtPP];
  TDirectory* dirPbPb_p[nJtPbPb];
  
  for(Int_t tI = 0; tI < nJtPbPb; ++tI){
    std::string tempStr = jetPbPbList.at(tI);
    if(tempStr.find("/") != std::string::npos) tempStr.replace(tempStr.find("/"), tempStr.size() - tempStr.find("/"), "");
    dirPbPb_p[tI] = (TDirectory*)outFile_p->mkdir(tempStr.c_str());
    dirPbPb_p[tI]->cd();

    for(Int_t cI = 0; cI < nCentBins; ++cI){		
      for(Int_t idI = 0; idI < nID; ++idI){
	for(Int_t mI = 0; mI < nResponseMod; ++mI){
	  for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
	    for(Int_t sI = 0; sI < nSyst; ++sI){
	      for(Int_t bI = 0; bI < nBayes; ++bI){
		std::string newName = jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idI][mI][aI][sI][bI]->GetName();
		newName.replace(newName.find("_h"), 2, "_Rescaled_h");
		jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idI][mI][aI][sI][bI]->Write(newName.c_str(), TObject::kOverwrite);
	      }
	    }
	  }
	}
      }
    }
  }

  

  for(Int_t tI = 0; tI < nJtPP; ++tI){
    std::string tempStr = jetPPList.at(tI);
    if(tempStr.find("/") != std::string::npos) tempStr.replace(tempStr.find("/"), tempStr.size() - tempStr.find("/"), "");
    dirPP_p[tI] = (TDirectory*)outFile_p->mkdir(tempStr.c_str());
    dirPP_p[tI]->cd();

    for(Int_t idI = 0; idI < nID; ++idI){
      for(Int_t mI = 0; mI < nResponseMod; ++mI){
	for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
	  for(Int_t sI = 0; sI < nSyst; ++sI){
	    for(Int_t bI = 0; bI < nBayes; ++bI){
	      std::string newName = jtPtUnfolded_RecoTrunc_PP_h[tI][idI][mI][aI][sI][bI]->GetName();
	      newName.replace(newName.find("_h"), 2, "_Rescaled_h");
	      jtPtUnfolded_RecoTrunc_PP_h[tI][idI][mI][aI][sI][bI]->Write(newName.c_str(), TObject::kOverwrite);
	    }
	  }
	}
      }
    }
  }

  std::cout << "Closing files..." << std::endl;

  outFile_p->Close();
  delete outFile_p;

  
  std::cout << "Plotting... " << std::endl;
  const std::string plotID = "LightMUAndCHID";
  const std::string plotAbsEtaStr = "AbsEta0p0to2p0";
  const Int_t plotBayesVal = 4;

  Int_t idPos = -1;
  Int_t absEtaPos = -1;
  Int_t bayesPos = -1;

  for(Int_t idI = 0; idI < nID; ++idI){
    if(isStrSame(idStr.at(idI), plotID)){
      idPos = idI;
      break;
    }
  }

  for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
    const std::string jtAbsEtaStr = "AbsEta" + prettyString(jtAbsEtaBinsLow.at(aI), 1, true) + "to" + prettyString(jtAbsEtaBinsHi.at(aI), 1, true);	
    
    if(isStrSame(jtAbsEtaStr, plotAbsEtaStr)){
      absEtaPos = aI;
      break;
    }
  }

  for(Int_t bI = 0; bI < nBayes; ++bI){
    if(bayesVal.at(bI) == plotBayesVal){
      bayesPos = bI;
      break;
    }
  }
  
  if(idPos == -1 || absEtaPos == -1 || bayesPos == -1){
    std::cout << "Plotting pos not found, skip plotting, return 1" << std::endl;

    inFilePP_p->Close();
    delete inFilePP_p;

    inFilePbPb_p->Close();
    delete inFilePbPb_p;

    return 1;
  }

  const std::string idNameStr = idStr.at(idPos);
  const std::string plotBayesStr = "Bayes" + std::to_string(plotBayesVal);

  for(Int_t tI = 0; tI < nJtPbPb; ++tI){
    const Int_t rVal = getRVal(jetPbPbList.at(tI));
    const Int_t ppPos = jetPPMatchedToPbPb.at(tI);
    const std::string rValStr = std::to_string(rVal);

    for(Int_t mI = 0; mI < nResponseMod; ++mI){
      const std::string responseStr = prettyString(responseMod.at(mI), 2, true);

      TCanvas* spectCanv_p = new TCanvas("spectCanv_p", "", 450, 450);
      spectCanv_p->SetTopMargin(0.01);
      spectCanv_p->SetRightMargin(0.01);
      spectCanv_p->SetLeftMargin(0.12);
      spectCanv_p->SetBottomMargin(0.12);
      
      Double_t minVal = getMinGTZero(jtPtUnfolded_RecoTrunc_PP_h[ppPos][idPos][mI][absEtaPos][0][bayesPos]);
      Double_t maxVal = getMax(jtPtUnfolded_RecoTrunc_PP_h[ppPos][idPos][mI][absEtaPos][0][bayesPos]);

      for(Int_t cI = 0; cI < nCentBins; ++cI){
	scaleCentralAndErrorValues(jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idPos][mI][absEtaPos][0][bayesPos], centBinsScalingFact.at(cI));

	Double_t tempMinVal = getMinGTZero(jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idPos][mI][absEtaPos][0][bayesPos]);
	Double_t tempMaxVal = getMax(jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idPos][mI][absEtaPos][0][bayesPos]);
	if(minVal > tempMinVal) minVal = tempMinVal;
	if(maxVal < tempMaxVal) maxVal = tempMaxVal;
      }

      jtPtUnfolded_RecoTrunc_PP_h[ppPos][idPos][mI][absEtaPos][0][bayesPos]->SetMaximum(maxVal*10.);
      jtPtUnfolded_RecoTrunc_PP_h[ppPos][idPos][mI][absEtaPos][0][bayesPos]->SetMinimum(minVal/10.);

      jtPtUnfolded_RecoTrunc_PP_h[ppPos][idPos][mI][absEtaPos][0][bayesPos]->SetMarkerColor(kPalette.getColor(getColorPosFromCent("", true)));
      jtPtUnfolded_RecoTrunc_PP_h[ppPos][idPos][mI][absEtaPos][0][bayesPos]->SetLineColor(kPalette.getColor(getColorPosFromCent("", true)));
      jtPtUnfolded_RecoTrunc_PP_h[ppPos][idPos][mI][absEtaPos][0][bayesPos]->SetMarkerStyle(getStyleFromCent("", true));
      jtPtUnfolded_RecoTrunc_PP_h[ppPos][idPos][mI][absEtaPos][0][bayesPos]->SetMarkerSize(1);

      jtPtUnfolded_RecoTrunc_PP_h[ppPos][idPos][mI][absEtaPos][0][bayesPos]->DrawCopy("HIST E1 P");

      for(Int_t cI = 0; cI < nCentBins; ++cI){
	const std::string centStr = "Cent" + std::to_string(centBinsLow.at(cI)) + "to" + std::to_string(centBinsHi.at(cI));

	jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idPos][mI][absEtaPos][0][bayesPos]->SetMarkerColor(kPalette.getColor(getColorPosFromCent(centStr, false)));
	jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idPos][mI][absEtaPos][0][bayesPos]->SetLineColor(kPalette.getColor(getColorPosFromCent(centStr, false)));
	jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idPos][mI][absEtaPos][0][bayesPos]->SetMarkerStyle(getStyleFromCent(centStr, false));
	jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idPos][mI][absEtaPos][0][bayesPos]->SetMarkerSize(1);

	jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idPos][mI][absEtaPos][0][bayesPos]->DrawCopy("HIST E1 P SAME");
      }

      gPad->SetLogy();
      gPad->SetLogx();
      gStyle->SetOptStat(0);

      std::string saveName = "spectra_" + jetPbPbList.at(tI) + "_R" + rValStr + "_" + idNameStr + "_" + plotAbsEtaStr + "_" + plotBayesStr + "_" + responseStr + "_" + dateStr + ".pdf";
      saveName = "pdfDir/" + dateStr + "/" + saveName;

      spectCanv_p->SaveAs(saveName.c_str());

      delete spectCanv_p;

      for(Int_t cI = 0; cI < nCentBins; ++cI){
	scaleCentralAndErrorValues(jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idPos][mI][absEtaPos][0][bayesPos], 1./centBinsScalingFact.at(cI));
      }
    }
  }


  inFilePP_p->Close();
  delete inFilePP_p;

  inFilePbPb_p->Close();
  delete inFilePbPb_p;

  std::cout << "Job complete" << std::endl;

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
