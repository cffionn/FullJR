//cpp dependencies
#include <iostream>
#include <string>
#include <fstream>
#include <vector>

//ROOT dependencies
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TDatime.h"
#include "TDirectory.h"

//RooUnfold dependencies
#include "src/RooUnfoldResponse.h"
#include "src/RooUnfoldBayes.h"

//Local dependencies
#include "MainAnalysis/include/cutPropagator.h"

//Non-local FullJR dependencies (Utility, etc.)
#include "Utility/include/goodGlobalSelection.h"
#include "Utility/include/mntToXRootdFileString.h"
#include "Utility/include/plotUtilities.h"
#include "Utility/include/returnRootFileContentsList.h"

int unfoldRawData(const std::string inDataFileName, const std::string inResponseName)
{
  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  TFile* responseFile_p = new TFile(inResponseName.c_str(), "READ");
  std::vector<std::string> responseJetDirList = returnRootFileContentsList(responseFile_p, "TDirectoryFile", "JetAnalyzer");
  std::cout << "Printing " << responseJetDirList.size() << " response jets..." << std::endl;
  for(unsigned int jI = 0; jI < responseJetDirList.size(); ++jI){
    std::cout << " " << jI << "/" << responseJetDirList.size() << ": " << responseJetDirList.at(jI) << std::endl;
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

  const Int_t isDataPP = cutPropData.GetIsPP();
  const Int_t isResponsePP = cutPropResponse.GetIsPP();
  
  const Int_t nCentBins = cutPropData.GetNCentBins();
  std::vector<Int_t> centBinsLow = cutPropData.GetCentBinsLow();
  std::vector<Int_t> centBinsHi = cutPropData.GetCentBinsHi();

  const Int_t nJtPtBins = cutPropData.GetNJtPtBins();
  std::vector<Double_t> jtPtBinsTemp = cutPropData.GetJtPtBins();

  const Int_t nJtAbsEtaBins = cutPropData.GetNJtAbsEtaBins();
  std::vector<Double_t> jtAbsEtaBinsLowTemp = cutPropData.GetJtAbsEtaBinsLow();
  std::vector<Double_t> jtAbsEtaBinsHiTemp = cutPropData.GetJtAbsEtaBinsHi();

  const Int_t nID = cutPropData.GetNID();
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

  const Int_t nDataJet = dataJetDirList.size();

  std::string outFileName = inDataFileName;
  while(outFileName.find("/") != std::string::npos){outFileName.replace(0, outFileName.find("/")+1, "");}
  if(outFileName.find(".txt") != std::string::npos) outFileName.replace(outFileName.find(".txt"), std::string(".txt").size(), "");
  else if(outFileName.find(".root") != std::string::npos) outFileName.replace(outFileName.find(".root"), std::string(".root").size(), "");

  std::string outFileName2 = inResponseName;
  while(outFileName2.find("/") != std::string::npos){outFileName2.replace(0, outFileName2.find("/")+1, "");}
  if(outFileName2.find(".txt") != std::string::npos) outFileName2.replace(outFileName2.find(".txt"), std::string(".txt").size(), "");
  else if(outFileName2.find(".root") != std::string::npos) outFileName2.replace(outFileName2.find(".root"), std::string(".root").size(), "");

  
  while(outFileName.size() > 20){outFileName = outFileName.substr(0,outFileName.size()-1);}
  while(outFileName2.size() > 20){outFileName2 = outFileName2.substr(0,outFileName2.size()-1);}
  outFileName = "output/" + outFileName + "_" + outFileName2 + "_UnfoldRawData_" + dateStr + ".root";

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");

  const Int_t nBayes = 20;
  TDirectory* dir_p[nDataJet];
  TH1D* jtPtUnfolded_h[nDataJet][nCentBins][nID][nJtAbsEtaBins][nBayes];
  TH1D* jtPtUnfolded_RecoTrunc_h[nDataJet][nCentBins][nID][nJtAbsEtaBins][nBayes];

  for(Int_t jI = 0; jI < nDataJet; ++jI){
    std::string tempStr = dataJetDirList.at(jI);
    if(tempStr.find("/") != std::string::npos) tempStr.replace(tempStr.find("/"), tempStr.size() - tempStr.find("/"), "");
    dir_p[jI] = (TDirectory*)outFile_p->mkdir(tempStr.c_str());

    for(Int_t cI = 0; cI < nCentBins; ++cI){
      std::string centStr = "PP";
      if(!isDataPP) centStr = "Cent" + std::to_string(centBinsLow.at(cI)) + "to" + std::to_string(centBinsHi.at(cI));
      
      for(Int_t idI = 0; idI < nID; ++idI){
	for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
	  const std::string jtAbsEtaStr = "AbsEta" + prettyString(jtAbsEtaBinsLow[aI], 1, true) + "to" + prettyString(jtAbsEtaBinsHi[aI], 1, true);

	  for(Int_t bI = 0; bI < nBayes; ++bI){
	    std::string bayesStr = "Bayes" + std::to_string(bI+1);

	    jtPtUnfolded_h[jI][cI][idI][aI][bI] = new TH1D(("jtPtUnfolded_" + tempStr + "_" + centStr + "_" + idStr.at(idI) + "_" + jtAbsEtaStr + "_" + bayesStr + "_h").c_str(), ";Unfolded Jet p_{T};Counts", nJtPtBins, jtPtBins);
	    jtPtUnfolded_RecoTrunc_h[jI][cI][idI][aI][bI] = new TH1D(("jtPtUnfolded_RecoTrunc_" + tempStr + "_" + centStr + "_" + idStr.at(idI) + "_" + jtAbsEtaStr + "_" + bayesStr + "_h").c_str(), ";Unfolded Jet p_{T};Counts", nJtPtBins, jtPtBins);
	  }
	}
      }
    }
  }

  responseFile_p = new TFile(inResponseName.c_str(), "READ");
  RooUnfoldResponse* rooResponse_RecoTrunc_h[nDataJet][nCentBins][nID][nJtAbsEtaBins];

  for(Int_t jI = 0; jI < nDataJet; ++jI){
    std::string tempStr = dataJetDirList.at(jI);
    if(tempStr.find("/") != std::string::npos) tempStr.replace(tempStr.find("/"), tempStr.size() - tempStr.find("/"), "");

    for(Int_t cI = 0; cI < nCentBins; ++cI){
      std::string centStr = "PP";
      if(!isResponsePP) centStr = "Cent" + std::to_string(centBinsLow.at(cI)) + "to" + std::to_string(centBinsHi.at(cI));

      for(Int_t idI = 0; idI < nID; ++idI){
	for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
	  const std::string jtAbsEtaStr = "AbsEta" + prettyString(jtAbsEtaBinsLow[aI], 1, true) + "to" + prettyString(jtAbsEtaBinsHi[aI], 1, true);

	  rooResponse_RecoTrunc_h[jI][cI][idI][aI] = (RooUnfoldResponse*)responseFile_p->Get((tempStr + "/rooResponse_" + tempStr + "_" + centStr + "_" + idStr.at(idI) + "_" + jtAbsEtaStr + "_RecoTrunc_h").c_str());
	}
      }
    }
  }


  dataFile_p = new TFile(inDataFileName.c_str(), "READ");
  TH1D* jtPtRaw_RecoTrunc_h[nDataJet][nCentBins][nID][nJtAbsEtaBins];

  for(Int_t jI = 0; jI < nDataJet; ++jI){
    std::string tempStr = dataJetDirList.at(jI);
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

  for(Int_t jI = 0; jI < nDataJet; ++jI){
    for(Int_t cI = 0; cI < nCentBins; ++cI){
      for(Int_t idI = 0; idI < nID; ++idI){
	for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
	  for(Int_t bI = 0; bI < nBayes; ++bI){	    
	    RooUnfoldBayes bayes(rooResponse_RecoTrunc_h[jI][cI][idI][aI], jtPtRaw_RecoTrunc_h[jI][cI][idI][aI], 1+bI, false, "name");
	    bayes.SetVerbose(0);	    
	    TH1D* unfold_h = (TH1D*)bayes.Hreco(RooUnfold::kCovToy);	  

	    for(Int_t bIX = 0; bIX < unfold_h->GetNbinsX(); ++bIX){
	      jtPtUnfolded_RecoTrunc_h[jI][cI][idI][aI][bI]->SetBinContent(bIX+1, unfold_h->GetBinContent(bIX+1));
	      jtPtUnfolded_RecoTrunc_h[jI][cI][idI][aI][bI]->SetBinError(bIX+1, unfold_h->GetBinError(bIX+1));
	    }	    
	  }
	}
      }
    }
  }

  dataFile_p->Close();
  delete dataFile_p;

  responseFile_p->Close();
  delete responseFile_p;

  outFile_p->cd();

  for(Int_t jI = 0; jI < nDataJet; ++jI){
    outFile_p->cd();
    dir_p[jI]->cd();

    for(Int_t cI = 0; cI < nCentBins; ++cI){
      for(Int_t idI = 0; idI < nID; ++idI){

	for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){

	  for(Int_t bI = 0; bI < nBayes; ++bI){
	    jtPtUnfolded_h[jI][cI][idI][aI][bI]->Write("", TObject::kOverwrite);
	    delete jtPtUnfolded_h[jI][cI][idI][aI][bI];
	    
	    jtPtUnfolded_RecoTrunc_h[jI][cI][idI][aI][bI]->Write("", TObject::kOverwrite);
	    delete jtPtUnfolded_RecoTrunc_h[jI][cI][idI][aI][bI];
	  }
	}
      }
    }
  }

  outFile_p->cd();
  TDirectory* cutDir_p = (TDirectory*)outFile_p->mkdir("cutDir");
  if(!cutPropData.WriteAllVarToFile(outFile_p, cutDir_p)) std::cout << "Warning: Cut writing has failed" << std::endl;

  outFile_p->Close();
  delete outFile_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 3){
    std::cout << "Usage: ./bin/unfoldRawData.exe <inDataFileName> <inResponseName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += unfoldRawData(argv[1], argv[2]);

  std::cout << "Job complete. Return " << retVal << "." << std::endl;
  return retVal;
}
