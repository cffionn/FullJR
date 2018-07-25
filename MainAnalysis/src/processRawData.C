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
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLegend.h"

//Local dependencies
#include "MainAnalysis/include/cutPropagator.h"

//Non-local FullJR dependencies (Utility, etc.)
#include "Utility/include/goodGlobalSelection.h"
#include "Utility/include/histDefUtility.h"
#include "Utility/include/mntToXRootdFileString.h"
#include "Utility/include/plotUtilities.h"
#include "Utility/include/returnRootFileContentsList.h"
#include "Utility/include/vanGoghPalette.h"


int processRawData(const std::string inDataFileName, const std::string inResponseName, bool isDataPP = false)
{
  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  vanGoghPalette vg;
  const Int_t nStyles = 3;
  const Int_t styles[nStyles] = {20, 21, 34};
  const Int_t colors[nStyles] = {1, vg.getColor(0), vg.getColor(1)};
  const Double_t yPadFrac = 0.35;

  TFile* responseFile_p = new TFile(inResponseName.c_str(), "READ");
  std::vector<std::string> jetDirList = returnRootFileContentsList(responseFile_p, "TDirectoryFile", "JetAnalyzer");

  std::cout << "Validating " << jetDirList.size() << " jets..." << std::endl;
  for(unsigned int jI = 0; jI < jetDirList.size(); ++jI){
    std::cout << " " << jI << "/" << jetDirList.size() << ": " << jetDirList.at(jI) << std::endl;
  }

  cutPropagator cutProp;
  cutProp.Clean();
  cutProp.GetAllVarFromFile(responseFile_p);
  Int_t nCentBinsTemp = cutProp.GetNCentBins();
  std::vector<Int_t> centBinsLow = cutProp.GetCentBinsLow();
  std::vector<Int_t> centBinsHi = cutProp.GetCentBinsHi();

  Float_t jtAbsEtaMaxTemp = cutProp.GetJtAbsEtaMax();

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

  if(nCentBinsTemp < 0) std::cout << "nCentBins less than 0. please check input file. return 1" << std::endl;
  if(nIDTemp < 0) std::cout << "nID less than 0. please check input file. return 1" << std::endl;
  if(nJtPtBinsTemp < 0) std::cout << "nJtPtBinsTemp less than 0. please check input file. return 1" << std::endl;
  if(nJtAbsEtaBinsTemp < 0) std::cout << "nJtAbsEtaBinsTemp less than 0. please check input file. return 1" << std::endl;
  if(jtAbsEtaMaxTemp < -98) std::cout << "jtAbsEtaMaxTemp less than -98. please check input file. return 1" << std::endl;

  if(nCentBinsTemp < 0 || nJtPtBinsTemp < 0 || nJtAbsEtaBinsTemp < 0 || jtAbsEtaMaxTemp < -98 || nIDTemp < 0){
    responseFile_p->Close();
    delete responseFile_p;
    return 1;
  }

  const Float_t jtAbsEtaMax = jtAbsEtaMaxTemp;
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

  std::cout << "Raw data from file: \'" << inDataFileName << "\'" << std::endl;

  std::vector<std::string> fileList;

  if(inDataFileName.find(".root") != std::string::npos){
    fileList.push_back(inDataFileName);
  }
  else if(inDataFileName.find(".txt") != std::string::npos){
    std::ifstream file(inDataFileName.c_str());
    std::string tempStr;

    while(std::getline(file, tempStr)){
      while(tempStr.find(" ") != std::string::npos){tempStr.replace(tempStr.find(" "), 1, "");}
      if(tempStr.size() == 0) continue;
      if(tempStr.find(".root") != std::string::npos) fileList.push_back(tempStr);
      else{
        if(tempStr.substr(0, std::string("ISPP=").size()).find("ISPP=") != std::string::npos){
          tempStr.replace(0,tempStr.find("=")+1, "");
          isDataPP = std::stoi(tempStr);
        }
        else std::cout << "WARNING: Line in \'" << inDataFileName << "\', \'" << tempStr << "\' is invalid. check input" << std::endl;
      }
    }

    file.close();
  }
  else{
    std::cout << "Given inDataFileName \'" << inDataFileName << "\' is invalid. return 1" << std::endl;
    return 1;
  }

  if(fileList.size() == 0){
    std::cout << "Given inDataFileName \'" << inDataFileName << "\' is gives no valid root files. return 1" << std::endl;
    return 1;
  }
  else if(isDataPP != isResponsePP){
    std::cout << "Response isDataPP == " << isResponsePP << ", data isDataPP == " << isDataPP << " aren't equivalent. return 1" << std::endl;
    return 1;
  }

  std::cout << "Checking for matched inputs..." << std::endl;
  TFile* inDataFile_p = TFile::Open(mntToXRootdFileString(fileList.at(0)).c_str(), "READ");
  std::vector<std::string> dataTreeList = returnRootFileContentsList(inDataFile_p, "TTree", "JetAna");

  unsigned int pos = 0;
  while(dataTreeList.size() > pos){   
    bool isFound = false;
    std::string tempJet = dataTreeList.at(pos);
    tempJet.replace(tempJet.find("/"), tempJet.size() - tempJet.find("/"), "");

    for(unsigned int jI = 0; jI < jetDirList.size(); ++jI){
      if(jetDirList.at(jI).find(tempJet) != std::string::npos && tempJet.size() == jetDirList.at(jI).size()){
	isFound = true;
	break;
      }
    }

    if(isFound) ++pos;
    else dataTreeList.erase(dataTreeList.begin() + pos);
  }

  inDataFile_p->Close();
  delete inDataFile_p;
  inDataFile_p = NULL;

  const Int_t nDataJet = dataTreeList.size();

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
  outFileName = "output/" + outFileName + "_" + outFileName2 + "_ProcessRawData_" + dateStr + ".root";

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  TDirectory* dir_p[nDataJet];
  TH1D* jtPtRaw_h[nDataJet][nCentBins][nID][nJtAbsEtaBins];
  TH1D* jtPtRaw_RecoTrunc_h[nDataJet][nCentBins][nID][nJtAbsEtaBins];

  for(Int_t jI = 0; jI < nDataJet; ++jI){
    std::string tempStr = dataTreeList.at(jI);
    tempStr.replace(tempStr.find("/"), tempStr.size() - tempStr.find("/"), "");
    dir_p[jI] = (TDirectory*)outFile_p->mkdir(tempStr.c_str());
    std::cout << " " << jI << "/" << nDataJet << ": " << tempStr << std::endl;

    for(Int_t cI = 0; cI < nCentBins; ++cI){
      std::string centStr = "PP";
      if(!isDataPP) centStr = "Cent" + std::to_string(centBinsLow.at(cI)) + "to" + std::to_string(centBinsHi.at(cI));
      
      for(Int_t idI = 0; idI < nID; ++idI){
	for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
	  const std::string jtAbsEtaStr = "AbsEta" + prettyString(jtAbsEtaBinsLow[aI], 1, true) + "to" + prettyString(jtAbsEtaBinsHi[aI], 1, true);
	  jtPtRaw_h[jI][cI][idI][aI] = new TH1D(("jtPtRaw_" + tempStr + "_" + centStr + "_" + idStr.at(idI) + "_" + jtAbsEtaStr + "_h").c_str(), ";Raw Jet p_{T};Counts", nJtPtBins, jtPtBins);
	  jtPtRaw_RecoTrunc_h[jI][cI][idI][aI] = new TH1D(("jtPtRaw_RecoTrunc_" + tempStr + "_" + centStr + "_" + idStr.at(idI) + "_" + jtAbsEtaStr + "_h").c_str(), ";Raw Jet p_{T};Counts", nJtPtBins, jtPtBins);

	  setSumW2({jtPtRaw_h[jI][cI][idI][aI], jtPtRaw_RecoTrunc_h[jI][cI][idI][aI]});
	  centerTitles({jtPtRaw_h[jI][cI][idI][aI], jtPtRaw_RecoTrunc_h[jI][cI][idI][aI]});
	}
      }
    }
  }

  std::cout << "Processing " << fileList.size() << " files..., isDataPP == " << isDataPP << std::endl;
  std::cout << "Jets: ";
  for(unsigned int jI = 0; jI < dataTreeList.size(); ++jI){
    std::cout << dataTreeList.at(jI) << ", ";
  }
  std::cout << std::endl;


  //Jet Var  
  const Int_t nMaxJet = 500;
  Int_t nref_[nDataJet];
  Float_t jtpt_[nDataJet][nMaxJet];
  Float_t jteta_[nDataJet][nMaxJet];
  Float_t jtphi_[nDataJet][nMaxJet];
  Float_t jtPfCHMF_[nDataJet][nMaxJet];
  Float_t jtPfMUMF_[nDataJet][nMaxJet];

  Float_t vz_;
  Float_t hiHF_;
  Int_t hiBin_;
  unsigned int run_, lumi_;
  unsigned long long evt_;

  Int_t HBHENoiseFilterResultRun2Loose_ = -1;
  Int_t pprimaryVertexFilter_ = -1;
  Int_t pBeamScrapingFilter_ = -1;
  Int_t phfCoincFilter3_ = -1;
  Int_t pclusterCompatibilityFilter_ = -1;

  goodGlobalSelection globalSel;
  globalSel.setIsPbPb(!isDataPP);

  for(unsigned int fI = 0; fI < fileList.size(); ++fI){
    std::cout << "File " << fI << "/" << fileList.size() << ": " << fileList.at(fI) << std::endl;

    TFile* inDataFile_p = TFile::Open(mntToXRootdFileString(fileList.at(fI)).c_str(), "READ");
    TTree* jetTrees_p[nDataJet];
    TTree* hiTree_p = (TTree*)inDataFile_p->Get("hiEvtAnalyzer/HiTree");
    TTree* skimTree_p = (TTree*)inDataFile_p->Get("skimanalysis/HltTree");

    for(Int_t jI = 0; jI < nDataJet; ++jI){
      jetTrees_p[jI] = NULL;
      jetTrees_p[jI] = (TTree*)inDataFile_p->Get(dataTreeList.at(jI).c_str());

      jetTrees_p[jI]->SetBranchStatus("*", 0);
      jetTrees_p[jI]->SetBranchStatus("nref", 1);
      jetTrees_p[jI]->SetBranchStatus("jtpt", 1);
      jetTrees_p[jI]->SetBranchStatus("jtphi", 1);
      jetTrees_p[jI]->SetBranchStatus("jteta", 1);
      jetTrees_p[jI]->SetBranchStatus("jtPfCHMF", 1);
      jetTrees_p[jI]->SetBranchStatus("jtPfMUMF", 1);

      jetTrees_p[jI]->SetBranchAddress("nref", &(nref_[jI]));
      jetTrees_p[jI]->SetBranchAddress("jtpt", jtpt_[jI]);
      jetTrees_p[jI]->SetBranchAddress("jtphi", jtphi_[jI]);
      jetTrees_p[jI]->SetBranchAddress("jteta", jteta_[jI]);
      jetTrees_p[jI]->SetBranchAddress("jtPfCHMF", jtPfCHMF_[jI]);
      jetTrees_p[jI]->SetBranchAddress("jtPfMUMF", jtPfMUMF_[jI]);
    }
    
    hiTree_p->SetBranchStatus("*", 0);
    hiTree_p->SetBranchStatus("hiBin", 1);
    hiTree_p->SetBranchStatus("vz", 1);
    hiTree_p->SetBranchStatus("hiHF", 1);
    hiTree_p->SetBranchStatus("run", 1);
    hiTree_p->SetBranchStatus("lumi", 1);
    hiTree_p->SetBranchStatus("evt", 1);
    
    hiTree_p->SetBranchAddress("hiBin", &hiBin_);
    hiTree_p->SetBranchAddress("vz", &vz_);
    hiTree_p->SetBranchAddress("hiHF", &hiHF_);
    hiTree_p->SetBranchAddress("run", &run_);
    hiTree_p->SetBranchAddress("lumi", &lumi_);
    hiTree_p->SetBranchAddress("evt", &evt_);
    
    
    skimTree_p->SetBranchStatus("*", 0);
    skimTree_p->SetBranchStatus("HBHENoiseFilterResultRun2Loose", 1);
    skimTree_p->SetBranchAddress("HBHENoiseFilterResultRun2Loose", &HBHENoiseFilterResultRun2Loose_);
    
    if(!isDataPP){
      skimTree_p->SetBranchStatus("pprimaryVertexFilter", 1);
      skimTree_p->SetBranchStatus("phfCoincFilter3", 1);
      skimTree_p->SetBranchStatus("pclusterCompatibilityFilter", 1);
      
      skimTree_p->SetBranchAddress("pprimaryVertexFilter", &pprimaryVertexFilter_);
      skimTree_p->SetBranchAddress("phfCoincFilter3", &phfCoincFilter3_);
      skimTree_p->SetBranchAddress("pclusterCompatibilityFilter", &pclusterCompatibilityFilter_);
    }
    else{
      skimTree_p->SetBranchStatus("pBeamScrapingFilter", 1);
      skimTree_p->SetBranchStatus("pPAprimaryVertexFilter", 1);
      
      skimTree_p->SetBranchAddress("pBeamScrapingFilter", &pBeamScrapingFilter_);
      skimTree_p->SetBranchAddress("pPAprimaryVertexFilter", &pprimaryVertexFilter_);
    }
    

    const Int_t nEntries = hiTree_p->GetEntries();
    const Int_t printInterval = TMath::Max(1, nEntries/20);

    for(Int_t entry = 0; entry < nEntries; ++entry){
      if(nEntries >= 50000 && entry%printInterval == 0) std::cout << " Entry: " << entry << "/" << nEntries << std::endl;

      hiTree_p->GetEntry(entry);
      skimTree_p->GetEntry(entry);

      globalSel.setVz(vz_);
      globalSel.setHiHF(hiHF_);
      globalSel.setPprimaryVertexFilter(pprimaryVertexFilter_);
      globalSel.setPBeamScrapingFilter(pBeamScrapingFilter_);
      globalSel.setPhfCoincFilter3(phfCoincFilter3_);
      globalSel.setHBHENoiseFilterResultRun2Loose(HBHENoiseFilterResultRun2Loose_);
      globalSel.setPclusterCompatibilityFilter(pclusterCompatibilityFilter_);

      if(!globalSel.isGood()) continue;

      Int_t centPos = -1;
      for(Int_t cI = 0; cI < nCentBins; ++cI){
	if(hiBin_/2 >= centBinsLow[cI] && hiBin_/2 < centBinsHi[cI]){
	  centPos = cI;
	  break;
	}
      }

      if(centPos == -1) continue;

      for(Int_t tI = 0; tI < nDataJet; ++tI){jetTrees_p[tI]->GetEntry(entry);}
      

      for(Int_t tI = 0; tI < nDataJet; ++tI){

	for(Int_t jI = 0; jI < nref_[tI]; ++jI){
	  if(TMath::Abs(jteta_[tI][jI]) > jtAbsEtaMax) continue;
	  
	  std::vector<int> jtAbsEtaPoses;
	  for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
	    if(TMath::Abs(jteta_[tI][jI]) >= jtAbsEtaBinsLow[aI] && TMath::Abs(jteta_[tI][jI]) < jtAbsEtaBinsHi[aI]){
	      jtAbsEtaPoses.push_back(aI);
	    }
	  }

	  /*
	  std::cout << "JtAbsEtaPoses.size(): " << jtAbsEtaPoses.size() << std::endl;
	  if(jtAbsEtaPoses.size() != 0) std::cout << " " << jtAbsEtaPoses.at(0) << ", " << jtAbsEtaPoses.at(jtAbsEtaPoses.size()-1) << std::endl;
	  */

	  std::vector<int> idPoses;
	  for(Int_t idI = 0; idI < nID; ++idI){
	    if(jtPfCHMFCutLow.at(idI) > jtPfCHMF_[tI][jI]) continue;
	    if(jtPfCHMFCutHi.at(idI) <= jtPfCHMF_[tI][jI]) continue;
	    if(jtPfMUMFCutLow.at(idI) > jtPfMUMF_[tI][jI]) continue;
	    if(jtPfMUMFCutHi.at(idI) <= jtPfMUMF_[tI][jI]) continue;

	    idPoses.push_back(idI);
	  }

	  for(unsigned int aI = 0; aI < jtAbsEtaPoses.size(); ++aI){
	    for(unsigned int idI = 0; idI < idPoses.size(); ++idI){
	      bool goodReco = (jtpt_[tI][jI] >= jtPtBins[0] && jtpt_[tI][jI] < jtPtBins[nJtPtBins]);
	      bool goodRecoTrunc = (jtpt_[tI][jI] >= jtPtBins[1] && jtpt_[tI][jI] < jtPtBins[nJtPtBins-1]);
	      
	      if(goodReco){
		jtPtRaw_h[tI][centPos][idPoses.at(idI)][jtAbsEtaPoses.at(aI)]->Fill(jtpt_[tI][jI]);
		if(goodRecoTrunc) jtPtRaw_RecoTrunc_h[tI][centPos][idPoses.at(idI)][jtAbsEtaPoses.at(aI)]->Fill(jtpt_[tI][jI]);
	      }
	    }
	  }

	}
      }
    }

    inDataFile_p->Close();
    delete inDataFile_p;
  }


  outFile_p->cd();

  for(Int_t jI = 0; jI < nDataJet; ++jI){
    std::string tempStr = dataTreeList.at(jI);
    tempStr.replace(tempStr.find("/"), tempStr.size() - tempStr.find("/"), "");

    for(Int_t cI = 0; cI < nCentBins; ++cI){
      std::string centStr = "PP";
      std::string centStr2 = "PP";
      if(!isDataPP){
	centStr = "Cent" + std::to_string(centBinsLow.at(cI)) + "to" + std::to_string(centBinsHi.at(cI));
	centStr2 =  std::to_string(centBinsLow.at(cI)) + "-" + std::to_string(centBinsHi.at(cI)) + "%";
      }

      for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
	const std::string jtAbsEtaStr = "AbsEta" + prettyString(jtAbsEtaBinsLow[aI], 1, true) + "to" + prettyString(jtAbsEtaBinsHi[aI], 1, true);
	const std::string jtAbsEtaStr2 = prettyString(jtAbsEtaBinsLow[aI], 1, true) + "<|#eta|<" + prettyString(jtAbsEtaBinsHi[aI], 1, true);
	
	TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);
	canv_p->SetTopMargin(0.01);
	canv_p->SetRightMargin(0.01);
	canv_p->SetLeftMargin(0.01);
	canv_p->SetBottomMargin(0.01);
	
	TPad* pads[2];
	canv_p->cd();
	pads[0] = new TPad("pad0", "", 0.0, yPadFrac, 1.0, 1.0);
	pads[0]->SetLeftMargin(0.12);
	pads[0]->SetTopMargin(0.01);
	pads[0]->SetBottomMargin(0.001);
	pads[0]->SetRightMargin(0.01);
	pads[0]->Draw();

	canv_p->cd();
	pads[1] = new TPad("pad1", "", 0.0, 0.0, 1.0, yPadFrac);
	pads[1]->Draw();
	pads[1]->SetLeftMargin(0.12);
	pads[1]->SetTopMargin(0.001);
	pads[1]->SetBottomMargin(0.12/yPadFrac);
	pads[1]->SetRightMargin(0.01);

	canv_p->cd();
	pads[0]->cd();

	double max = -1;
	double min = 10000000;

	double maxRat = -1;
	double minRat = 10000000;

	for(Int_t idI = 0; idI < nID; ++idI){
	  jtPtRaw_h[jI][cI][idI][aI]->SetMarkerColor(colors[idI]);
	  jtPtRaw_h[jI][cI][idI][aI]->SetLineColor(colors[idI]);
	  jtPtRaw_h[jI][cI][idI][aI]->SetMarkerStyle(styles[idI]);
	  jtPtRaw_h[jI][cI][idI][aI]->SetMarkerSize(1.);	  

	  jtPtRaw_h[jI][cI][idI][aI]->GetXaxis()->SetTitleFont(43);
	  jtPtRaw_h[jI][cI][idI][aI]->GetXaxis()->SetTitleSize(14);
	  jtPtRaw_h[jI][cI][idI][aI]->GetYaxis()->SetTitleFont(43);
	  jtPtRaw_h[jI][cI][idI][aI]->GetYaxis()->SetTitleSize(14);

	  jtPtRaw_h[jI][cI][idI][aI]->GetXaxis()->SetLabelFont(43);
	  jtPtRaw_h[jI][cI][idI][aI]->GetXaxis()->SetLabelSize(14);
	  jtPtRaw_h[jI][cI][idI][aI]->GetYaxis()->SetLabelFont(43);
	  jtPtRaw_h[jI][cI][idI][aI]->GetYaxis()->SetLabelSize(14);

	  jtPtRaw_h[jI][cI][idI][aI]->GetXaxis()->SetNdivisions(505);
	  jtPtRaw_h[jI][cI][idI][aI]->GetYaxis()->SetNdivisions(505);

	  jtPtRaw_h[jI][cI][idI][aI]->GetXaxis()->SetTitleOffset(3.0);
	  jtPtRaw_h[jI][cI][idI][aI]->GetYaxis()->SetTitleOffset(2.0);

	  for(Int_t bIX = 0; bIX < jtPtRaw_h[jI][cI][idI][aI]->GetNbinsX(); ++bIX){
	    if(jtPtRaw_h[jI][cI][idI][aI]->GetBinContent(bIX+1) > max) max = jtPtRaw_h[jI][cI][idI][aI]->GetBinContent(bIX+1);
	    if(jtPtRaw_h[jI][cI][idI][aI]->GetBinContent(bIX+1) < min && jtPtRaw_h[jI][cI][idI][aI]->GetBinContent(bIX+1) > 0) min = jtPtRaw_h[jI][cI][idI][aI]->GetBinContent(bIX+1);

	    if(jtPtRaw_h[jI][cI][idI][aI]->GetBinContent(bIX+1)/jtPtRaw_h[jI][cI][0][aI]->GetBinContent(bIX+1) > maxRat) maxRat = jtPtRaw_h[jI][cI][idI][aI]->GetBinContent(bIX+1)/jtPtRaw_h[jI][cI][0][aI]->GetBinContent(bIX+1);
	    if(jtPtRaw_h[jI][cI][idI][aI]->GetBinContent(bIX+1)/jtPtRaw_h[jI][cI][0][aI]->GetBinContent(bIX+1) < minRat) minRat = jtPtRaw_h[jI][cI][idI][aI]->GetBinContent(bIX+1)/jtPtRaw_h[jI][cI][0][aI]->GetBinContent(bIX+1);
	  }
	}

	TLegend* leg_p = new TLegend(0.55, 0.6, 0.9, 0.8);
	leg_p->SetBorderSize(0.0);
	leg_p->SetFillStyle(0);
	leg_p->SetFillColor(0);

	for(Int_t idI = 0; idI < nID; ++idI){
	  jtPtRaw_h[jI][cI][idI][aI]->SetMaximum(max*5.);
	  jtPtRaw_h[jI][cI][idI][aI]->SetMinimum(min/5.);

	  std::string id = idStr.at(idI) + ", N=" + std::to_string((int)(jtPtRaw_h[jI][cI][idI][aI]->GetEntries()));
	  leg_p->AddEntry(jtPtRaw_h[jI][cI][idI][aI], id.c_str(), "P L");

	  if(idI == 0) jtPtRaw_h[jI][cI][idI][aI]->DrawCopy("HIST E1 P");
	  else jtPtRaw_h[jI][cI][idI][aI]->DrawCopy("HIST E1 P SAME");
	}

	TLatex* label_p = new TLatex();
	label_p->SetTextFont(43);
	label_p->SetTextSize(14);
	label_p->SetNDC();

	label_p->DrawLatex(0.5, 0.94, tempStr.c_str());
	label_p->DrawLatex(0.5, 0.88, centStr2.c_str());
	label_p->DrawLatex(0.5, 0.82, jtAbsEtaStr2.c_str());
	
	leg_p->Draw("SAME");

	delete label_p;

	bool doLogX = false;
	if(jtPtRaw_h[jI][cI][0][aI]->GetBinWidth(1)*3 < jtPtRaw_h[jI][cI][0][aI]->GetBinWidth(jtPtRaw_h[jI][cI][0][aI]->GetNbinsX()-1)) doLogX = true;

	gPad->SetLogy();
	if(doLogX) gPad->SetLogx();
	gStyle->SetOptStat(0);

	canv_p->cd();
	pads[1]->cd();

	for(Int_t idI = 1; idI < nID; ++idI){
	  TH1D* clone_p = (TH1D*)jtPtRaw_h[jI][cI][idI][aI]->Clone("temp");
	  setSumW2(clone_p);
	  centerTitles(clone_p);

	  clone_p->Divide(jtPtRaw_h[jI][cI][0][aI]);

	  if(maxRat < 1.05) maxRat = 1.05;
	  if(minRat > 0.95) minRat = 0.95;

	  double interval = maxRat - minRat;
	  maxRat += interval/5.;
	  minRat -= interval/5.;

          clone_p->SetMaximum(maxRat);
          clone_p->SetMinimum(minRat);

	  clone_p->GetYaxis()->SetTitle("Ratio");

          if(idI == 1) clone_p->DrawCopy("HIST E1 P");
          else clone_p->DrawCopy("HIST E1 P SAME");

	  delete clone_p;
        }

	if(doLogX) gPad->SetLogx();
	gStyle->SetOptStat(0);

	canv_p->SaveAs(("pdfDir/jtPtRaw_" + tempStr + "_" + centStr + "_" + jtAbsEtaStr + "_" + dateStr + ".pdf").c_str());

	delete pads[0];
	delete pads[1];
	delete canv_p;
	delete leg_p;
      }
    }
  }  


  for(Int_t jI = 0; jI < nDataJet; ++jI){
    outFile_p->cd();
    dir_p[jI]->cd();

    for(Int_t cI = 0; cI < nCentBins; ++cI){
      for(Int_t idI = 0; idI < nID; ++idI){

	for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
	  jtPtRaw_h[jI][cI][idI][aI]->Write("", TObject::kOverwrite);
	  delete jtPtRaw_h[jI][cI][idI][aI];
	  
	  jtPtRaw_RecoTrunc_h[jI][cI][idI][aI]->Write("", TObject::kOverwrite);
	  delete jtPtRaw_RecoTrunc_h[jI][cI][idI][aI];
	}
      }
    }
  }

  outFile_p->cd();
  TDirectory* cutDir_p = (TDirectory*)outFile_p->mkdir("cutDir");

  cutPropagator cutPropOut;
  cutPropOut.Clean();
  cutPropOut.SetIsPP(isDataPP);
  cutPropOut.SetJtAbsEtaMax(jtAbsEtaMax);
  cutPropOut.SetNJtPtBins(nJtPtBins);
  cutPropOut.SetJtPtBins(nJtPtBins+1, jtPtBins);
  cutPropOut.SetNJtAbsEtaBins(nJtAbsEtaBins);
  cutPropOut.SetJtAbsEtaBinsLow(nJtAbsEtaBins, jtAbsEtaBinsLow);
  cutPropOut.SetJtAbsEtaBinsHi(nJtAbsEtaBins, jtAbsEtaBinsHi);
  cutPropOut.SetNPthats(0);
  cutPropOut.SetPthats({});
  cutPropOut.SetPthatWeights({});
  cutPropOut.SetNCentBins(nCentBins);
  cutPropOut.SetCentBinsLow(centBinsLow);
  cutPropOut.SetCentBinsHi(centBinsHi);
  cutPropOut.SetNID(nID);
  cutPropOut.SetIdStr(idStr);
  cutPropOut.SetJtPfCHMFCutLow(jtPfCHMFCutLow);
  cutPropOut.SetJtPfCHMFCutHi(jtPfCHMFCutHi);
  cutPropOut.SetJtPfMUMFCutLow(jtPfMUMFCutLow);
  cutPropOut.SetJtPfMUMFCutHi(jtPfMUMFCutHi);

  if(!cutPropOut.WriteAllVarToFile(outFile_p, cutDir_p)) std::cout << "Warning: Cut writing has failed" << std::endl;

  outFile_p->Close();
  delete outFile_p;

  responseFile_p->Close();
  delete responseFile_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 3 && argc != 4){
    std::cout << "Usage: ./bin/processRawData.exe <inDataFileName> <inResponseName> <isPP-opt>" << std::endl;
    return 1;
  }

  int retVal = 0;
  if(argc == 3) retVal += processRawData(argv[1], argv[2]);
  else retVal += processRawData(argv[1], argv[2], std::stoi(argv[3]));
  return retVal;
}
