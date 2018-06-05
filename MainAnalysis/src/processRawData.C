#include <iostream>
#include <string>
#include <fstream>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TDatime.h"

#include "Utility/include/goodGlobalSelection.h"
#include "Utility/include/mntToXRootdFileString.h"
#include "Utility/include/plotUtilities.h"
#include "Utility/include/returnRootFileContentsList.h"

int processRawData(const std::string inDataFileName, const std::string inResponseName, bool isDataPP = false)
{
  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  TFile* responseFile_p = new TFile(inResponseName.c_str(), "READ");
  std::vector<std::string> jetDirList = returnRootFileContentsList(responseFile_p, "TDirectoryFile", "JetAnalyzer");
  std::vector<std::string> cutDirList = returnRootFileContentsList(responseFile_p, "TNamed", "");

  std::cout << "Validating " << jetDirList.size() << " jets..." << std::endl;
  for(unsigned int jI = 0; jI < jetDirList.size(); ++jI){
    std::cout << " " << jI << "/" << jetDirList.size() << ": " << jetDirList.at(jI) << std::endl;
  }

  Int_t nCentBins = -1;
  std::vector<Int_t> centBinsLow;
  std::vector<Int_t> centBinsHi;

  Float_t jtAbsEtaMaxTemp = -99.;

  Int_t nJtPtBinsTemp = -1;
  std::vector<Double_t> jtPtBinsTemp;

  Int_t nJtAbsEtaBinsTemp = -1;
  std::vector<Double_t> jtAbsEtaBinsLowTemp;
  std::vector<Double_t> jtAbsEtaBinsHiTemp;

  bool isResponsePP = false;

  std::cout << "Using cuts: " << std::endl;
  for(unsigned int cI = 0; cI < cutDirList.size(); ++cI){
    std::cout << " " << cI << "/" << cutDirList.size() << ": " << cutDirList.at(cI) << std::endl;

    std::string tempStr = cutDirList.at(cI);
    while(tempStr.find("/") != std::string::npos){tempStr.replace(0, tempStr.find("/")+1, "");}

    if(tempStr.find("nCentBins") != std::string::npos && tempStr.size() == std::string("nCentBins").size()) nCentBins = std::stoi(((TNamed*)responseFile_p->Get(cutDirList.at(cI).c_str()))->GetTitle());
    else if(tempStr.find("centBinsLow") != std::string::npos && tempStr.size() == std::string("centBinsLow").size()){
      std::string tempStr2 = ((TNamed*)responseFile_p->Get(cutDirList.at(cI).c_str()))->GetTitle();
      while(tempStr2.find(",") != std::string::npos){
        centBinsLow.push_back(std::stoi(tempStr2.substr(0, tempStr2.find(","))));
        tempStr2.replace(0, tempStr2.find(",")+1, "");
      }
      if(tempStr2.size() != 0) centBinsLow.push_back(std::stoi(tempStr2));
    }
    else if(tempStr.find("centBinsHi") != std::string::npos && tempStr.size() == std::string("centBinsHi").size()){
      std::string tempStr2 = ((TNamed*)responseFile_p->Get(cutDirList.at(cI).c_str()))->GetTitle();
      while(tempStr2.find(",") != std::string::npos){
        centBinsHi.push_back(std::stoi(tempStr2.substr(0, tempStr2.find(","))));
        tempStr2.replace(0, tempStr2.find(",")+1, "");
      }
      if(tempStr2.size() != 0) centBinsHi.push_back(std::stoi(tempStr2));
    }
    else if(tempStr.find("jtAbsEtaMax") != std::string::npos && tempStr.size() == std::string("jtAbsEtaMax").size()) jtAbsEtaMaxTemp = std::stof(((TNamed*)responseFile_p->Get(cutDirList.at(cI).c_str()))->GetTitle());
    else if(tempStr.find("nJtPtBins") != std::string::npos && tempStr.size() == std::string("nJtPtBins").size()) nJtPtBinsTemp = std::stoi(((TNamed*)responseFile_p->Get(cutDirList.at(cI).c_str()))->GetTitle());
    else if(tempStr.find("nJtAbsEtaBins") != std::string::npos && tempStr.size() == std::string("nJtAbsEtaBins").size()) nJtAbsEtaBinsTemp = std::stoi(((TNamed*)responseFile_p->Get(cutDirList.at(cI).c_str()))->GetTitle());
    else if(tempStr.find("isPP") != std::string::npos && tempStr.size() == std::string("isPP").size()) isResponsePP = std::stoi(((TNamed*)responseFile_p->Get(cutDirList.at(cI).c_str()))->GetTitle());
    else if(tempStr.find("jtPtBins") != std::string::npos && tempStr.size() == std::string("jtPtBins").size()){
      std::string tempStr2 = ((TNamed*)responseFile_p->Get(cutDirList.at(cI).c_str()))->GetTitle();
      while(tempStr2.find(",") != std::string::npos){
        jtPtBinsTemp.push_back(std::stod(tempStr2.substr(0, tempStr2.find(","))));
        tempStr2.replace(0, tempStr2.find(",")+1, "");
      }
      if(tempStr2.size() != 0) jtPtBinsTemp.push_back(std::stod(tempStr2));
    }
    else if(tempStr.find("jtAbsEtaBinsLow") != std::string::npos && tempStr.size() == std::string("jtAbsEtaBinsLow").size()){
      std::string tempStr2 = ((TNamed*)responseFile_p->Get(cutDirList.at(cI).c_str()))->GetTitle();
      while(tempStr2.find(",") != std::string::npos){
        jtAbsEtaBinsLowTemp.push_back(std::stod(tempStr2.substr(0, tempStr2.find(","))));
        tempStr2.replace(0, tempStr2.find(",")+1, "");
      }
      if(tempStr2.size() != 0) jtAbsEtaBinsLowTemp.push_back(std::stod(tempStr2));
    }
    else if(tempStr.find("jtAbsEtaBinsHi") != std::string::npos && tempStr.size() == std::string("jtAbsEtaBinsHi").size()){
      std::string tempStr2 = ((TNamed*)responseFile_p->Get(cutDirList.at(cI).c_str()))->GetTitle();
      while(tempStr2.find(",") != std::string::npos){
        jtAbsEtaBinsHiTemp.push_back(std::stod(tempStr2.substr(0, tempStr2.find(","))));
        tempStr2.replace(0, tempStr2.find(",")+1, "");
      }
      if(tempStr2.size() != 0) jtAbsEtaBinsHiTemp.push_back(std::stod(tempStr2));
    }
  }

  if(nCentBins < 0) std::cout << "nCentBins less than 0. please check input file. return 1" << std::endl;
  if(nJtPtBinsTemp < 0) std::cout << "nJtPtBinsTemp less than 0. please check input file. return 1" << std::endl;
  if(nJtAbsEtaBinsTemp < 0) std::cout << "nJtAbsEtaBinsTemp less than 0. please check input file. return 1" << std::endl;
  if(jtAbsEtaMaxTemp < -98) std::cout << "jtAbsEtaMaxTemp less than -98. please check input file. return 1" << std::endl;

  if(nCentBins < 0 || nJtPtBinsTemp < 0 || nJtAbsEtaBinsTemp < 0 || jtAbsEtaMaxTemp < -98){
    responseFile_p->Close();
    delete responseFile_p;
    return 1;
  }

  const Float_t jtAbsEtaMax = jtAbsEtaMaxTemp;

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

  outFileName = "output/" + outFileName + "_" + outFileName2 + "_ProcessRawData_" + dateStr + ".root";

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  TDirectory* dir_p[nDataJet];
  TH1D* jtPtRaw_h[nDataJet][nCentBins][nJtAbsEtaBins];
  TH1D* jtPtRaw_RecoTrunc_h[nDataJet][nCentBins][nJtAbsEtaBins];

  for(Int_t jI = 0; jI < nDataJet; ++jI){
    std::string tempStr = dataTreeList.at(jI);
    tempStr.replace(tempStr.find("/"), tempStr.size() - tempStr.find("/"), "");
    dir_p[jI] = (TDirectory*)outFile_p->mkdir(tempStr.c_str());

    for(Int_t cI = 0; cI < nCentBins; ++cI){
      const std::string centStr = "Cent" + std::to_string(centBinsLow.at(cI)) + "to" + std::to_string(centBinsHi.at(cI));

      for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
        const std::string jtAbsEtaStr = "AbsEta" + prettyString(jtAbsEtaBinsLow[aI], 1, true) + "to" + prettyString(jtAbsEtaBinsHi[aI], 1, true);
	jtPtRaw_h[jI][cI][aI] = new TH1D(("jtPtRaw_" + tempStr + "_" + centStr + "_" + jtAbsEtaStr + "_h").c_str(), ";Raw Jet p_{T};Counts", nJtPtBins, jtPtBins);
	jtPtRaw_RecoTrunc_h[jI][cI][aI] = new TH1D(("jtPtRaw_RecoTrunc_" + tempStr + "_" + centStr + "_" + jtAbsEtaStr + "_h").c_str(), ";Raw Jet p_{T};Counts", nJtPtBins, jtPtBins);
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

      jetTrees_p[jI]->SetBranchAddress("nref", &(nref_[jI]));
      jetTrees_p[jI]->SetBranchAddress("jtpt", jtpt_[jI]);
      jetTrees_p[jI]->SetBranchAddress("jtphi", jtphi_[jI]);
      jetTrees_p[jI]->SetBranchAddress("jteta", jteta_[jI]);
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

	  for(unsigned int aI = 0; aI < jtAbsEtaPoses.size(); ++aI){
	    bool goodReco = (jtpt_[tI][jI] >= jtPtBins[0] && jtpt_[tI][jI] < jtPtBins[nJtPtBins]);
	    bool goodRecoTrunc = (jtpt_[tI][jI] >= jtPtBins[1] && jtpt_[tI][jI] < jtPtBins[nJtPtBins-1]);

	    if(goodReco){
	      jtPtRaw_h[tI][centPos][aI]->Fill(jtpt_[tI][jI]);
	      if(goodRecoTrunc) jtPtRaw_RecoTrunc_h[tI][centPos][aI]->Fill(jtpt_[tI][jI]);
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
    outFile_p->cd();
    dir_p[jI]->cd();

    for(Int_t cI = 0; cI < nCentBins; ++cI){

      for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
	jtPtRaw_h[jI][cI][aI]->Write("", TObject::kOverwrite);
	delete jtPtRaw_h[jI][cI][aI];

	jtPtRaw_RecoTrunc_h[jI][cI][aI]->Write("", TObject::kOverwrite);
	delete jtPtRaw_RecoTrunc_h[jI][cI][aI];
      }
    }
  }

  outFile_p->Close();
  delete outFile_p;

  responseFile_p->Close();
  delete responseFile_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 3){
    std::cout << "Usage: ./bin/processRawData.exe <inDataFileName> <inResponseName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += processRawData(argv[1], argv[2]);
  return retVal;
}
