 #include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TDatime.h"

#include "Utility/include/checkMakeDir.h"
#include "Utility/include/getLinBins.h"
#include "Utility/include/histDefUtility.h"
#include "Utility/include/plotUtilities.h"
#include "Utility/include/returnRootFileContentsList.h"
#include "Utility/include/specialHYDJETEventExclude.h"

int checkGenHydSpectra(const std::string inName)
{
  std::vector<std::string> fileList;      
  if(inName.find(".root") != std::string::npos)fileList.push_back(inName);
  else if(inName.find(".txt") != std::string::npos){
    std::ifstream file(inName.c_str());
    std::string tempStr;                  
    while(std::getline(file, tempStr)){
      while(tempStr.find(" ") != std::string::npos){tempStr.replace(tempStr.find(" "), 1, "");}
      if(tempStr.size() == 0) continue;
      if(tempStr.find(".root") != std::string::npos) fileList.push_back(tempStr);
    }
  }

  if(fileList.size() == 0){
    std::cout << "inName \'" << inName << "\' provides no valid root files, return 1." << std::endl;
    return 1;
  }

  TFile* inFile_p = new TFile(fileList.at(0).c_str(), "READ");
  std::vector<std::string> jetList = returnRootFileContentsList(inFile_p, "TTree", "JetAnalyzer");
  inFile_p->Close();
  delete inFile_p;
  inFile_p = NULL;

  const int nRAlgoPerma = 7;
  const int rAlgoPerma[nRAlgoPerma] = {2, 3, 4, 5, 6, 8, 10};
  bool isRAlgoPerma[nRAlgoPerma];
  for(Int_t i = 0; i < nRAlgoPerma; ++i){isRAlgoPerma[i] = false;}

  std::vector<std::string> jetAlgos;
  std::vector<int> jetRs;
  Int_t r4Pos = -1;

  for(unsigned int rI = 0; rI < jetList.size(); ++rI){
    for(Int_t i = 0; i < nRAlgoPerma; ++i){
      if(isRAlgoPerma[i]) continue;

      if(jetList.at(rI).find("akCs" + std::to_string(rAlgoPerma[i])) != std::string::npos){
	jetAlgos.push_back(jetList.at(rI));
	jetRs.push_back(rAlgoPerma[i]);
	isRAlgoPerma[i] = true;
	if(i == 2) r4Pos = jetAlgos.size()-1;
      }
      else if(jetList.at(rI).find("akPu" + std::to_string(rAlgoPerma[i])) != std::string::npos){
	jetAlgos.push_back(jetList.at(rI));
	jetRs.push_back(rAlgoPerma[i]);
	isRAlgoPerma[i] = true;
	if(i == 2) r4Pos = jetAlgos.size()-1;
      }
    }
  }

  if(r4Pos == -1){
    std::cout << "R4Pos: " << r4Pos << ", return 1" << std::endl;
    return 1;
  }

  const Int_t nJetAlgos = jetAlgos.size();

  std::cout << "Processing " << nJetAlgos << " jet algos..." << std::endl;
  for(Int_t jI = 0; jI < nJetAlgos; ++jI){
    std::cout << " " << jI << "/" << nJetAlgos << ": " << jetRs.at(jI) << " with " << jetAlgos.at(jI) << std::endl;
  }
  
  const Int_t nCentBins = 5;
  const Int_t centBinsLow[nCentBins] = {0, 10, 30, 50, 0};
  const Int_t centBinsHi[nCentBins] = {10, 30, 50, 90, 90};  

  const Int_t nJtPtBins = 200;
  Double_t jtPtBins[nJtPtBins+1];
  const Float_t jtPtLow = 100;
  const Float_t jtPtHi = 700;
  getLinBins(jtPtLow, jtPtHi, nJtPtBins, jtPtBins);

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  std::string outFileName = inName;
  while(outFileName.find("/") != std::string::npos){outFileName.replace(0, outFileName.find("/")+1, "");}
  if(outFileName.find(".root") != std::string::npos) outFileName.replace(outFileName.find(".root"), std::string(".root").size(), "");
  else if(outFileName.find(".txt") != std::string::npos) outFileName.replace(outFileName.find(".txt"), std::string(".txt").size(), "");
  outFileName = "output/" + outFileName + "_checkGenHydSpectra_" + dateStr + ".root";
  checkMakeDir("output");

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");

  TH1D* genPt_h[nJetAlgos][nCentBins];
  TH1D* genPt_Exclude_h[nJetAlgos][nCentBins];
  
  for(Int_t jI = 0; jI < nJetAlgos; ++jI){
    std::string genStr = "GenR" + std::to_string(jetRs.at(jI));

    for(Int_t cI = 0; cI < nCentBins; ++cI){
      std::string centStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);

      genPt_h[jI][cI] = new TH1D(("genPt_" + genStr + "_" + centStr + "_h").c_str(), ";Gen. p_{T} (subid != 0);Counts", nJtPtBins, jtPtBins);
      genPt_Exclude_h[jI][cI] = new TH1D(("genPt_" + genStr + "_" + centStr + "_Exclude_h").c_str(), ";Gen. p_{T} (subid != 0);Counts", nJtPtBins, jtPtBins);
      centerTitles({genPt_h[jI][cI], genPt_Exclude_h[jI][cI]});
      setSumW2({genPt_h[jI][cI], genPt_Exclude_h[jI][cI]});
    }
  }

  specialHYDJETEventExclude specialSel;

  Int_t totalEvt = 0;
  
  std::cout << "Processing " << fileList.size() << " files..." << std::endl;
  for(unsigned int fI = 0; fI < fileList.size(); ++fI){
    std::cout << " " << fI << "/" << fileList.size() << ": \'" << fileList.at(fI) << "\'" << std::endl;

    inFile_p = new TFile(fileList.at(fI).c_str(), "READ");
    TTree* hiTree_p = (TTree*)inFile_p->Get("hiEvtAnalyzer/HiTree");
    TTree* jetTrees_p[nJetAlgos];
    for(Int_t jI = 0; jI < nJetAlgos; ++jI){
      jetTrees_p[jI] = (TTree*)inFile_p->Get(jetAlgos.at(jI).c_str());
    }

    Int_t hiBin_;
    
    hiTree_p->SetBranchStatus("*", 0);
    hiTree_p->SetBranchStatus("hiBin", 1);

    hiTree_p->SetBranchAddress("hiBin", &hiBin_);

    const Int_t nMaxJets = 500;
    Int_t ngen_[nJetAlgos];
    Float_t genpt_[nJetAlgos][nMaxJets];
    Float_t genphi_[nJetAlgos][nMaxJets];
    Float_t geneta_[nJetAlgos][nMaxJets];
    Int_t gensubid_[nJetAlgos][nMaxJets];

    for(Int_t jI = 0; jI < nJetAlgos; ++jI){
      jetTrees_p[jI]->SetBranchStatus("*", 0);
      jetTrees_p[jI]->SetBranchStatus("ngen", 1);
      jetTrees_p[jI]->SetBranchStatus("genpt", 1);
      jetTrees_p[jI]->SetBranchStatus("genphi", 1);
      jetTrees_p[jI]->SetBranchStatus("geneta", 1);
      jetTrees_p[jI]->SetBranchStatus("gensubid", 1);

      jetTrees_p[jI]->SetBranchAddress("ngen", &(ngen_[jI]));
      jetTrees_p[jI]->SetBranchAddress("genpt", (genpt_[jI]));
      jetTrees_p[jI]->SetBranchAddress("genphi", (genphi_[jI]));
      jetTrees_p[jI]->SetBranchAddress("geneta", (geneta_[jI]));
      jetTrees_p[jI]->SetBranchAddress("gensubid", (gensubid_[jI]));
    }

    const Int_t nEntries = hiTree_p->GetEntries();
    const Int_t printInterval = TMath::Max(1, nEntries/20);

    for(Int_t entry = 0; entry < nEntries; ++entry){
      if(nEntries >= 50000 && entry%printInterval == 0) std::cout << " Entry: " << entry << "/" << nEntries << std::endl;

      hiTree_p->GetEntry(entry);
      for(Int_t jI = 0; jI < nJetAlgos; ++jI){
	jetTrees_p[jI]->GetEntry(entry);
      }

      std::vector<Int_t> centPos;
      for(Int_t cI = 0; cI < nCentBins; ++cI){
	if(hiBin_/2 >= centBinsLow[cI] && hiBin_/2 < centBinsHi[cI]){
	  centPos.push_back(cI);
	  //	  break;
	}
      }

      if(centPos.size() == 0) continue;
      
      ++totalEvt;
      
      bool isBadEvent = specialSel.CheckEventBadJet(ngen_[r4Pos], genpt_[r4Pos], genphi_[r4Pos], geneta_[r4Pos], gensubid_[r4Pos]);

      Int_t subid = -1;
      Float_t leadingPt = -1;
      Float_t subleadingPt = -1;
    
      for(Int_t jI = 0; jI < nJetAlgos; ++jI){
	for(Int_t gI = 0; gI < ngen_[jI]; ++gI){       
	  if(gensubid_[jI][gI] == 0) continue;

	  //	  if(r4Pos == jI && genpt_[jI][gI] > 200 && !isBadEvent) std::cout << "  " << genpt_[jI][gI] << ", " << entry << std::endl;

	  if(r4Pos == jI){
	    if(gensubid_[jI][gI] != subid){
	      if(!isBadEvent && leadingPt > 150 && subleadingPt < 50) std::cout << "Large asymm jet p_{T}: " << leadingPt << ", " << subleadingPt << ", " << subid << ", " << ", CORRECT, " << entry << std::endl;
	      subid = gensubid_[jI][gI];
	      leadingPt = genpt_[jI][gI];
	      subleadingPt = -1;
	    }
	    else{
	      if(genpt_[jI][gI] > leadingPt){
		subleadingPt = leadingPt;
		leadingPt = genpt_[jI][gI];		
	      }
	      else if(genpt_[jI][gI] > subleadingPt) subleadingPt = genpt_[jI][gI];
	    }
	  }

	  if(genpt_[jI][gI] >= jtPtHi) continue;
	  if(genpt_[jI][gI] < jtPtLow) continue;

	  for(unsigned int cI = 0; cI < centPos.size(); ++cI){
	    genPt_h[jI][centPos.at(cI)]->Fill(genpt_[jI][gI]);
	    if(!isBadEvent) genPt_Exclude_h[jI][centPos.at(cI)]->Fill(genpt_[jI][gI]);
	  }
	}
      }

      if(!isBadEvent && leadingPt > 150 && subleadingPt < 50) std::cout << "Large asymm jet p_{T} (END): " << leadingPt << ", " << subleadingPt << ", " << subid << ", " << entry << std::endl;



      centPos.clear();
    }

    inFile_p->Close();
    delete inFile_p;
    inFile_p = NULL;
  }

  specialSel.PrintExcludedNumbers();
  std::cout << "Total events, fraction: " << totalEvt << ", " << ((double)specialSel.TotalExcludedNumbers())/((double)totalEvt)<< std::endl;


  outFile_p->cd();
  for(Int_t jI = 0; jI < nJetAlgos; ++jI){
    for(Int_t cI = 0; cI < nCentBins; ++cI){
      genPt_h[jI][cI]->Write("", TObject::kOverwrite);
      delete genPt_h[jI][cI];

      genPt_Exclude_h[jI][cI]->Write("", TObject::kOverwrite);
      delete genPt_Exclude_h[jI][cI];
    }
  }
  outFile_p->Close();
  delete outFile_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: checkGenHydSpectra.exe <inName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += checkGenHydSpectra(argv[1]);
  return retVal;
}
