#include <iostream>
#include <string>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TDatime.h"
#include "TMath.h"
#include "TRandom3.h"

#include "src/RooUnfoldResponse.h"
#include "src/RooUnfoldBayes.h"

#include "Utility/include/checkMakeDir.h"
#include "Utility/include/getLinBins.h"
#include "Utility/include/plotUtilities.h"
#include "Utility/include/histDefUtility.h"
#include "Utility/include/returnRootFileContentsList.h"

#include "MCExcample/include/atlasRAA.h"

void doUnfold(TH1D ** inUnfoldHist_p, TH1* inRawHist_p, RooUnfoldResponse* response_p, const Int_t bI, const std::string cloneStr)
{
  RooUnfoldBayes* temp = new RooUnfoldBayes(response_p, inRawHist_p, bI+1);
  temp->SetVerbose(0);
  (*inUnfoldHist_p) = (TH1D*)temp->Hreco(RooUnfold::kCovToy)->Clone(cloneStr.c_str());
  centerTitles(*inUnfoldHist_p);
  setSumW2(*inUnfoldHist_p);
  delete temp;

  return;
}

int makeRawRAASpectra(const std::string inFileName)
{
  if(!checkFile(inFileName)){
    std::cout << "Given inFileName \'" << inFileName << "\' is invalid. return 1" << std::endl;
    return 1;
  }
  else if(inFileName.find(".root") == std::string::npos){
    std::cout << "Given inFileName \'" << inFileName << "\' does not contain \'.root\'. return 1" << std::endl;
    return 1;
  }


  TRandom3* randGen_p = new TRandom3(0);

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  checkMakeDir("output");
  std::string outFileName = inFileName;
  while(outFileName.find("/") != std::string::npos){outFileName.replace(0, outFileName.find("/")+1, "");}
  while(outFileName.find(".root") != std::string::npos){outFileName.replace(outFileName.find(".root"), 5, "");}
  outFileName = "output/" + outFileName + "_RawRAA_" + dateStr + ".root";

  const Int_t nRVals = 2;
  const Double_t rVals[nRVals] = {0.4, 0.8};
  const Int_t nRCFakes[nRVals] = {(Int_t)(2.*2.*2./(rVals[0]*rVals[0])), (Int_t)(2.*2.*2./(rVals[1]*rVals[1]))};
  const Float_t cVal = .06;
  const Float_t sVal = 1.0;
  const Float_t rcWidth[nRVals] = {20., 40.};
  const Int_t nBayes = 6;

  const Int_t nLogFactor = 100;
  const Float_t logFactorLow = 2.;
  const Float_t logFactorHi = 7.;
  Double_t logFactor[nLogFactor+1];
  getLinBins(logFactorLow, logFactorHi, nLogFactor, logFactor);
  Double_t logFactorBins[nLogFactor+2];

  Double_t delta = (logFactor[1] - logFactor[0])/2.;
  for(Int_t lI = 0; lI < nLogFactor+1; ++lI){
    logFactorBins[lI] = logFactor[lI] - delta;
  }
  logFactorBins[nLogFactor+1] = logFactor[nLogFactor] + delta;


  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");

  atlasRAA raa;
  TH1D* atlasRAA_R0p4_h=NULL;
  raa.SetHistogram(&atlasRAA_R0p4_h, "atlasRAA_R0p4_h");
  
  const Int_t nAtlasBins = atlasRAA_R0p4_h->GetNbinsX();
  Double_t atlasBins[nAtlasBins+1];
  for(Int_t bI = 0; bI < nAtlasBins+1; ++bI){
    atlasBins[bI] = atlasRAA_R0p4_h->GetBinLowEdge(bI+1);
  }

  TH1D* genJtPt_AtlasBinnedR0p4_h = new TH1D("genJtPt_AtlasBinnedR0p4_h", ";;", nAtlasBins, atlasBins);
  TH1D* genJtPt_AtlasBinnedR0p4_LogLoss_h[nLogFactor+1];
  TH1D* raa_AtlasBinnedR0p4_LogLoss_h[nLogFactor+1];
  TH1D* atlasChi2_h = new TH1D("atlasChi2_h", ";;", nLogFactor+1, logFactorBins);

  centerTitles(genJtPt_AtlasBinnedR0p4_h);
  setSumW2(genJtPt_AtlasBinnedR0p4_h);
  
  for(Int_t lI = 0; lI < nLogFactor+1; ++lI){
    genJtPt_AtlasBinnedR0p4_LogLoss_h[lI] = new TH1D(("genJtPt_AtlasBinnedR0p4_LogLoss" + prettyString(logFactor[lI], 2, true) + "_h").c_str(), ";Gen Jet p_{T};Counts", nAtlasBins, atlasBins);

    raa_AtlasBinnedR0p4_LogLoss_h[lI] = new TH1D(("raa_AtlasBinnedR0p4_LogLoss" + prettyString(logFactor[lI], 2, true) + "_h").c_str(), ";Gen Jet p_{T};R_{AA}", nAtlasBins, atlasBins);

    centerTitles({genJtPt_AtlasBinnedR0p4_LogLoss_h[lI], raa_AtlasBinnedR0p4_LogLoss_h[lI]});
    setSumW2({genJtPt_AtlasBinnedR0p4_LogLoss_h[lI], raa_AtlasBinnedR0p4_LogLoss_h[lI]});
  }
  

  TH1D* genJtPt_h[nRVals];
  TH1D* genJtPt_LogLoss_h[nRVals][nLogFactor+1];
  TH1D* genJtPt_GoodRecoFull_h[nRVals];
  TH1D* recoJtPt_h[nRVals];
  TH1D* recoJtPt_RCFake_h[nRVals];
  TH1D* recoJtPt_RCFakeCheck_h[nRVals];
  TH1D* recoJtPt_NoGen_h[nRVals];
  TH1D* genJtPt_Parallel_h[nRVals];
  TH1D* genJtPt_GoodRecoFull_Parallel_h[nRVals];
  TH1D* recoJtPt_Parallel_h[nRVals];
  TH1D* recoJtPt_RCFake_Parallel_h[nRVals];
  TH1D* recoJtPt_RCFakeCheck_Parallel_h[nRVals];
  TH1D* recoJtPt_NoGen_Parallel_h[nRVals];
  TH2D* recoVsGenJtPt_h[nRVals];

  TH2D* recoVsGenJtPtRes0p9_h[nRVals];
  TH2D* recoVsGenJtPtRes0p85_h[nRVals];
  TH2D* recoVsGenJtPtRes1p1_h[nRVals];
  TH2D* recoVsGenJtPtRes1p15_h[nRVals];

  RooUnfoldResponse* response_h[nRVals];
  RooUnfoldResponse* responseRes0p9_h[nRVals];
  RooUnfoldResponse* responseRes0p85_h[nRVals];
  RooUnfoldResponse* responseRes1p1_h[nRVals];
  RooUnfoldResponse* responseRes1p15_h[nRVals];
  RooUnfoldResponse* responseFake_h[nRVals];

  TH1D* unfoldedJtPt_h[nRVals][nBayes];
  TH1D* unfoldedJtPt_Parallel_h[nRVals][nBayes];
  TH1D* unfoldedJtPtRes0p9_h[nRVals][nBayes];
  TH1D* unfoldedJtPtRes0p9_Parallel_h[nRVals][nBayes];
  TH1D* unfoldedJtPtRes0p85_h[nRVals][nBayes];
  TH1D* unfoldedJtPtRes0p85_Parallel_h[nRVals][nBayes];
  TH1D* unfoldedJtPtRes1p1_h[nRVals][nBayes];
  TH1D* unfoldedJtPtRes1p1_Parallel_h[nRVals][nBayes];
  TH1D* unfoldedJtPtRes1p15_h[nRVals][nBayes];
  TH1D* unfoldedJtPtRes1p15_Parallel_h[nRVals][nBayes];
  TH1D* unfoldedJtPtFake_h[nRVals][nBayes];
  TH1D* unfoldedJtPtFake_Parallel_h[nRVals][nBayes];

  const Int_t nJtPtBins = 10;
  Float_t jtPtBins[nJtPtBins+1] = {100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100};

  /*
  const Int_t nJtPtBins = 12;
  Float_t jtPtBins[nJtPtBins+1] = {100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700};
  */

  for(Int_t rI = 0; rI < nRVals; ++rI){
    genJtPt_h[rI] = new TH1D(("genJtPt_R" + prettyString(rVals[rI], 1, true) + "_h").c_str(), ";Gen Jet p_{T};Counts", nJtPtBins, jtPtBins);
    for(Int_t lI = 0; lI < nLogFactor+1; ++lI){
      genJtPt_LogLoss_h[rI][lI] = new TH1D(("genJtPt_LogLoss" + prettyString(logFactor[lI], 2, true) + "_R" + prettyString(rVals[rI], 1, true) + "_h").c_str(), ";Gen Jet p_{T};Counts", nJtPtBins, jtPtBins);
    }

    genJtPt_GoodRecoFull_h[rI] = new TH1D(("genJtPt_GoodRecoFull_R" + prettyString(rVals[rI], 1, true) + "_h").c_str(), ";Gen Jet p_{T};Counts", nJtPtBins, jtPtBins);
    recoJtPt_h[rI] = new TH1D(("recoJtPt_R" + prettyString(rVals[rI], 1, true) + "_h").c_str(), ";Reco Jet p_{T};Counts", nJtPtBins, jtPtBins);
    recoJtPt_RCFake_h[rI] = new TH1D(("recoJtPt_R" + prettyString(rVals[rI], 1, true) + "_RCFake_h").c_str(), ";Reco Jet p_{T};Counts", nJtPtBins, jtPtBins);
    recoJtPt_RCFakeCheck_h[rI] = new TH1D(("recoJtPt_R" + prettyString(rVals[rI], 1, true) + "_RCFakeCheck_h").c_str(), ";Reco Jet p_{T};Counts", 200, -jtPtBins[nJtPtBins], jtPtBins[nJtPtBins]);
    recoJtPt_NoGen_h[rI] = new TH1D(("recoJtPt_R" + prettyString(rVals[rI], 1, true) + "_NoGen_h").c_str(), ";Reco Jet p_{T};Counts", nJtPtBins, jtPtBins);

    genJtPt_Parallel_h[rI] = new TH1D(("genJtPt_R" + prettyString(rVals[rI], 1, true) + "_Parallel_h").c_str(), ";Gen Jet p_{T};Counts", nJtPtBins, jtPtBins);
    genJtPt_GoodRecoFull_Parallel_h[rI] = new TH1D(("genJtPt_GoodRecoFull_R" + prettyString(rVals[rI], 1, true) + "_Parallel_h").c_str(), ";Gen Jet p_{T};Counts", nJtPtBins, jtPtBins);
    recoJtPt_Parallel_h[rI] = new TH1D(("recoJtPt_R" + prettyString(rVals[rI], 1, true) + "_Parallel_h").c_str(), ";Reco Jet p_{T};Counts", nJtPtBins, jtPtBins);
    recoJtPt_RCFake_Parallel_h[rI] = new TH1D(("recoJtPt_R" + prettyString(rVals[rI], 1, true) + "_RCFake_Parallel_h").c_str(), ";Reco Jet p_{T};Counts", nJtPtBins, jtPtBins);
    recoJtPt_RCFakeCheck_Parallel_h[rI] = new TH1D(("recoJtPt_R" + prettyString(rVals[rI], 1, true) + "_RCFakeCheck_Parallel_h").c_str(), ";Reco Jet p_{T};Counts", 200, -jtPtBins[nJtPtBins], jtPtBins[nJtPtBins]);
    recoJtPt_NoGen_Parallel_h[rI] = new TH1D(("recoJtPt_NoGen_R" + prettyString(rVals[rI], 1, true) + "_Parallel_h").c_str(), ";Reco Jet p_{T};Counts", nJtPtBins, jtPtBins);

    recoVsGenJtPt_h[rI] = new TH2D(("recoVsGenJtPt_R" + prettyString(rVals[rI], 1, true) + "_h").c_str(), ";Reco Jet p_{T};Gen. Jet p_{T}", nJtPtBins, jtPtBins, nJtPtBins, jtPtBins);
    recoVsGenJtPtRes0p9_h[rI] = new TH2D(("recoVsGenJtPt_Res0p9_R" + prettyString(rVals[rI], 1, true) + "_h").c_str(), ";Reco Jet p_{T};Gen. Jet p_{T}", nJtPtBins, jtPtBins, nJtPtBins, jtPtBins);
    recoVsGenJtPtRes0p85_h[rI] = new TH2D(("recoVsGenJtPt_Res0p85_R" + prettyString(rVals[rI], 1, true) + "_h").c_str(), ";Reco Jet p_{T};Gen. Jet p_{T}", nJtPtBins, jtPtBins, nJtPtBins, jtPtBins);
    recoVsGenJtPtRes1p1_h[rI] = new TH2D(("recoVsGenJtPt_Res1p1_R" + prettyString(rVals[rI], 1, true) + "_h").c_str(), ";Reco Jet p_{T};Gen. Jet p_{T}", nJtPtBins, jtPtBins, nJtPtBins, jtPtBins);
    recoVsGenJtPtRes1p15_h[rI] = new TH2D(("recoVsGenJtPt_Res1p15_R" + prettyString(rVals[rI], 1, true) + "_h").c_str(), ";Reco Jet p_{T};Gen. Jet p_{T}", nJtPtBins, jtPtBins, nJtPtBins, jtPtBins);

    response_h[rI] = new RooUnfoldResponse(("response_" + prettyString(rVals[rI], 1, true)).c_str(), "");
    response_h[rI]->Setup(recoJtPt_h[rI], genJtPt_h[rI]);

    responseRes0p9_h[rI] = new RooUnfoldResponse(("responseRes0p9_" + prettyString(rVals[rI], 1, true)).c_str(), "");
    responseRes0p9_h[rI]->Setup(recoJtPt_h[rI], genJtPt_h[rI]);

    responseRes0p85_h[rI] = new RooUnfoldResponse(("responseRes0p85_" + prettyString(rVals[rI], 1, true)).c_str(), "");
    responseRes0p85_h[rI]->Setup(recoJtPt_h[rI], genJtPt_h[rI]);

    responseRes1p1_h[rI] = new RooUnfoldResponse(("responseRes1p1_" + prettyString(rVals[rI], 1, true)).c_str(), "");
    responseRes1p1_h[rI]->Setup(recoJtPt_h[rI], genJtPt_h[rI]);

    responseRes1p15_h[rI] = new RooUnfoldResponse(("responseRes1p15_" + prettyString(rVals[rI], 1, true)).c_str(), "");
    responseRes1p15_h[rI]->Setup(recoJtPt_h[rI], genJtPt_h[rI]);

    responseFake_h[rI] = new RooUnfoldResponse(("responseFake_" + prettyString(rVals[rI], 1, true)).c_str(), "");
    responseFake_h[rI]->Setup(recoJtPt_h[rI], genJtPt_h[rI]);
    
    for(Int_t bI = 0; bI < nBayes; ++bI){
      unfoldedJtPt_h[rI][bI] = NULL;
      unfoldedJtPtRes0p9_h[rI][bI] = NULL;
      unfoldedJtPtRes0p85_h[rI][bI] = NULL;
      unfoldedJtPtRes1p1_h[rI][bI] = NULL;
      unfoldedJtPtRes1p15_h[rI][bI] = NULL;
      unfoldedJtPtFake_h[rI][bI] = NULL;

      unfoldedJtPt_Parallel_h[rI][bI] = NULL;
      unfoldedJtPtRes0p9_Parallel_h[rI][bI] = NULL;
      unfoldedJtPtRes0p85_Parallel_h[rI][bI] = NULL;
      unfoldedJtPtRes1p1_Parallel_h[rI][bI] = NULL;
      unfoldedJtPtRes1p15_Parallel_h[rI][bI] = NULL;
      unfoldedJtPtFake_Parallel_h[rI][bI] = NULL;
    }

    setSumW2({genJtPt_h[rI], genJtPt_GoodRecoFull_h[rI], recoJtPt_h[rI], recoJtPt_RCFake_h[rI], recoJtPt_RCFakeCheck_h[rI], recoJtPt_NoGen_h[rI], genJtPt_Parallel_h[rI], genJtPt_GoodRecoFull_Parallel_h[rI], recoJtPt_Parallel_h[rI], recoJtPt_RCFake_Parallel_h[rI], recoJtPt_RCFakeCheck_Parallel_h[rI], recoJtPt_NoGen_Parallel_h[rI]});
    centerTitles({genJtPt_h[rI], genJtPt_GoodRecoFull_h[rI], recoJtPt_h[rI], recoJtPt_RCFake_h[rI], recoJtPt_RCFakeCheck_h[rI], recoJtPt_NoGen_h[rI], genJtPt_Parallel_h[rI], genJtPt_GoodRecoFull_Parallel_h[rI], recoJtPt_Parallel_h[rI], recoJtPt_RCFake_Parallel_h[rI], recoJtPt_RCFakeCheck_Parallel_h[rI], recoJtPt_NoGen_Parallel_h[rI]});

    for(Int_t lI = 0; lI < nLogFactor+1; ++lI){
      setSumW2(genJtPt_LogLoss_h[rI][lI]);
      centerTitles(genJtPt_LogLoss_h[rI][lI]);
    }

    setSumW2({recoVsGenJtPt_h[rI], recoVsGenJtPtRes0p9_h[rI], recoVsGenJtPtRes1p1_h[rI], recoVsGenJtPtRes0p85_h[rI], recoVsGenJtPtRes1p15_h[rI]});
    centerTitles({recoVsGenJtPt_h[rI], recoVsGenJtPtRes0p9_h[rI], recoVsGenJtPtRes1p1_h[rI], recoVsGenJtPtRes0p85_h[rI], recoVsGenJtPtRes1p15_h[rI]});
  }

  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  std::vector<std::string> listOfTNamed = returnRootFileContentsList(inFile_p, "TNamed", "");

  unsigned int pos = 0;
  while(pos < listOfTNamed.size()){
    bool isUnique = true;
    for(unsigned int i = pos+1; i < listOfTNamed.size(); ++i){
      if(listOfTNamed.at(i).size() == listOfTNamed.at(pos).size() && listOfTNamed.at(i).find(listOfTNamed.at(pos)) != std::string::npos){
        isUnique = false;
        listOfTNamed.erase(listOfTNamed.begin()+i);
      }
    }

    if(isUnique) ++pos;
  }

  Int_t doRecoJetsPos = -1;

  std::cout << "Check tnamed: " << std::endl;
  for(unsigned int i = 0; i < listOfTNamed.size(); ++i){
    std::cout << " " << i << "/" << listOfTNamed.size() << ": " << listOfTNamed.at(i) << std::endl;

    std::string tempStr = listOfTNamed.at(i);
    while(tempStr.find("/") != std::string::npos){tempStr.replace(0, tempStr.find("/")+1, "");}

    if(tempStr.find("doRecoJets") != std::string::npos && tempStr.size() == std::string("doRecoJets").size()) doRecoJetsPos = i;
  }

  const Bool_t doRecoJets = std::stoi(std::string(((TNamed*)inFile_p->Get(listOfTNamed.at(doRecoJetsPos).c_str()))->GetTitle()));

  TTree* partonTree_p = (TTree*)inFile_p->Get("partonTree");
  TTree* inTree_p[nRVals];

  Float_t pthat_;
  Float_t pthatWeight_;

  const Int_t nJtMax = 500;
  Int_t nJt_[nRVals];
  Float_t genJtPt_[nRVals][nJtMax];
  Float_t genJtPhi_[nRVals][nJtMax];
  Float_t genJtEta_[nRVals][nJtMax];
  Float_t recoJtPt_[nRVals][nJtMax];
  Float_t recoJtPtRes0p9_[nRVals][nJtMax];
  Float_t recoJtPtRes1p1_[nRVals][nJtMax];
  Float_t recoJtPtRes0p85_[nRVals][nJtMax];
  Float_t recoJtPtRes1p15_[nRVals][nJtMax];
  
  partonTree_p->SetBranchStatus("*", 0);
  partonTree_p->SetBranchStatus("pthat", 1);
  partonTree_p->SetBranchStatus("pthatWeight", 1);

  partonTree_p->SetBranchAddress("pthat", &pthat_);
  partonTree_p->SetBranchAddress("pthatWeight", &pthatWeight_);

  for(Int_t rI = 0; rI < nRVals; ++rI){
    inTree_p[rI] = (TTree*)inFile_p->Get(("recoAndGenTreeR" + prettyString(rVals[rI], 1, true)).c_str());
    
    inTree_p[rI]->SetBranchAddress("nJt", &(nJt_[rI]));
    inTree_p[rI]->SetBranchAddress("genJtPt", genJtPt_[rI]);
    inTree_p[rI]->SetBranchAddress("genJtPhi", genJtPhi_[rI]);
    inTree_p[rI]->SetBranchAddress("genJtEta", genJtEta_[rI]);

    if(doRecoJets){
      inTree_p[rI]->SetBranchAddress("recoJtPt", recoJtPt_[rI]);
      inTree_p[rI]->SetBranchAddress("recoJtPt0p9Res", recoJtPtRes0p9_[rI]);
      inTree_p[rI]->SetBranchAddress("recoJtPt0p85Res", recoJtPtRes0p85_[rI]);
      inTree_p[rI]->SetBranchAddress("recoJtPt1p1Res", recoJtPtRes1p1_[rI]);
      inTree_p[rI]->SetBranchAddress("recoJtPt1p15Res", recoJtPtRes1p15_[rI]);
    }
  }

  const Int_t nEntries = TMath::Min((Int_t)1000000000, (Int_t)inTree_p[0]->GetEntries());
  Int_t fills = 0;
  Int_t paraFills = 0;

  std::cout << "Processing Entries: " << nEntries << "..." << std::endl;
  for(Int_t entry = 0; entry < nEntries; ++entry){
    if(entry%10000 == 0) std::cout << " Entry: " << entry << "/" << nEntries << std::endl;

    for(Int_t rI = 0; rI < nRVals; ++rI){
      partonTree_p->GetEntry(entry);
      inTree_p[rI]->GetEntry(entry);

      bool fillPara = randGen_p->Uniform(0.0, 1.0) >= 0.5;

      if(!doRecoJets){
	for(Int_t jI = 0; jI < nJt_[rI]; ++jI){
	  Double_t sigma = TMath::Sqrt(cVal*cVal + sVal*sVal/genJtPt_[rI][jI] + rcWidth[rI]*rcWidth[rI]/(genJtPt_[rI][jI]*genJtPt_[rI][jI]));

	  recoJtPt_[rI][jI] = genJtPt_[rI][jI]*randGen_p->Gaus(1.0, sigma);
	  recoJtPtRes0p9_[rI][jI] = genJtPt_[rI][jI]*randGen_p->Gaus(1.0, sigma*0.9);
	  recoJtPtRes0p85_[rI][jI] = genJtPt_[rI][jI]*randGen_p->Gaus(1.0, sigma*0.85);
	  recoJtPtRes1p1_[rI][jI] = genJtPt_[rI][jI]*randGen_p->Gaus(1.0, sigma*1.1);
	  recoJtPtRes1p15_[rI][jI] = genJtPt_[rI][jI]*randGen_p->Gaus(1.0, sigma*1.15);
	}
      }

      if(fillPara) paraFills++;
      else fills++;

      for(Int_t jI = 0; jI < nJt_[rI]; ++jI){
	bool goodRecoFull = recoJtPt_[rI][jI] >= jtPtBins[0] && recoJtPt_[rI][jI] < jtPtBins[nJtPtBins];
	bool goodRecoTrunc = recoJtPt_[rI][jI] >= jtPtBins[1] && recoJtPt_[rI][jI] < jtPtBins[nJtPtBins-1];

	bool goodRecoTrunc0p9 = recoJtPtRes0p9_[rI][jI] >= jtPtBins[1] && recoJtPtRes0p9_[rI][jI] < jtPtBins[nJtPtBins-1];
	bool goodRecoTrunc1p1 = recoJtPtRes1p1_[rI][jI] >= jtPtBins[1] && recoJtPtRes1p1_[rI][jI] < jtPtBins[nJtPtBins-1];
	bool goodRecoTrunc0p85 = recoJtPtRes0p85_[rI][jI] >= jtPtBins[1] && recoJtPtRes0p85_[rI][jI] < jtPtBins[nJtPtBins-1];
	bool goodRecoTrunc1p15 = recoJtPtRes1p15_[rI][jI] >= jtPtBins[1] && recoJtPtRes1p15_[rI][jI] < jtPtBins[nJtPtBins-1];

	bool goodGen = genJtPt_[rI][jI] >= jtPtBins[0] && genJtPt_[rI][jI] < jtPtBins[nJtPtBins];
	bool goodGenAtlas = genJtPt_[rI][jI] >= atlasBins[0] && genJtPt_[rI][jI] < atlasBins[nAtlasBins];

	if(!fillPara){
	  if(goodGen){
	    genJtPt_h[rI]->Fill(genJtPt_[rI][jI], pthatWeight_);
	    if(goodRecoFull) genJtPt_GoodRecoFull_h[rI]->Fill(genJtPt_[rI][jI], pthatWeight_);
	  }

	  if(goodGenAtlas) genJtPt_AtlasBinnedR0p4_h->Fill(genJtPt_[rI][jI], pthatWeight_);
	  
	  for(Int_t lI = 0; lI < nLogFactor+1; ++lI){
	    Double_t genJtPt_LogLoss_ = genJtPt_[rI][jI] - logFactor[lI]*TMath::Log(genJtPt_[rI][jI]);
	    bool goodGenLogLoss = genJtPt_LogLoss_ >= jtPtBins[0] && genJtPt_LogLoss_ < jtPtBins[nJtPtBins];
	    bool goodGenLogLossAtlas = genJtPt_LogLoss_ >= atlasBins[0] && genJtPt_LogLoss_ < atlasBins[nAtlasBins];
	    if(goodGenLogLoss) genJtPt_LogLoss_h[rI][lI]->Fill(genJtPt_LogLoss_, pthatWeight_);	    

	    if(goodGenLogLossAtlas) genJtPt_AtlasBinnedR0p4_LogLoss_h[lI]->Fill(genJtPt_LogLoss_, pthatWeight_);
	  }	  

	  if(goodRecoTrunc){
	    recoJtPt_h[rI]->Fill(recoJtPt_[rI][jI], pthatWeight_);
	    if(!goodGen) recoJtPt_NoGen_h[rI]->Fill(recoJtPt_[rI][jI], pthatWeight_);
	  }
	  
	  if(goodRecoTrunc && goodGen){
	    recoVsGenJtPt_h[rI]->Fill(recoJtPt_[rI][jI], genJtPt_[rI][jI], pthatWeight_);
	    response_h[rI]->Fill(recoJtPt_[rI][jI], genJtPt_[rI][jI], pthatWeight_);

	    responseFake_h[rI]->Fill(recoJtPt_[rI][jI], genJtPt_[rI][jI], pthatWeight_);
	  }
	  else if(goodGen){
	    response_h[rI]->Miss(genJtPt_[rI][jI], pthatWeight_);
	    responseFake_h[rI]->Miss(genJtPt_[rI][jI], pthatWeight_);
	  }
	  else if(goodRecoTrunc){
	    response_h[rI]->Fake(recoJtPt_[rI][jI], pthatWeight_);
	    responseFake_h[rI]->Fake(recoJtPt_[rI][jI], pthatWeight_);
	  }

	  if(goodRecoTrunc0p9 && goodGen){
	    recoVsGenJtPtRes0p9_h[rI]->Fill(recoJtPtRes0p9_[rI][jI], genJtPt_[rI][jI], pthatWeight_);
	    responseRes0p9_h[rI]->Fill(recoJtPtRes0p9_[rI][jI], genJtPt_[rI][jI], pthatWeight_);
	  }
	  else if(goodGen) responseRes0p9_h[rI]->Miss(genJtPt_[rI][jI], pthatWeight_);
	  else if(goodRecoTrunc0p9) responseRes0p9_h[rI]->Fake(recoJtPtRes0p9_[rI][jI], pthatWeight_);

	  if(goodRecoTrunc0p85 && goodGen){
	    recoVsGenJtPtRes0p85_h[rI]->Fill(recoJtPtRes0p85_[rI][jI], genJtPt_[rI][jI], pthatWeight_);
	    responseRes0p85_h[rI]->Fill(recoJtPtRes0p85_[rI][jI], genJtPt_[rI][jI], pthatWeight_);
	  }
	  else if(goodGen) responseRes0p85_h[rI]->Miss(genJtPt_[rI][jI], pthatWeight_);
	  else if(goodRecoTrunc0p85) responseRes0p85_h[rI]->Fake(recoJtPtRes0p85_[rI][jI], pthatWeight_);

	  if(goodRecoTrunc1p1 && goodGen){
	    recoVsGenJtPtRes1p1_h[rI]->Fill(recoJtPtRes1p1_[rI][jI], genJtPt_[rI][jI], pthatWeight_);
	    responseRes1p1_h[rI]->Fill(recoJtPtRes1p1_[rI][jI], genJtPt_[rI][jI], pthatWeight_);
	  }
	  else if(goodGen) responseRes1p1_h[rI]->Miss(genJtPt_[rI][jI], pthatWeight_);
	  else if(goodRecoTrunc1p1) responseRes1p1_h[rI]->Fake(recoJtPtRes1p1_[rI][jI], pthatWeight_);	  

	  if(goodRecoTrunc1p15 && goodGen){
	    recoVsGenJtPtRes1p15_h[rI]->Fill(recoJtPtRes1p15_[rI][jI], genJtPt_[rI][jI], pthatWeight_);
	    responseRes1p15_h[rI]->Fill(recoJtPtRes1p15_[rI][jI], genJtPt_[rI][jI], pthatWeight_);
	  }
	  else if(goodGen) responseRes1p15_h[rI]->Miss(genJtPt_[rI][jI], pthatWeight_);
	  else if(goodRecoTrunc1p15) responseRes1p15_h[rI]->Fake(recoJtPtRes1p15_[rI][jI], pthatWeight_);	  
	}
	else{
	  if(goodGen){
	    genJtPt_Parallel_h[rI]->Fill(genJtPt_[rI][jI], pthatWeight_);
	    if(goodRecoFull) genJtPt_GoodRecoFull_Parallel_h[rI]->Fill(genJtPt_[rI][jI], pthatWeight_);
	  }
	  if(goodRecoTrunc){
	    recoJtPt_Parallel_h[rI]->Fill(recoJtPt_[rI][jI], pthatWeight_);
	    if(!goodGen) recoJtPt_NoGen_Parallel_h[rI]->Fill(recoJtPt_[rI][jI], pthatWeight_);
	  }
	}
      }

      for(Int_t jI = 0; jI < nRCFakes[rI]-1; ++jI){
	Double_t fakeJtPt = randGen_p->Gaus(0, rcWidth[rI]);
        bool goodRecoTrunc = fakeJtPt >= jtPtBins[1] && fakeJtPt < jtPtBins[nJtPtBins-1];

	if(!fillPara) recoJtPt_RCFakeCheck_h[rI]->Fill(fakeJtPt, pthatWeight_);
	else recoJtPt_RCFakeCheck_Parallel_h[rI]->Fill(fakeJtPt, pthatWeight_);


	if(!goodRecoTrunc) continue;
	
	if(!fillPara){
	  recoJtPt_RCFake_h[rI]->Fill(fakeJtPt, pthatWeight_);
	  responseFake_h[rI]->Fake(fakeJtPt, pthatWeight_);
	}
	else recoJtPt_RCFake_Parallel_h[rI]->Fill(fakeJtPt, pthatWeight_);
      }
    }    
  }

  inFile_p->Close();
  delete inFile_p;

  outFile_p->cd();

  for(Int_t rI = 0; rI < nRVals; ++rI){

    for(Int_t bI = 0; bI < nBayes;  ++bI){
      doUnfold(&(unfoldedJtPt_h[rI][bI]), recoJtPt_h[rI], response_h[rI], bI+1, std::string("unfoldedJtPt_R" + prettyString(rVals[rI], 1, true) + "_Bayes" + std::to_string(bI+1) + "_h"));
      doUnfold(&(unfoldedJtPtFake_h[rI][bI]), recoJtPt_h[rI], responseFake_h[rI], bI+1, std::string("unfoldedJtPtFake_R" + prettyString(rVals[rI], 1, true) + "_Bayes" + std::to_string(bI+1) + "_h"));
      doUnfold(&(unfoldedJtPtRes0p9_h[rI][bI]), recoJtPt_h[rI], responseRes0p9_h[rI], bI+1, std::string("unfoldedJtPtRes0p9_R" + prettyString(rVals[rI], 1, true) + "_Bayes" + std::to_string(bI+1) + "_h"));
      doUnfold(&(unfoldedJtPtRes0p85_h[rI][bI]), recoJtPt_h[rI], responseRes0p85_h[rI], bI+1, std::string("unfoldedJtPtRes0p85_R" + prettyString(rVals[rI], 1, true) + "_Bayes" + std::to_string(bI+1) + "_h"));
      doUnfold(&(unfoldedJtPtRes1p1_h[rI][bI]), recoJtPt_h[rI], responseRes1p1_h[rI], bI+1, std::string("unfoldedJtPtRes1p1_R" + prettyString(rVals[rI], 1, true) + "_Bayes" + std::to_string(bI+1) + "_h"));
      doUnfold(&(unfoldedJtPtRes1p15_h[rI][bI]), recoJtPt_h[rI], responseRes1p15_h[rI], bI+1, std::string("unfoldedJtPtRes1p15_R" + prettyString(rVals[rI], 1, true) + "_Bayes" + std::to_string(bI+1) + "_h"));

      doUnfold(&(unfoldedJtPt_Parallel_h[rI][bI]), recoJtPt_Parallel_h[rI], response_h[rI], bI+1, std::string("unfoldedJtPt_R" + prettyString(rVals[rI], 1, true) + "_Bayes" + std::to_string(bI+1) + "_Parallel_h"));
      doUnfold(&(unfoldedJtPtFake_Parallel_h[rI][bI]), recoJtPt_Parallel_h[rI], responseFake_h[rI], bI+1, std::string("unfoldedJtPtFake_R" + prettyString(rVals[rI], 1, true) + "_Bayes" + std::to_string(bI+1) + "_Parallel_h"));
      doUnfold(&(unfoldedJtPtRes0p9_Parallel_h[rI][bI]), recoJtPt_Parallel_h[rI], responseRes0p9_h[rI], bI+1, std::string("unfoldedJtPtRes0p9_R" + prettyString(rVals[rI], 1, true) + "_Bayes" + std::to_string(bI+1) + "_Parallel_h"));
      doUnfold(&(unfoldedJtPtRes0p85_Parallel_h[rI][bI]), recoJtPt_Parallel_h[rI], responseRes0p85_h[rI], bI+1, std::string("unfoldedJtPtRes0p85_R" + prettyString(rVals[rI], 1, true) + "_Bayes" + std::to_string(bI+1) + "_Parallel_h"));
      doUnfold(&(unfoldedJtPtRes1p1_Parallel_h[rI][bI]), recoJtPt_Parallel_h[rI], responseRes1p1_h[rI], bI+1, std::string("unfoldedJtPtRes1p1_R" + prettyString(rVals[rI], 1, true) + "_Bayes" + std::to_string(bI+1) + "_Parallel_h"));
      doUnfold(&(unfoldedJtPtRes1p15_Parallel_h[rI][bI]), recoJtPt_Parallel_h[rI], responseRes1p15_h[rI], bI+1, std::string("unfoldedJtPtRes1p15_R" + prettyString(rVals[rI], 1, true) + "_Bayes" + std::to_string(bI+1) + "_Parallel_h"));
    }
  }

  atlasRAA_R0p4_h->Write("", TObject::kOverwrite);
  genJtPt_AtlasBinnedR0p4_h->Scale(1./(Double_t)fills);
  genJtPt_AtlasBinnedR0p4_h->Write("", TObject::kOverwrite);

  for(Int_t lI = 0; lI < nLogFactor+1; ++lI){
    genJtPt_AtlasBinnedR0p4_LogLoss_h[lI]->Scale(1./(Double_t)fills);
    genJtPt_AtlasBinnedR0p4_LogLoss_h[lI]->Write("", TObject::kOverwrite);

    raa_AtlasBinnedR0p4_LogLoss_h[lI]->Divide(genJtPt_AtlasBinnedR0p4_LogLoss_h[lI], genJtPt_AtlasBinnedR0p4_h);   
    raa_AtlasBinnedR0p4_LogLoss_h[lI]->Write("", TObject::kOverwrite);

    //    atlasChi2_h->SetBinContent(lI+1, raa_AtlasBinnedR0p4_LogLoss_h[lI]->Chi2Test(atlasRAA_R0p4_h, "WW CHI2/NDF"));

    Double_t chi2 = 0.;
    for(Int_t bI = 0; bI < atlasRAA_R0p4_h->GetNbinsX(); ++bI){
      chi2 += (atlasRAA_R0p4_h->GetBinContent(bI+1) - raa_AtlasBinnedR0p4_LogLoss_h[lI]->GetBinContent(bI+1))*(atlasRAA_R0p4_h->GetBinContent(bI+1) - raa_AtlasBinnedR0p4_LogLoss_h[lI]->GetBinContent(bI+1))/atlasRAA_R0p4_h->GetBinContent(bI+1);
    }
    Double_t chi2NDF = chi2/(Double_t)atlasRAA_R0p4_h->GetNbinsX();

    atlasChi2_h->SetBinContent(lI+1, chi2NDF);
    atlasChi2_h->SetBinError(lI+1, 0.0);
  }

  atlasChi2_h->Write("", TObject::kOverwrite);

  for(Int_t rI = 0; rI < nRVals; ++rI){
    genJtPt_h[rI]->Scale(1./(Double_t)fills);

    for(Int_t lI = 0; lI < nLogFactor+1; ++lI){
      genJtPt_LogLoss_h[rI][lI]->Scale(1./(Double_t)fills);
    }

    genJtPt_GoodRecoFull_h[rI]->Scale(1./(Double_t)fills);
    recoJtPt_h[rI]->Scale(1./(Double_t)fills);
    recoJtPt_RCFake_h[rI]->Scale(1./(Double_t)fills);
    recoJtPt_RCFakeCheck_h[rI]->Scale(1./(Double_t)fills);
    recoJtPt_NoGen_h[rI]->Scale(1./(Double_t)fills);

    genJtPt_Parallel_h[rI]->Scale(1./(Double_t)paraFills);
    genJtPt_GoodRecoFull_Parallel_h[rI]->Scale(1./(Double_t)paraFills);
    recoJtPt_Parallel_h[rI]->Scale(1./(Double_t)paraFills);
    recoJtPt_RCFake_Parallel_h[rI]->Scale(1./(Double_t)paraFills);
    recoJtPt_RCFakeCheck_Parallel_h[rI]->Scale(1./(Double_t)paraFills);
    recoJtPt_NoGen_Parallel_h[rI]->Scale(1./(Double_t)paraFills);

    genJtPt_h[rI]->Write("", TObject::kOverwrite);

    for(Int_t lI = 0; lI < nLogFactor+1; ++lI){
      genJtPt_LogLoss_h[rI][lI]->Write("", TObject::kOverwrite);
    }

    genJtPt_GoodRecoFull_h[rI]->Write("", TObject::kOverwrite);
    recoJtPt_h[rI]->Write("", TObject::kOverwrite);
    recoJtPt_RCFake_h[rI]->Write("", TObject::kOverwrite);
    recoJtPt_RCFakeCheck_h[rI]->Write("", TObject::kOverwrite);
    recoJtPt_NoGen_h[rI]->Write("", TObject::kOverwrite);
    genJtPt_Parallel_h[rI]->Write("", TObject::kOverwrite);
    genJtPt_GoodRecoFull_Parallel_h[rI]->Write("", TObject::kOverwrite);
    recoJtPt_Parallel_h[rI]->Write("", TObject::kOverwrite);
    recoJtPt_RCFake_Parallel_h[rI]->Write("", TObject::kOverwrite);
    recoJtPt_RCFakeCheck_Parallel_h[rI]->Write("", TObject::kOverwrite);
    recoJtPt_NoGen_Parallel_h[rI]->Write("", TObject::kOverwrite);

    recoVsGenJtPt_h[rI]->Write("", TObject::kOverwrite);
    response_h[rI]->Write("", TObject::kOverwrite);

    recoVsGenJtPtRes0p9_h[rI]->Write("", TObject::kOverwrite);
    responseRes0p9_h[rI]->Write("", TObject::kOverwrite);

    recoVsGenJtPtRes0p85_h[rI]->Write("", TObject::kOverwrite);
    responseRes0p85_h[rI]->Write("", TObject::kOverwrite);

    recoVsGenJtPtRes1p1_h[rI]->Write("", TObject::kOverwrite);
    responseRes1p1_h[rI]->Write("", TObject::kOverwrite);

    recoVsGenJtPtRes1p15_h[rI]->Write("", TObject::kOverwrite);
    responseRes1p15_h[rI]->Write("", TObject::kOverwrite);

    for(Int_t bI = 0; bI < nBayes;  ++bI){
      unfoldedJtPt_h[rI][bI]->Scale(1./(Double_t)fills);
      unfoldedJtPt_h[rI][bI]->Write("", TObject::kOverwrite);

      unfoldedJtPtFake_h[rI][bI]->Scale(1./(Double_t)fills);
      unfoldedJtPtFake_h[rI][bI]->Write("", TObject::kOverwrite);

      unfoldedJtPtRes0p9_h[rI][bI]->Scale(1./(Double_t)fills);
      unfoldedJtPtRes0p9_h[rI][bI]->Write("", TObject::kOverwrite);

      unfoldedJtPtRes0p85_h[rI][bI]->Scale(1./(Double_t)fills);
      unfoldedJtPtRes0p85_h[rI][bI]->Write("", TObject::kOverwrite);

      unfoldedJtPtRes1p1_h[rI][bI]->Scale(1./(Double_t)fills);
      unfoldedJtPtRes1p1_h[rI][bI]->Write("", TObject::kOverwrite);

      unfoldedJtPtRes1p15_h[rI][bI]->Scale(1./(Double_t)fills);
      unfoldedJtPtRes1p15_h[rI][bI]->Write("", TObject::kOverwrite);

      delete unfoldedJtPt_h[rI][bI];
      delete unfoldedJtPtFake_h[rI][bI];
      delete unfoldedJtPtRes0p9_h[rI][bI];
      delete unfoldedJtPtRes0p85_h[rI][bI];
      delete unfoldedJtPtRes1p1_h[rI][bI];
      delete unfoldedJtPtRes1p15_h[rI][bI];
    }

    for(Int_t bI = 0; bI < nBayes;  ++bI){
      unfoldedJtPt_Parallel_h[rI][bI]->Scale(1./(Double_t)paraFills);
      unfoldedJtPt_Parallel_h[rI][bI]->Write("", TObject::kOverwrite);

      unfoldedJtPtFake_Parallel_h[rI][bI]->Scale(1./(Double_t)paraFills);
      unfoldedJtPtFake_Parallel_h[rI][bI]->Write("", TObject::kOverwrite);

      unfoldedJtPtRes0p9_Parallel_h[rI][bI]->Scale(1./(Double_t)fills);
      unfoldedJtPtRes0p9_Parallel_h[rI][bI]->Write("", TObject::kOverwrite);

      unfoldedJtPtRes0p85_Parallel_h[rI][bI]->Scale(1./(Double_t)fills);
      unfoldedJtPtRes0p85_Parallel_h[rI][bI]->Write("", TObject::kOverwrite);

      unfoldedJtPtRes1p1_Parallel_h[rI][bI]->Scale(1./(Double_t)fills);
      unfoldedJtPtRes1p1_Parallel_h[rI][bI]->Write("", TObject::kOverwrite);

      unfoldedJtPtRes1p15_Parallel_h[rI][bI]->Scale(1./(Double_t)fills);
      unfoldedJtPtRes1p15_Parallel_h[rI][bI]->Write("", TObject::kOverwrite);
     
      delete unfoldedJtPt_Parallel_h[rI][bI];
      delete unfoldedJtPtFake_Parallel_h[rI][bI];
      delete unfoldedJtPtRes0p9_Parallel_h[rI][bI];
      delete unfoldedJtPtRes0p85_Parallel_h[rI][bI];
      delete unfoldedJtPtRes1p1_Parallel_h[rI][bI];
      delete unfoldedJtPtRes1p15_Parallel_h[rI][bI];
    }

    delete genJtPt_h[rI];


    for(Int_t lI = 0; lI < nLogFactor+1; ++lI){
      delete genJtPt_LogLoss_h[rI][lI];
    }

    delete genJtPt_GoodRecoFull_h[rI];
    delete recoJtPt_h[rI];
    delete recoJtPt_RCFake_h[rI];
    delete recoJtPt_RCFakeCheck_h[rI];
    delete recoJtPt_NoGen_h[rI];
    delete genJtPt_Parallel_h[rI];
    delete genJtPt_GoodRecoFull_Parallel_h[rI];
    delete recoJtPt_Parallel_h[rI];
    delete recoJtPt_RCFake_Parallel_h[rI];
    delete recoJtPt_RCFakeCheck_Parallel_h[rI];
    delete recoJtPt_NoGen_Parallel_h[rI];
    delete recoVsGenJtPt_h[rI];
    delete response_h[rI]; 
    delete responseFake_h[rI]; 
    delete recoVsGenJtPtRes0p9_h[rI];
    delete recoVsGenJtPtRes0p85_h[rI];
    delete responseRes0p9_h[rI]; 
    delete responseRes0p85_h[rI]; 
    delete recoVsGenJtPtRes1p1_h[rI];
    delete recoVsGenJtPtRes1p15_h[rI];
    delete responseRes1p1_h[rI]; 
    delete responseRes1p15_h[rI]; 
  }

  delete atlasRAA_R0p4_h;
  delete atlasChi2_h;

  delete genJtPt_AtlasBinnedR0p4_h;
  for(Int_t lI = 0; lI < nLogFactor+1; ++lI){
    delete genJtPt_AtlasBinnedR0p4_LogLoss_h[lI];
    delete raa_AtlasBinnedR0p4_LogLoss_h[lI];
  }

  outFile_p->Close();
  delete outFile_p;  

  delete randGen_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/makeRawRAASpectra.exe <inFileName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += makeRawRAASpectra(argv[1]);
  return retVal;
}
