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
#include "Utility/include/plotUtilities.h"
#include "Utility/include/histDefUtility.h"

int makeRawRAASpectra(const std::string inFileName)
{
  TRandom3* randGen_p = new TRandom3(0);

  if(!checkFile(inFileName)){
    std::cout << "Given inFileName \'" << inFileName << "\' is invalid. return 1" << std::endl;
    return 1;
  }
  else if(inFileName.find(".root") == std::string::npos){
    std::cout << "Given inFileName \'" << inFileName << "\' does not contain \'.root\'. return 1" << std::endl;
    return 1;
  }

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
  const Int_t nBayes = 6;

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  TH1D* genJtPt_h[nRVals];
  TH1D* genJtPt_GoodRecoFull_h[nRVals];
  TH1D* recoJtPt_h[nRVals];
  TH1D* recoJtPt_NoGen_h[nRVals];
  TH1D* genJtPt_Parallel_h[nRVals];
  TH1D* genJtPt_GoodRecoFull_Parallel_h[nRVals];
  TH1D* recoJtPt_Parallel_h[nRVals];
  TH1D* recoJtPt_NoGen_Parallel_h[nRVals];
  TH2D* recoVsGenJtPt_h[nRVals];
  RooUnfoldResponse* response_h[nRVals];
  RooUnfoldResponse* response0p9Res_h[nRVals];
  RooUnfoldResponse* response1p1Res_h[nRVals];
  TH1D* unfoldedJtPt_h[nRVals][nBayes];
  TH1D* unfoldedJtPt_Parallel_h[nRVals][nBayes];
  TH1D* unfoldedJtPt_0p9Res_h[nRVals][nBayes];
  TH1D* unfoldedJtPt_Parallel_0p9Res_h[nRVals][nBayes];
  TH1D* unfoldedJtPt_1p1Res_h[nRVals][nBayes];
  TH1D* unfoldedJtPt_Parallel_1p1Res_h[nRVals][nBayes];

  const Int_t nJtPtBins = 12;
  Float_t jtPtBins[nJtPtBins+1] = {100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700};
  
  for(Int_t rI = 0; rI < nRVals; ++rI){
    genJtPt_h[rI] = new TH1D(("genJtPt_R" + prettyString(rVals[rI], 1, true) + "_h").c_str(), ";Gen Jet p_{T};Counts", nJtPtBins, jtPtBins);
    genJtPt_GoodRecoFull_h[rI] = new TH1D(("genJtPt_GoodRecoFull_R" + prettyString(rVals[rI], 1, true) + "_h").c_str(), ";Gen Jet p_{T};Counts", nJtPtBins, jtPtBins);
    recoJtPt_h[rI] = new TH1D(("recoJtPt_R" + prettyString(rVals[rI], 1, true) + "_h").c_str(), ";Reco Jet p_{T};Counts", nJtPtBins, jtPtBins);
    recoJtPt_NoGen_h[rI] = new TH1D(("recoJtPt_R" + prettyString(rVals[rI], 1, true) + "_NoGen_h").c_str(), ";Reco Jet p_{T};Counts", nJtPtBins, jtPtBins);

    genJtPt_Parallel_h[rI] = new TH1D(("genJtPt_R" + prettyString(rVals[rI], 1, true) + "_Parallel_h").c_str(), ";Gen Jet p_{T};Counts", nJtPtBins, jtPtBins);
    genJtPt_GoodRecoFull_Parallel_h[rI] = new TH1D(("genJtPt_GoodRecoFull_R" + prettyString(rVals[rI], 1, true) + "_Parallel_h").c_str(), ";Gen Jet p_{T};Counts", nJtPtBins, jtPtBins);
    recoJtPt_Parallel_h[rI] = new TH1D(("recoJtPt_R" + prettyString(rVals[rI], 1, true) + "_Parallel_h").c_str(), ";Reco Jet p_{T};Counts", nJtPtBins, jtPtBins);
    recoJtPt_NoGen_Parallel_h[rI] = new TH1D(("recoJtPt_NoGen_R" + prettyString(rVals[rI], 1, true) + "_Parallel_h").c_str(), ";Reco Jet p_{T};Counts", nJtPtBins, jtPtBins);

    recoVsGenJtPt_h[rI] = new TH2D(("recoVsGenJtPt_R" + prettyString(rVals[rI], 1, true) + "_h").c_str(), ";Reco Jet p_{T};Gen. Jet p_{T}", nJtPtBins, jtPtBins, nJtPtBins, jtPtBins);

    response_h[rI] = new RooUnfoldResponse(("response_" + prettyString(rVals[rI], 1, true)).c_str(), "");
    response_h[rI]->Setup(recoJtPt_h[rI], genJtPt_h[rI]);
    
    for(Int_t bI = 0; bI < nBayes; ++bI){
      unfoldedJtPt_h[rI][bI] = NULL;
      unfoldedJtPt_Parallel_h[rI][bI] = NULL;
    }

    setSumW2({genJtPt_h[rI], genJtPt_GoodRecoFull_h[rI], recoJtPt_h[rI], recoJtPt_NoGen_h[rI], genJtPt_Parallel_h[rI], genJtPt_GoodRecoFull_Parallel_h[rI], recoJtPt_Parallel_h[rI], recoJtPt_NoGen_Parallel_h[rI]});
    centerTitles({genJtPt_h[rI], genJtPt_GoodRecoFull_h[rI], recoJtPt_h[rI], recoJtPt_NoGen_h[rI], genJtPt_Parallel_h[rI], genJtPt_GoodRecoFull_Parallel_h[rI], recoJtPt_Parallel_h[rI], recoJtPt_NoGen_Parallel_h[rI]});

    setSumW2(recoVsGenJtPt_h[rI]);
    centerTitles(recoVsGenJtPt_h[rI]);
  }

  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  TTree* inTree_p[nRVals];

  Float_t pthat_[nRVals];
  Float_t pthatWeight_[nRVals];

  const Int_t nJtMax = 500;
  Int_t nJt_[nRVals];
  Float_t genJtPt_[nRVals][nJtMax];
  Float_t genJtPhi_[nRVals][nJtMax];
  Float_t genJtEta_[nRVals][nJtMax];
  Float_t recoJtPt_[nRVals][nJtMax];
  
  for(Int_t rI = 0; rI < nRVals; ++rI){
    inTree_p[rI] = (TTree*)inFile_p->Get(("recoAndGenTreeR" + prettyString(rVals[rI], 1, true)).c_str());
    
    inTree_p[rI]->SetBranchAddress("pthat", &(pthat_[rI]));
    inTree_p[rI]->SetBranchAddress("pthatWeight", &(pthatWeight_[rI]));
    inTree_p[rI]->SetBranchAddress("nJt", &(nJt_[rI]));
    inTree_p[rI]->SetBranchAddress("genJtPt", genJtPt_[rI]);
    inTree_p[rI]->SetBranchAddress("genJtPhi", genJtPhi_[rI]);
    inTree_p[rI]->SetBranchAddress("genJtEta", genJtEta_[rI]);
    inTree_p[rI]->SetBranchAddress("recoJtPt", recoJtPt_[rI]);
  }

  const Int_t nEntries = TMath::Min((Int_t)100000000, (Int_t)inTree_p[0]->GetEntries());
  Int_t fills = 0;
  Int_t paraFills = 0;

  std::cout << "Processing Entries: " << nEntries << "..." << std::endl;
  for(Int_t entry = 0; entry < nEntries; ++entry){
    if(entry%10000 == 0) std::cout << " Entry: " << entry << "/" << nEntries << std::endl;

    for(Int_t rI = 0; rI < nRVals; ++rI){
      inTree_p[rI]->GetEntry(entry);

      bool fillPara = randGen_p->Uniform(0.0, 1.0) >= 0.5;

      if(fillPara) paraFills++;
      else fills++;

      for(Int_t jI = 0; jI < nJt_[rI]; ++jI){
	bool goodRecoFull = recoJtPt_[rI][jI] >= jtPtBins[0] && recoJtPt_[rI][jI] < jtPtBins[nJtPtBins];
	bool goodRecoTrunc = recoJtPt_[rI][jI] >= jtPtBins[1] && recoJtPt_[rI][jI] < jtPtBins[nJtPtBins-1];
	bool goodGen = genJtPt_[rI][jI] >= jtPtBins[0] && genJtPt_[rI][jI] < jtPtBins[nJtPtBins];

	if(!fillPara){
	  if(goodGen){
	    genJtPt_h[rI]->Fill(genJtPt_[rI][jI], pthatWeight_[rI]);
	    if(goodRecoFull) genJtPt_GoodRecoFull_h[rI]->Fill(genJtPt_[rI][jI], pthatWeight_[rI]);
	  }
	  if(goodRecoTrunc){
	    recoJtPt_h[rI]->Fill(recoJtPt_[rI][jI], pthatWeight_[rI]);
	    if(!goodGen) recoJtPt_NoGen_h[rI]->Fill(recoJtPt_[rI][jI], pthatWeight_[rI]);
	  }
	  
	  if(goodRecoTrunc && goodGen){
	    recoVsGenJtPt_h[rI]->Fill(recoJtPt_[rI][jI], genJtPt_[rI][jI], pthatWeight_[rI]);
	    response_h[rI]->Fill(recoJtPt_[rI][jI], genJtPt_[rI][jI], pthatWeight_[rI]);
	  }
	  else if(goodGen) response_h[rI]->Miss(genJtPt_[rI][jI], pthatWeight_[rI]);
	  else if(goodRecoTrunc) response_h[rI]->Fake(recoJtPt_[rI][jI], pthatWeight_[rI]);
	}
	else{
	  if(goodGen){
	    genJtPt_Parallel_h[rI]->Fill(genJtPt_[rI][jI], pthatWeight_[rI]);
	    if(goodRecoFull) genJtPt_GoodRecoFull_Parallel_h[rI]->Fill(genJtPt_[rI][jI], pthatWeight_[rI]);
	  }
	  if(goodRecoTrunc){
	    recoJtPt_Parallel_h[rI]->Fill(recoJtPt_[rI][jI], pthatWeight_[rI]);
	    if(!goodGen) recoJtPt_NoGen_Parallel_h[rI]->Fill(recoJtPt_[rI][jI], pthatWeight_[rI]);
	  }
	}
      }
    }    
  }

  inFile_p->Close();
  delete inFile_p;

  outFile_p->cd();

  for(Int_t rI = 0; rI < nRVals; ++rI){

    for(Int_t bI = 0; bI < nBayes;  ++bI){
      RooUnfoldBayes* temp = new RooUnfoldBayes(response_h[rI], recoJtPt_h[rI], bI+1);
      temp->SetVerbose(0);
      unfoldedJtPt_h[rI][bI] = (TH1D*)temp->Hreco(RooUnfold::kCovToy)->Clone(("unfoldedJtPt_R" + prettyString(rVals[rI], 1, true) + "_Bayes" + std::to_string(bI+1) + "_h").c_str());
      centerTitles(unfoldedJtPt_h[rI][bI]);
      setSumW2(unfoldedJtPt_h[rI][bI]);

      delete temp;
    }

    for(Int_t bI = 0; bI < nBayes;  ++bI){
      RooUnfoldBayes* temp = new RooUnfoldBayes(response_h[rI], recoJtPt_Parallel_h[rI], bI+1);
      temp->SetVerbose(0);
      unfoldedJtPt_Parallel_h[rI][bI] = (TH1D*)temp->Hreco(RooUnfold::kCovToy)->Clone(("unfoldedJtPt_R" + prettyString(rVals[rI], 1, true) + "_Bayes" + std::to_string(bI+1) + "_Parallel_h").c_str());
      centerTitles(unfoldedJtPt_Parallel_h[rI][bI]);
      setSumW2(unfoldedJtPt_Parallel_h[rI][bI]);

      delete temp;
    }

    genJtPt_h[rI]->Scale(1./(Double_t)fills);
    genJtPt_GoodRecoFull_h[rI]->Scale(1./(Double_t)fills);
    recoJtPt_h[rI]->Scale(1./(Double_t)fills);
    recoJtPt_NoGen_h[rI]->Scale(1./(Double_t)fills);

    genJtPt_Parallel_h[rI]->Scale(1./(Double_t)paraFills);
    genJtPt_GoodRecoFull_Parallel_h[rI]->Scale(1./(Double_t)paraFills);
    recoJtPt_Parallel_h[rI]->Scale(1./(Double_t)paraFills);
    recoJtPt_NoGen_Parallel_h[rI]->Scale(1./(Double_t)paraFills);

    genJtPt_h[rI]->Write("", TObject::kOverwrite);
    genJtPt_GoodRecoFull_h[rI]->Write("", TObject::kOverwrite);
    recoJtPt_h[rI]->Write("", TObject::kOverwrite);
    recoJtPt_NoGen_h[rI]->Write("", TObject::kOverwrite);
    genJtPt_Parallel_h[rI]->Write("", TObject::kOverwrite);
    genJtPt_GoodRecoFull_Parallel_h[rI]->Write("", TObject::kOverwrite);
    recoJtPt_Parallel_h[rI]->Write("", TObject::kOverwrite);
    recoJtPt_NoGen_Parallel_h[rI]->Write("", TObject::kOverwrite);
    recoVsGenJtPt_h[rI]->Write("", TObject::kOverwrite);
    response_h[rI]->Write("", TObject::kOverwrite);

    for(Int_t bI = 0; bI < nBayes;  ++bI){
      unfoldedJtPt_h[rI][bI]->Scale(1./(Double_t)fills);
      unfoldedJtPt_h[rI][bI]->Write("", TObject::kOverwrite);

      delete unfoldedJtPt_h[rI][bI];
    }

    for(Int_t bI = 0; bI < nBayes;  ++bI){
      unfoldedJtPt_Parallel_h[rI][bI]->Scale(1./(Double_t)paraFills);
      unfoldedJtPt_Parallel_h[rI][bI]->Write("", TObject::kOverwrite);

      delete unfoldedJtPt_Parallel_h[rI][bI];
    }

    delete genJtPt_h[rI];
    delete genJtPt_GoodRecoFull_h[rI];
    delete recoJtPt_h[rI];
    delete recoJtPt_NoGen_h[rI];
    delete genJtPt_Parallel_h[rI];
    delete genJtPt_GoodRecoFull_Parallel_h[rI];
    delete recoJtPt_Parallel_h[rI];
    delete recoJtPt_NoGen_Parallel_h[rI];
    delete recoVsGenJtPt_h[rI];
    delete response_h[rI]; 
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
