//cpp dependencies
#include <iostream>

//ROOT dependencies
#include "TFile.h"
#include "TH2D.h"
#include "TMath.h"
#include "TRandom3.h"

//Non-local dependencies
#include "Utility/include/checkMakeDir.h"
#include "Utility/include/histDefUtility.h"

int comboTest()
{
  TRandom3* randGen_p = new TRandom3();

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  checkMakeDir("output/" + dateStr);
  std::string outFileName = "output/" + dateStr + "/comboTest_" + dateStr + ".root";
  std::string outFileName1 = "output/" + dateStr + "/comboTest1_" + dateStr + ".root";
  std::string outFileName2 = "output/" + dateStr + "/comboTest2_" + dateStr + ".root";

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");

  const Int_t nBins = 5;
  const Double_t lowVal = 0;
  const Double_t highVal = 5;

  TH2D* allWeight_p = new TH2D("allWeight_h", ";X;Y", nBins, lowVal, highVal, nBins, lowVal, highVal);
  TH2D* allWeightErr_p = new TH2D("allWeightErr_h", ";X;Y", nBins, lowVal, highVal, nBins, lowVal, highVal);
  TH2D* oneWeight_p = new TH2D("oneWeight_h", ";X;Y", nBins, lowVal, highVal, nBins, lowVal, highVal);
  TH2D* twoWeight_p = new TH2D("twoWeight_h", ";X;Y", nBins, lowVal, highVal, nBins, lowVal, highVal);
  TH2D* comboWeight_p = new TH2D("comboWeight_h", ";X;Y", nBins, lowVal, highVal, nBins, lowVal, highVal);
  TH2D* comboWeightErr_p = new TH2D("comboWeightErr_h", ";X;Y", nBins, lowVal, highVal, nBins, lowVal, highVal);

  TH2D* ratioWeightErr_p = new TH2D("ratioWeightErr_h", ";X;Y", nBins, lowVal, highVal, nBins, lowVal, highVal);

  std::vector<TH1*> tempVect = {allWeight_p, allWeightErr_p, oneWeight_p, twoWeight_p, comboWeight_p, comboWeightErr_p};
  centerTitles(tempVect);
  setSumW2(tempVect);

  const Double_t weight1 = 0.21;
  const Double_t weight2 = 0.37;
  const Double_t probWeight1 = 0.7;
  const Double_t fillHist1 = 0.5;
  const Int_t nFills = 2000;

  for(Int_t fI = 0; fI < nFills; ++fI){
    Double_t xVal = randGen_p->Uniform(lowVal, highVal);
    Double_t yVal = randGen_p->Uniform(lowVal, highVal);
    Double_t weight = weight2;
    if(randGen_p->Uniform(0, 1) < probWeight1) weight = weight1;

    if(randGen_p->Uniform(0, 1) < fillHist1) oneWeight_p->Fill(xVal, yVal, weight);
    else twoWeight_p->Fill(xVal, yVal, weight);

    allWeight_p->Fill(xVal, yVal, weight);
  }

  for(Int_t bIX = 0; bIX < comboWeight_p->GetXaxis()->GetNbins(); ++bIX){
    for(Int_t bIY = 0; bIY < comboWeight_p->GetYaxis()->GetNbins(); ++bIY){
      Double_t val = oneWeight_p->GetBinContent(bIX+1, bIY+1) + twoWeight_p->GetBinContent(bIX+1, bIY+1);

      Double_t err = TMath::Sqrt(oneWeight_p->GetBinError(bIX+1, bIY+1)*oneWeight_p->GetBinError(bIX+1, bIY+1) + twoWeight_p->GetBinError(bIX+1, bIY+1)*twoWeight_p->GetBinError(bIX+1, bIY+1));

      comboWeight_p->SetBinContent(bIX+1, bIY+1, val);
      comboWeight_p->SetBinError(bIX+1, bIY+1, err);

      comboWeightErr_p->SetBinContent(bIX+1, bIY+1, err);
      comboWeightErr_p->SetBinError(bIX+1, bIY+1, 0);

      allWeightErr_p->SetBinContent(bIX+1, bIY+1, allWeight_p->GetBinError(bIX+1, bIY+1));
      allWeightErr_p->SetBinError(bIX+1, bIY+1, 0);
    }
  }

  outFile_p->cd();
  
  allWeight_p->Write("", TObject::kOverwrite);
  allWeightErr_p->Write("", TObject::kOverwrite);
  oneWeight_p->Write("", TObject::kOverwrite);
  twoWeight_p->Write("", TObject::kOverwrite);
  comboWeight_p->Write("", TObject::kOverwrite);
  comboWeightErr_p->Write("", TObject::kOverwrite);

  ratioWeightErr_p->Divide(allWeightErr_p, comboWeightErr_p);
  ratioWeightErr_p->Write("", TObject::kOverwrite);

  TFile* outFile1_p = new TFile(outFileName1.c_str(), "RECREATE");
  outFile1_p->cd();

  oneWeight_p->Write("", TObject::kOverwrite);

  outFile1_p->Close();
  delete outFile1_p;

  TFile* outFile2_p = new TFile(outFileName2.c_str(), "RECREATE");
  outFile2_p->cd();

  twoWeight_p->Write("oneWeight_h", TObject::kOverwrite);

  outFile2_p->Close();
  delete outFile2_p;

  delete allWeight_p;
  delete allWeightErr_p;
  delete oneWeight_p;
  delete twoWeight_p;
  delete comboWeight_p;
  delete comboWeightErr_p;
  delete ratioWeightErr_p;

  outFile_p->Close();
  delete outFile_p;

  delete randGen_p;

  return 0;
}

int main()
{
  int retVal = 0;
  retVal += comboTest();
  return retVal;
}
