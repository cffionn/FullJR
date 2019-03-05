#ifndef MACROHISTTOSUBSETHIST_H
#define MACROHISTTOSUBSETHIST_H

#include <vector>

#include "TH1D.h"
#include "TMath.h"
#include "TRandom3.h"

bool macroHistToSubsetHist(TH1D* macroHist_p, TH1D* subsetHist_p, bool doSumW2 = false)
{
  std::vector<double> macroBins, subsetBins;

  for(Int_t bIX = 0; bIX < macroHist_p->GetNbinsX()+1; ++bIX){
    macroBins.push_back(macroHist_p->GetBinLowEdge(bIX+1));
  }
  for(Int_t bIX = 0; bIX < subsetHist_p->GetNbinsX()+1; ++bIX){
    subsetBins.push_back(subsetHist_p->GetBinLowEdge(bIX+1));
  }

  //Check that subset bins Delta > 1
  bool allSubsetLargeDelta = true;
  for(unsigned int sI = 0; sI < subsetBins.size()-1; ++sI){
    if(subsetBins[sI+1] - subsetBins[sI] < 1){
      allSubsetLargeDelta = false;
      break;
    } 
  }

  if(!allSubsetLargeDelta){
    std::cout << "Error in macroHistToSubsetHist: Not all bin boundaries in subset hist have delta > 1" << std::endl;
    std::cout << " Subset hist: ";
    for(auto const & subsetVal : subsetBins){
      std::cout << subsetVal << ",";
    }
    std::cout << std::endl;
    std::cout << "return false" << std::endl;
    return false;
  }

  //Check that macrobins are well contained by subset bins
  bool allBinsGood = true;
  for(auto const & subVal : subsetBins){
    bool binHasMatch = false;
    for(auto const & macroVal : macroBins){
      if(TMath::Abs(macroVal - subVal) < 0.5){
	binHasMatch = true;
	break;
      }
    }

    if(!binHasMatch){
      allBinsGood = false;
      break;
    }
  }

  if(!allBinsGood){
    std::cout << "Error in macroHistToSubsetHist: Not all bin boundaries in subset hist have macro hist match" << std::endl;
    std::cout << " Macro name: " << macroHist_p->GetName() << std::endl;
    std::cout << " subset name: " << subsetHist_p->GetName() << std::endl;
    std::cout << " Macro hist: ";
    for(auto const & macroVal : macroBins){
      std::cout << macroVal << ",";
    }
    std::cout << std::endl;
    std::cout << " Subset hist: ";
    for(auto const & subsetVal : subsetBins){
      std::cout << subsetVal << ",";
    }
    std::cout << std::endl;
    std::cout << "return false" << std::endl;
  }
  else{
    for(unsigned int sI = 0; sI < subsetBins.size(); ++sI){
      Double_t val = 0;
      Double_t err = 0;
      
      for(Int_t bIX = 0; bIX < macroHist_p->GetNbinsX(); ++bIX){
	if(subsetBins[sI] <= macroHist_p->GetBinCenter(bIX+1) && macroHist_p->GetBinCenter(bIX+1) < subsetBins[sI+1]){
	  val += macroHist_p->GetBinContent(bIX+1);
	  err = TMath::Sqrt(err*err + macroHist_p->GetBinError(bIX+1)*macroHist_p->GetBinError(bIX+1));
	}
      }

      subsetHist_p->SetBinContent(sI+1, val);
      if(doSumW2) subsetHist_p->SetBinError(sI+1, err);
      else subsetHist_p->SetBinError(sI+1, TMath::Sqrt(val));
    }

  }
 
  return allBinsGood;
}


void test()
{
  TRandom3* randGen_p = new TRandom3(0);

  TH1D* test1_h = new TH1D("test1_h", "", 10, 0, 10);
  TH1D* test2_h = new TH1D("test2_h", "", 5, 0, 10);

  for(unsigned int i = 0; i < 500; ++i){
    Double_t val = randGen_p->Gaus(5, 1.0);

    test1_h->Fill(val);
  }

  macroHistToSubsetHist(test1_h, test2_h);


  test1_h->Print("ALL");
  test2_h->Print("ALL");

  delete test1_h;
  delete test2_h;

  delete randGen_p;
      
  return;
}

#endif
