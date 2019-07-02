//cpp
#include <iostream>
#include <string>

//ROOT
#include "TDatime.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TRandom3.h"

//Local dependencies
#include "MainAnalysis/include/macroHistToSubsetHist.h"
#include "Utility/include/checkMakeDir.h"
#include "Utility/include/cppWatch.h"
#include "Utility/include/histDefUtility.h"
#include "Utility/include/plotUtilities.h"

//RooUnfoldResponse
#include "src/RooUnfoldBayes.h"
#include "src/RooUnfoldResponse.h"

void constructResponse(Double_t par0, Double_t par1, Double_t par2, Int_t nMax, TH2D* outResponse_p)
{
  TRandom3* randGen_p = new TRandom3(0);

  Int_t minValX = outResponse_p->GetXaxis()->GetBinLowEdge(1);
  Int_t maxValX = outResponse_p->GetXaxis()->GetBinLowEdge(outResponse_p->GetXaxis()->GetNbins()+1);

  Int_t minValY = outResponse_p->GetYaxis()->GetBinLowEdge(1);
  Int_t maxValY = outResponse_p->GetYaxis()->GetBinLowEdge(outResponse_p->GetYaxis()->GetNbins()+1);

  Double_t minValXD = minValX;
  Double_t maxValXD = maxValX;
  Double_t minValYD = minValY;
  Double_t maxValYD = maxValY;

  //  Int_t nBinsX = maxValX - minValX;
  Int_t nBinsY = maxValY - minValY;
  const Int_t nMaxBins = 2000;
  //  Double_t binsX[nMaxBins];
  Double_t binsY[nMaxBins];

  Double_t counts[nMaxBins];
  Double_t weights[nMaxBins];

  //  for(Int_t bI = 0; bI < nBinsX+1; ++bI){
  //    binsX[bI] = minValX + bI;
  //  }
  for(Int_t bI = 0; bI < nBinsY+1; ++bI){
    binsY[bI] = minValY + bI;
  }
  
  for(Int_t bI = 0; bI < nBinsY; ++bI){
    Double_t yVal = (binsY[bI] + binsY[bI+1])/2.;
    Double_t res = TMath::Sqrt(par0*par0 + par1*par1/yVal + par2*par2/(yVal*yVal));
    Double_t norm = 1./TMath::Sqrt(2.*TMath::Pi()*res*res);

    TF1* gaus_p = new TF1("gaus_p", "gaus", 1 - 10.*res, 1 + 10.*res);
    gaus_p->SetParameters(norm, 1, res);

    std::cout << "GAUS: " << yVal << ", " << res << ", " << norm << std::endl;
    
    weights[bI] = 1.;
    if(1-10.*res <= minValXD/yVal) weights[bI] -= gaus_p->Integral(1 - 10.*res, minValXD/yVal);
    if(1+10.*res >= maxValXD/yVal) weights[bI] -= gaus_p->Integral(maxValXD/yVal, 1 + 10.*res);

    counts[bI] = 0;

    delete gaus_p;
  }

  
  for(Int_t bI = 0; bI < nBinsY; ++bI){
    Double_t yVal = (binsY[bI] + binsY[bI+1])/2.;
    std::cout << "WEIGHTS: " << yVal << ", " << weights[bI] << std::endl;
  }

  Double_t uniformMin = TMath::Power(minValYD/maxValYD, 4);
  
  for(Int_t nI = 0; nI < nMax; ++nI){
    Double_t prob = randGen_p->Uniform(uniformMin, 1.);
    Double_t genPt = minValYD/TMath::Power(prob, 1./4.);

    Int_t pos = (Int_t)(genPt - minValY);
    //    std::cout << pos << ", " << genPt << ", " << prob << ", " << uniformMin << ", " << minVal << std::endl;
    ++counts[pos];
  }

  std::cout << "COUNTS: " << std::endl;
  for(Int_t bI = 0; bI < nBinsY; ++bI){
    //    std::cout << " " << bI << ": " << counts[bI] << std::endl;

    Double_t yVal = (binsY[bI] + binsY[bI+1])/2.;
    Double_t res = TMath::Sqrt(par0*par0 + par1*par1/yVal + par2*par2/(yVal*yVal));
    Double_t norm = 1./TMath::Sqrt(2.*TMath::Pi()*res*res);

    TF1* gaus_p = new TF1("gaus_p", "gaus", minValXD/yVal, maxValXD/yVal);
    gaus_p->SetParameters(norm, 1, res);

    for(Int_t nI = 0; nI < counts[bI]; ++nI){      
      Double_t xVal = yVal*gaus_p->GetRandom(minValXD/yVal, maxValXD/yVal);
      
      outResponse_p->Fill(xVal, yVal, weights[bI]);
    }

    delete gaus_p;
  } 
  delete randGen_p;  

  return;
}


void doUnfold(TH2D* matrix_p, TH1D* prior_p, TH1D* data_p, Int_t nBayes, Int_t bayesVals[], std::string tagStr, TH1D* unfolded_p[])
{
  TH2D* matrixClone_p = (TH2D*)matrix_p->Clone("matrixClone_h");

  const Int_t nMaxBins = 100;
  Double_t recoBins[nMaxBins+1];
  Double_t genBins[nMaxBins+1];
  Int_t nRecoBins = -1;
  Int_t nGenBins = -1;
  
  for(Int_t bIY = 0; bIY < matrixClone_p->GetYaxis()->GetNbins()+1; ++bIY){
    genBins[bIY] = matrixClone_p->GetYaxis()->GetBinLowEdge(bIY+1);
    ++nGenBins;
  }

  for(Int_t bIX = 0; bIX < matrixClone_p->GetXaxis()->GetNbins()+1; ++bIX){
    recoBins[bIX] = matrixClone_p->GetXaxis()->GetBinLowEdge(bIX+1);
    ++nRecoBins;
  }

  TH1D* gen_p = new TH1D("gen_h", "", nGenBins, genBins);
  TH1D* reco_p = new TH1D("reco_h", "", nRecoBins, recoBins);

  macroHistToSubsetHistY(matrixClone_p, gen_p, true);

  if(tagStr.find("5p") != std::string::npos){
    std::cout << "CHECKStart: " << std::endl;
    gen_p->Print("ALL");
  }
  
  for(Int_t bIY = 0; bIY < matrixClone_p->GetYaxis()->GetNbins(); ++bIY){
    Double_t total = 0.0;
    
    for(Int_t bIX = 0; bIX < matrixClone_p->GetXaxis()->GetNbins(); ++bIX){
      total += matrixClone_p->GetBinContent(bIX+1, bIY+1);
    }

    if(total <= TMath::Power(10, -20)) continue;

    Double_t scale = prior_p->GetBinContent(bIY+1)/total;

    for(Int_t bIX = 0; bIX < matrixClone_p->GetXaxis()->GetNbins(); ++bIX){
      Double_t val = matrixClone_p->GetBinContent(bIX+1, bIY+1)*scale;
      Double_t err = matrixClone_p->GetBinError(bIX+1, bIY+1)*scale;

      matrixClone_p->SetBinContent(bIX+1, bIY+1, val);
      matrixClone_p->SetBinError(bIX+1, bIY+1, err);
    }
  }

  macroHistToSubsetHistX(matrixClone_p, reco_p, true);
  macroHistToSubsetHistY(matrixClone_p, gen_p, true);

  if(tagStr.find("5p") != std::string::npos){
    std::cout << "CHECK: " << std::endl;
    prior_p->Print("ALL");
    gen_p->Print("ALL");
  }
  
  RooUnfoldResponse* rooRes_p = NULL;
  RooUnfoldBayes* rooBayes_p = NULL;

  for(Int_t bI = 0; bI < nBayes; ++bI){   
    //    std::cout << "Check integrals: " << reco_p->Integral() << ", " << gen_p->Integral() << ", " << matrixClone_p->Integral() << std::endl;
    
    rooRes_p = new RooUnfoldResponse(reco_p, gen_p, matrixClone_p);
    rooBayes_p = new RooUnfoldBayes(rooRes_p, data_p, bayesVals[bI], false, "temp");
    rooBayes_p->SetVerbose(-1);
    rooBayes_p->SetNToys(1000);

    TH1D* tempUnfold_h = (TH1D*)rooBayes_p->Hreco(RooUnfold::kCovToy);

    unfolded_p[bI] = new TH1D((tagStr + "_unfoldedBayes" + std::to_string(bayesVals[bI]) + "_h").c_str(), ";pT;Counts", nGenBins, genBins);

    for(Int_t tI = 0; tI < tempUnfold_h->GetNbinsX(); ++tI){
      unfolded_p[bI]->SetBinContent(tI+1, tempUnfold_h->GetBinContent(tI+1));
      unfolded_p[bI]->SetBinError(tI+1, tempUnfold_h->GetBinError(tI+1));
    }

    delete tempUnfold_h;
    delete rooBayes_p;
    delete rooRes_p;
  }   

  delete gen_p;
  delete reco_p;
  delete matrixClone_p;

  return;
}

int toyATLASCMS(int seed, int number)  
{
  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  TRandom3* randGen_p = new TRandom3(seed);
  
  const Double_t atlasPar0 = 0.0552816;
  const Double_t atlasPar1 = 0.501155;
  const Double_t atlasPar2 = 14.5295;

  const Double_t cmsPar0 = 0.0443812;
  const Double_t cmsPar1 = 1.49998;
  const Double_t cmsPar2 = 14.1569;

  const Int_t nGenBins = 11;
  Double_t genBins[nGenBins+1] = {20,80,100,160,220,280,350,420,500,630,800,1400};

  /*
  const Int_t nRecoBins = 25;
  Double_t recoBins[nRecoBins+1] = {80,90,100,110,120,140,160,180,200,220,250,280,310,350,400,450,500,560,630,710,800,890,1000,1120,1260,1410};
  */
  const Int_t nRecoBins = 9;
  Double_t recoBins[nRecoBins+1] = {80,100,160,220,280,350,420,500,630,800};

  const Int_t nPowers = 4;
  const Double_t powers[nPowers] = {2.,3.,4.,5.};

  TH2D* atlasResponse_h[nPowers];
  TH1D* atlasResReco_h[nPowers];
  TH1D* atlasResGen_h[nPowers];

  TH2D* cmsResponse_h[nPowers];
  TH1D* cmsResReco_h[nPowers];
  TH1D* cmsResGen_h[nPowers];
  
  TH1D* genSpectraPower_p[nPowers];
  TH1D* genSpectraATLASTruncPower_p[nPowers];
  TH1D* recoATLASSpectraPower_p[nPowers];
  TH1D* genSpectraCMSTruncPower_p[nPowers];
  TH1D* recoCMSSpectraPower_p[nPowers];

  for(Int_t pI = 0; pI < nPowers; ++pI){
    std::string powStr = "Power" + prettyString(powers[pI], 1, true);
    
    atlasResponse_h[pI] = new TH2D(("atlasResponse" + powStr  + "_h").c_str(), ";Reco;Gen", nRecoBins, recoBins, nGenBins, genBins);
    atlasResReco_h[pI] = new TH1D(("atlasResReco" + powStr  + "_h").c_str(), ";Reco;Counts", nRecoBins, recoBins);
    atlasResGen_h[pI] = new TH1D(("atlasResGen" + powStr  + "_h").c_str(), ";Gen;Counts", nGenBins, genBins);
    
    cmsResponse_h[pI] = new TH2D(("cmsResponse" + powStr  + "_h").c_str(), ";Reco;Gen", nRecoBins, recoBins, nGenBins, genBins);
    cmsResReco_h[pI] = new TH1D(("cmsResReco" + powStr  + "_h").c_str(), ";Reco;Counts", nRecoBins, recoBins);
    cmsResGen_h[pI] = new TH1D(("cmsResGen" + powStr  + "_h").c_str(), ";Gen;Counts", nGenBins, genBins);
    
    genSpectraPower_p[pI] = new TH1D(("genSpectra" + powStr + "_h").c_str(), ";gen;counts", nGenBins, genBins);
    genSpectraATLASTruncPower_p[pI] = new TH1D(("genSpectraATLASTrunc" + powStr + "_h").c_str(), ";gen;counts", nGenBins, genBins);
    recoATLASSpectraPower_p[pI] = new TH1D(("recoATLASSpectra" + powStr + "_h").c_str(), ";recoATLAS;counts", nRecoBins, recoBins);
    genSpectraCMSTruncPower_p[pI] = new TH1D(("genSpectraCMSTrunc" + powStr + "_h").c_str(), ";gen;counts", nGenBins, genBins);
    recoCMSSpectraPower_p[pI] = new TH1D(("recoCMSSpectra" + powStr + "_h").c_str(), ";recoCMS;counts", nRecoBins, recoBins);

    setSumW2({atlasResponse_h[pI], atlasResReco_h[pI], atlasResGen_h[pI], cmsResponse_h[pI], cmsResReco_h[pI], cmsResGen_h[pI], genSpectraPower_p[pI], genSpectraATLASTruncPower_p[pI], recoATLASSpectraPower_p[pI], genSpectraCMSTruncPower_p[pI], recoCMSSpectraPower_p[pI]});
  }

  ULong64_t nFill = (ULong64_t)number;
  Double_t fillWeight = 0.01;

  cppWatch response1, response2;

  response1.start();

  for(Int_t pI = 0; pI < nPowers; ++pI){
    for(ULong64_t gI = 0; gI < nFill; ++gI){
      Double_t prob = randGen_p->Uniform(0, 1);
      Double_t genPt = 20./TMath::Power(prob, 1./powers[pI]);
      Double_t genPtToFill = genPt;
    
      //    if(genPt < genBins[0]) genPtToFill = genBins[0]+1;
      //    if(genPt >= genBins[nGenBins]) genPtToFill = genBins[nGenBins]-1;
      
      if(genPt < genBins[0]) continue;
      if(genPt >= genBins[nGenBins]) continue;
      
      Double_t resATLAS = TMath::Sqrt(atlasPar0*atlasPar0 + atlasPar1*atlasPar1/genPt + atlasPar2*atlasPar2/(genPt*genPt));
      Double_t resCMS = TMath::Sqrt(cmsPar0*cmsPar0 + cmsPar1*cmsPar1/genPt + cmsPar2*cmsPar2/(genPt*genPt));
      
      Double_t recoATLAS = genPt*randGen_p->Gaus(1., resATLAS);
      Double_t recoCMS = genPt*randGen_p->Gaus(1., resCMS);
      
      if(recoATLAS >= recoBins[0] && recoATLAS < recoBins[nRecoBins]){
	atlasResponse_h[pI]->Fill(recoATLAS, genPtToFill, fillWeight);
	atlasResReco_h[pI]->Fill(recoATLAS, fillWeight);
	atlasResGen_h[pI]->Fill(genPtToFill, fillWeight);
      }
      
      if(recoCMS >= recoBins[0] && recoCMS < recoBins[nRecoBins]){
	cmsResponse_h[pI]->Fill(recoCMS, genPtToFill, fillWeight);
	cmsResReco_h[pI]->Fill(recoCMS, fillWeight);
	cmsResGen_h[pI]->Fill(genPtToFill, fillWeight);
      }
    }
  }
  
  response1.stop();

  /*
  response2.start();
  constructResponse(atlasPar0, atlasPar1, atlasPar2, nFill, atlas2Response_h);
  constructResponse(cmsPar0, cmsPar1, cmsPar2, nFill, cms2Response_h);
  response2.stop();
  */
  std::cout << "TIMING: " << std::endl;
  std::cout << " First construction: " << response1.totalWall() << ", " << response1.totalCPU() << std::endl;
  std::cout << " Second construction: " << response2.totalWall() << ", " << response2.totalCPU() << std::endl;

  for(Int_t pI = 0; pI < nPowers; ++pI){
    for(ULong64_t gI = 0; gI < nFill; ++gI){
      Double_t prob = randGen_p->Uniform(0, 1);
      Double_t genPt = 20./TMath::Power(prob, 1./powers[pI]);
      Double_t genPtToFill = genPt;
      
      //    if(genPt < genBins[0]) genPtToFill = genBins[0]+1;
      //    if(genPt >= genBins[nGenBins]) genPtToFill = genBins[nGenBins]-1;
      
      if(genPt < genBins[0]) continue;
      if(genPt >= genBins[nGenBins]) continue;
      
      genSpectraPower_p[pI]->Fill(genPtToFill, fillWeight);
      
      Double_t resATLAS = TMath::Sqrt(atlasPar0*atlasPar0 + atlasPar1*atlasPar1/genPt + atlasPar2*atlasPar2/(genPt*genPt));
      Double_t resCMS = TMath::Sqrt(cmsPar0*cmsPar0 + cmsPar1*cmsPar1/genPt + cmsPar2*cmsPar2/(genPt*genPt));
      
      Double_t recoATLAS = genPt*randGen_p->Gaus(1., resATLAS);
      Double_t recoCMS = genPt*randGen_p->Gaus(1., resCMS);
      
      if(recoATLAS >= recoBins[0] && recoATLAS < recoBins[nRecoBins]){
	recoATLASSpectraPower_p[pI]->Fill(recoATLAS, fillWeight);
	genSpectraATLASTruncPower_p[pI]->Fill(genPtToFill, fillWeight);
      }
      
      if(recoCMS >= recoBins[0] && recoCMS < recoBins[nRecoBins]){
	recoCMSSpectraPower_p[pI]->Fill(recoCMS, fillWeight);
	genSpectraCMSTruncPower_p[pI]->Fill(genPtToFill, fillWeight);
      }      
    }
  }
  
  
  std::cout << "STARTED UNFOLDING" << std::endl;
  //Unfolding testing

  checkMakeDir("output");
  checkMakeDir("output/" + dateStr);
  
  TFile* outFile_p = new TFile(("output/" + dateStr + "/toyCMSATLAS_NFill" + std::to_string(nFill) + "_Seed" + std::to_string(seed) + "_" + dateStr + ".root").c_str(), "RECREATE");

  for(Int_t pI = 0; pI < nPowers; ++pI){
    atlasResponse_h[pI]->Write("", TObject::kOverwrite);
    atlasResReco_h[pI]->Write("", TObject::kOverwrite);
    atlasResGen_h[pI]->Write("", TObject::kOverwrite);
    cmsResponse_h[pI]->Write("", TObject::kOverwrite);
    cmsResReco_h[pI]->Write("", TObject::kOverwrite);
    cmsResGen_h[pI]->Write("", TObject::kOverwrite);
    
    genSpectraPower_p[pI]->Write("", TObject::kOverwrite);
    genSpectraATLASTruncPower_p[pI]->Write("", TObject::kOverwrite);
    recoATLASSpectraPower_p[pI]->Write("", TObject::kOverwrite);
    genSpectraCMSTruncPower_p[pI]->Write("", TObject::kOverwrite);
    recoCMSSpectraPower_p[pI]->Write("", TObject::kOverwrite);
  }
  

  for(Int_t pI = 0; pI < nPowers; ++pI){
    delete atlasResponse_h[pI];
    delete atlasResReco_h[pI];
    delete atlasResGen_h[pI];
    delete cmsResponse_h[pI];
    delete cmsResReco_h[pI];
    delete cmsResGen_h[pI];
    
    delete genSpectraPower_p[pI];
    delete genSpectraATLASTruncPower_p[pI];
    delete recoATLASSpectraPower_p[pI];
    delete genSpectraCMSTruncPower_p[pI];
    delete recoCMSSpectraPower_p[pI];
  }
  
  outFile_p->Close();
  delete outFile_p;

  delete randGen_p;
  
  return 0;
}


int main(int argc, char* argv[])
{
  if(argc != 3){
    std::cout << "Usage: ./bin/toyATLASCMS.exe <seed> <number>" << std::endl;
    return 1;
  }
  
  int retVal = 0;
  retVal += toyATLASCMS(std::stoi(argv[1]), std::stoi(argv[2]));
  return retVal;
}
