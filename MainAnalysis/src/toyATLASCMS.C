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

int toyATLASCMS(int seed, ULong64_t number)  
{
  cppWatch total;
  total.start();

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  TRandom3* randGen_p = new TRandom3(seed);
  
  const Double_t atlasPar0 = 0.0552816;
  const Double_t atlasPar1 = 0.501155;
  const Double_t atlasPar2 = 14.5295;

  const Double_t cmsR4Par0 = 0.0443812;
  const Double_t cmsR4Par1 = 1.49998;
  const Double_t cmsR4Par2 = 14.1569;

  const Double_t cmsPPR4Par0 = 0.0443812;
  const Double_t cmsPPR4Par1 = 1.49998;
  const Double_t cmsPPR4Par2 = 0.0;

  const Double_t cmsR10Par0 = 0.0749348;
  const Double_t cmsR10Par1 = 1.5;
  const Double_t cmsR10Par2 = 30.;
  /*
  const Int_t nGenBins = 11;
  Double_t genBins[nGenBins+1] = {20,80,100,160,220,280,350,420,500,630,800,1400};
  */
  const Int_t nMaxBins = 1000;
  /*
  const Int_t nGenBins = 30;
  Double_t genBins[nGenBins+1] = {20,30,40,50,60,70,80,90,100,110,120,140,160,180,200,220,250,280,310,350,400,450,500,560,630,710,800,890,1000,1120,1260};
  */
  Int_t lowGenVal = 10;
  Int_t hiGenVal = 1500;
  Int_t nGenBins = (hiGenVal - lowGenVal)/5.;
  Double_t genBins[nMaxBins+1];
  for(Int_t gI = 0; gI < nGenBins+1; ++gI){
    genBins[gI] = lowGenVal + 5.*gI;
  }

  /*
  const Int_t nRecoBins = 25;
  Double_t recoBins[nRecoBins+1] = {80,90,100,110,120,140,160,180,200,220,250,280,310,350,400,450,500,560,630,710,800,890,1000,1120,1260,1410};
  */
  Int_t lowRecoVal = 40;
  Int_t hiRecoVal = 1500;
  Int_t nRecoBins = (hiRecoVal - lowRecoVal)/5.;
  Double_t recoBins[nMaxBins+1];
  for(Int_t gI = 0; gI < nRecoBins+1; ++gI){
    recoBins[gI] = lowRecoVal + 5.*gI;
  }
  
  /*
  const Int_t nRecoBins = 9;
  Double_t recoBins[nRecoBins+1] = {80,100,160,220,280,350,420,500,630,800};
  */
  const Int_t nPowers = 4;
  const Double_t powers[nPowers] = {6, 7, 8 ,9};

  const Int_t nRandomGen = 5;
  const Int_t nTarget = 3;
  const Double_t targetPoint = 500.;
  const Double_t errTarget = 200./20.; 

  TH2D* atlasResponse_h[nPowers];
  TH1D* atlasResReco_h[nPowers];
  TH1D* atlasResGen_h[nPowers];

  TH2D* cmsR4Response_h[nPowers];
  TH1D* cmsR4ResReco_h[nPowers];
  TH1D* cmsR4ResGen_h[nPowers];
  TH2D* cmsPPR4Response_h[nPowers];
  TH1D* cmsPPR4ResReco_h[nPowers];
  TH1D* cmsPPR4ResGen_h[nPowers];
  TH2D* cmsR10Response_h[nPowers];
  TH1D* cmsR10ResReco_h[nPowers];
  TH1D* cmsR10ResGen_h[nPowers];
  TH2D* cmsR10SmearResponse_h[nPowers];
  TH1D* cmsR10SmearResReco_h[nPowers];
  TH1D* cmsR10SmearResGen_h[nPowers];
  
  TH1D* genSpectraPower_p[nPowers];
  TH1D* genSpectraATLASTruncPower_p[nPowers];
  TH1D* recoATLASSpectraPower_p[nPowers];
  TH1D* genSpectraCMSR4TruncPower_p[nPowers];
  TH1D* recoCMSR4SpectraPower_p[nPowers];
  TH1D* genSpectraCMSPPR4TruncPower_p[nPowers];
  TH1D* recoCMSPPR4SpectraPower_p[nPowers];
  TH1D* genSpectraCMSR10TruncPower_p[nPowers];
  TH1D* recoCMSR10SpectraPower_p[nPowers];
  TH1D* genSpectraCMSR10SmearTruncPower_p[nPowers];
  TH1D* recoCMSR10SmearSpectraPower_p[nPowers];

  TH2D* recoGenSpectraATLASTruncPower_p[nPowers];
  TH2D* recoGenSpectraCMSR4TruncPower_p[nPowers];
  TH2D* recoGenSpectraCMSPPR4TruncPower_p[nPowers];
  TH2D* recoGenSpectraCMSR10TruncPower_p[nPowers];
  TH2D* recoGenSpectraCMSR10SmearTruncPower_p[nPowers];


  TH1D* genSpectraATLASTruncPower_Rand_p[nPowers][nRandomGen][nTarget];
  TH1D* genSpectraCMSR4TruncPower_Rand_p[nPowers][nRandomGen][nTarget];
  TH1D* genSpectraCMSPPR4TruncPower_Rand_p[nPowers][nRandomGen][nTarget];
  TH1D* genSpectraCMSR10TruncPower_Rand_p[nPowers][nRandomGen][nTarget];
  TH1D* genSpectraCMSR10SmearTruncPower_Rand_p[nPowers][nRandomGen][nTarget];

  TH2D* recoGenSpectraATLASTruncPower_Rand_p[nPowers][nRandomGen][nTarget];
  TH2D* recoGenSpectraCMSR4TruncPower_Rand_p[nPowers][nRandomGen][nTarget];
  TH2D* recoGenSpectraCMSPPR4TruncPower_Rand_p[nPowers][nRandomGen][nTarget];
  TH2D* recoGenSpectraCMSR10TruncPower_Rand_p[nPowers][nRandomGen][nTarget];
  TH2D* recoGenSpectraCMSR10SmearTruncPower_Rand_p[nPowers][nRandomGen][nTarget];

  for(Int_t pI = 0; pI < nPowers; ++pI){
    std::string powStr = "Power" + prettyString(powers[pI], 1, true);
    
    atlasResponse_h[pI] = new TH2D(("atlasResponse" + powStr  + "_h").c_str(), ";Reco;Gen", nRecoBins, recoBins, nGenBins, genBins);
    atlasResReco_h[pI] = new TH1D(("atlasResReco" + powStr  + "_h").c_str(), ";Reco;Counts", nRecoBins, recoBins);
    atlasResGen_h[pI] = new TH1D(("atlasResGen" + powStr  + "_h").c_str(), ";Gen;Counts", nGenBins, genBins);
    
    cmsR4Response_h[pI] = new TH2D(("cmsR4Response" + powStr  + "_h").c_str(), ";Reco;Gen", nRecoBins, recoBins, nGenBins, genBins);
    cmsR4ResReco_h[pI] = new TH1D(("cmsR4ResReco" + powStr  + "_h").c_str(), ";Reco;Counts", nRecoBins, recoBins);
    cmsR4ResGen_h[pI] = new TH1D(("cmsR4ResGen" + powStr  + "_h").c_str(), ";Gen;Counts", nGenBins, genBins);

    cmsPPR4Response_h[pI] = new TH2D(("cmsPPR4Response" + powStr  + "_h").c_str(), ";Reco;Gen", nRecoBins, recoBins, nGenBins, genBins);
    cmsPPR4ResReco_h[pI] = new TH1D(("cmsPPR4ResReco" + powStr  + "_h").c_str(), ";Reco;Counts", nRecoBins, recoBins);
    cmsPPR4ResGen_h[pI] = new TH1D(("cmsPPR4ResGen" + powStr  + "_h").c_str(), ";Gen;Counts", nGenBins, genBins);

    cmsR10Response_h[pI] = new TH2D(("cmsR10Response" + powStr  + "_h").c_str(), ";Reco;Gen", nRecoBins, recoBins, nGenBins, genBins);
    cmsR10ResReco_h[pI] = new TH1D(("cmsR10ResReco" + powStr  + "_h").c_str(), ";Reco;Counts", nRecoBins, recoBins);
    cmsR10ResGen_h[pI] = new TH1D(("cmsR10ResGen" + powStr  + "_h").c_str(), ";Gen;Counts", nGenBins, genBins);

    cmsR10SmearResponse_h[pI] = new TH2D(("cmsR10SmearResponse" + powStr  + "_h").c_str(), ";Reco;Gen", nRecoBins, recoBins, nGenBins, genBins);
    cmsR10SmearResReco_h[pI] = new TH1D(("cmsR10SmearResReco" + powStr  + "_h").c_str(), ";Reco;Counts", nRecoBins, recoBins);
    cmsR10SmearResGen_h[pI] = new TH1D(("cmsR10SmearResGen" + powStr  + "_h").c_str(), ";Gen;Counts", nGenBins, genBins);
    
    genSpectraPower_p[pI] = new TH1D(("genSpectra" + powStr + "_h").c_str(), ";gen;counts", nGenBins, genBins);
    genSpectraATLASTruncPower_p[pI] = new TH1D(("genSpectraATLASTrunc" + powStr + "_h").c_str(), ";gen;counts", nGenBins, genBins);
    recoATLASSpectraPower_p[pI] = new TH1D(("recoATLASSpectra" + powStr + "_h").c_str(), ";recoATLAS;counts", nRecoBins, recoBins);
    genSpectraCMSR4TruncPower_p[pI] = new TH1D(("genSpectraCMSR4Trunc" + powStr + "_h").c_str(), ";gen;counts", nGenBins, genBins);
    recoCMSR4SpectraPower_p[pI] = new TH1D(("recoCMSR4Spectra" + powStr + "_h").c_str(), ";recoCMSR4;counts", nRecoBins, recoBins);

    genSpectraCMSPPR4TruncPower_p[pI] = new TH1D(("genSpectraCMSPPR4Trunc" + powStr + "_h").c_str(), ";gen;counts", nGenBins, genBins);
    recoCMSPPR4SpectraPower_p[pI] = new TH1D(("recoCMSPPR4Spectra" + powStr + "_h").c_str(), ";recoCMSPPR4;counts", nRecoBins, recoBins);

    genSpectraCMSR10TruncPower_p[pI] = new TH1D(("genSpectraCMSR10Trunc" + powStr + "_h").c_str(), ";gen;counts", nGenBins, genBins);
    recoCMSR10SpectraPower_p[pI] = new TH1D(("recoCMSR10Spectra" + powStr + "_h").c_str()
, ";recoCMSR10;counts", nRecoBins, recoBins);

    genSpectraCMSR10SmearTruncPower_p[pI] = new TH1D(("genSpectraCMSR10SmearTrunc" + powStr + "_h").c_str(), ";gen;counts", nGenBins, genBins);
    recoCMSR10SmearSpectraPower_p[pI] = new TH1D(("recoCMSR10SmearSpectra" + powStr + "_h").c_str()
, ";recoCMSR10Smear;counts", nRecoBins, recoBins);


    recoGenSpectraATLASTruncPower_p[pI] = new TH2D(("recoGenSpectraATLASTrunc" + powStr + "_h").c_str(), ";reco;gen", nRecoBins, recoBins, nGenBins, genBins);
    recoGenSpectraCMSR4TruncPower_p[pI] = new TH2D(("recoGenSpectraCMSR4Trunc" + powStr + "_h").c_str(), ";reco;gen", nRecoBins, recoBins, nGenBins, genBins);
    recoGenSpectraCMSPPR4TruncPower_p[pI] = new TH2D(("recoGenSpectraCMSPPR4Trunc" + powStr + "_h").c_str(), ";reco;gen", nRecoBins, recoBins, nGenBins, genBins);
    recoGenSpectraCMSR10TruncPower_p[pI] = new TH2D(("recoGenSpectraCMSR10Trunc" + powStr + "_h").c_str(), ";reco;gen", nRecoBins, recoBins, nGenBins, genBins);

    recoGenSpectraCMSR10SmearTruncPower_p[pI] = new TH2D(("recoGenSpectraCMSR10SmearTrunc" + powStr + "_h").c_str(), ";reco;gen", nRecoBins, recoBins, nGenBins, genBins);
  
    for(Int_t rI = 0; rI < nRandomGen; ++rI){
      for(Int_t tI = 0; tI < nTarget; ++tI){
	genSpectraATLASTruncPower_Rand_p[pI][rI][tI] = new TH1D(("genSpectraATLASTrunc" + powStr + "_Random" + std::to_string(rI) + "_Target" + std::to_string(tI) + "_h").c_str(), ";reco;gen", nGenBins, genBins);
	genSpectraCMSR4TruncPower_Rand_p[pI][rI][tI] = new TH1D(("genSpectraCMSR4Trunc" + powStr + "_Random" + std::to_string(rI) + "_Target" + std::to_string(tI) + "_h").c_str(), ";reco;gen", nGenBins, genBins);
	genSpectraCMSPPR4TruncPower_Rand_p[pI][rI][tI] = new TH1D(("genSpectraCMSPPR4Trunc" + powStr + "_Random" + std::to_string(rI) + "_Target" + std::to_string(tI) + "_h").c_str(), ";reco;gen", nGenBins, genBins);
	genSpectraCMSR10TruncPower_Rand_p[pI][rI][tI] = new TH1D(("genSpectraCMSR10Trunc" + powStr + "_Random" + std::to_string(rI) + "_Target" + std::to_string(tI) + "_h").c_str(), ";reco;gen", nGenBins, genBins);
	genSpectraCMSR10SmearTruncPower_Rand_p[pI][rI][tI] = new TH1D(("genSpectraCMSR10SmearTrunc" + powStr + "_Random" + std::to_string(rI) + "_Target" + std::to_string(tI) + "_h").c_str(), ";reco;gen", nGenBins, genBins);
	
	recoGenSpectraATLASTruncPower_Rand_p[pI][rI][tI] = new TH2D(("recoGenSpectraATLASTrunc" + powStr + "_Random" + std::to_string(rI) + "_Target" + std::to_string(tI) + "_h").c_str(), ";reco;gen", nRecoBins, recoBins, nGenBins, genBins);
	recoGenSpectraCMSR4TruncPower_Rand_p[pI][rI][tI] = new TH2D(("recoGenSpectraCMSR4Trunc" + powStr + "_Random" + std::to_string(rI) + "_Target" + std::to_string(tI) + "_h").c_str(), ";reco;gen", nRecoBins, recoBins, nGenBins, genBins);
	recoGenSpectraCMSPPR4TruncPower_Rand_p[pI][rI][tI] = new TH2D(("recoGenSpectraCMSPPR4Trunc" + powStr + "_Random" + std::to_string(rI) + "_Target" + std::to_string(tI) + "_h").c_str(), ";reco;gen", nRecoBins, recoBins, nGenBins, genBins);
	recoGenSpectraCMSR10TruncPower_Rand_p[pI][rI][tI] = new TH2D(("recoGenSpectraCMSR10Trunc" + powStr + "_Random" + std::to_string(rI) + "_Target" + std::to_string(tI) + "_h").c_str(), ";reco;gen", nRecoBins, recoBins, nGenBins, genBins);
	recoGenSpectraCMSR10SmearTruncPower_Rand_p[pI][rI][tI] = new TH2D(("recoGenSpectraCMSR10SmearTrunc" + powStr + "_Random" + std::to_string(rI) + "_Target" + std::to_string(tI) + "_h").c_str(), ";reco;gen", nRecoBins, recoBins, nGenBins, genBins);

	setSumW2({genSpectraATLASTruncPower_Rand_p[pI][rI][tI], genSpectraCMSR4TruncPower_Rand_p[pI][rI][tI], genSpectraCMSPPR4TruncPower_Rand_p[pI][rI][tI], genSpectraCMSR10TruncPower_Rand_p[pI][rI][tI], genSpectraCMSR10SmearTruncPower_Rand_p[pI][rI][tI],  recoGenSpectraATLASTruncPower_Rand_p[pI][rI][tI], recoGenSpectraCMSR4TruncPower_Rand_p[pI][rI][tI], recoGenSpectraCMSPPR4TruncPower_Rand_p[pI][rI][tI], recoGenSpectraCMSR10TruncPower_Rand_p[pI][rI][tI], recoGenSpectraCMSR10SmearTruncPower_Rand_p[pI][rI][tI]});
      }
    }

    setSumW2({atlasResponse_h[pI], atlasResReco_h[pI], atlasResGen_h[pI], cmsR4Response_h[pI], cmsR4ResReco_h[pI], cmsR4ResGen_h[pI], cmsPPR4Response_h[pI], cmsPPR4ResReco_h[pI], cmsPPR4ResGen_h[pI], cmsR10Response_h[pI], cmsR10ResReco_h[pI], cmsR10ResGen_h[pI], cmsR10SmearResponse_h[pI], cmsR10SmearResReco_h[pI], cmsR10SmearResGen_h[pI], genSpectraPower_p[pI], genSpectraATLASTruncPower_p[pI], recoATLASSpectraPower_p[pI], genSpectraCMSR4TruncPower_p[pI], recoCMSR4SpectraPower_p[pI], genSpectraCMSPPR4TruncPower_p[pI], recoCMSPPR4SpectraPower_p[pI], genSpectraCMSR10TruncPower_p[pI], recoCMSR10SpectraPower_p[pI], genSpectraCMSR10SmearTruncPower_p[pI], recoCMSR10SmearSpectraPower_p[pI], recoGenSpectraATLASTruncPower_p[pI], recoGenSpectraCMSR4TruncPower_p[pI], recoGenSpectraCMSPPR4TruncPower_p[pI], recoGenSpectraCMSR10TruncPower_p[pI], recoGenSpectraCMSR10SmearTruncPower_p[pI]});
  }

  ULong64_t nFill = (ULong64_t)number;
  Double_t fillWeight = 0.01;

  cppWatch response1, response2;

  response1.start();

  for(Int_t pI = 0; pI < nPowers; ++pI){
    for(ULong64_t gI = 0; gI < nFill; ++gI){
      Double_t prob = randGen_p->Uniform(TMath::Power(genBins[0]/genBins[nGenBins], 1), 1);
      Double_t genPt = genBins[0]/TMath::Power(prob, 1./1.);
      Double_t genPtToFill = genPt;
      Double_t weight2 = TMath::Power(genBins[0]/genPt, powers[pI]-2);
      fillWeight = weight2;
      //      fillWeight = 1.0;

      //    if(genPt < genBins[0]) genPtToFill = genBins[0]+1;
      //    if(genPt >= genBins[nGenBins]) genPtToFill = genBins[nGenBins]-1;
      
      if(genPt < genBins[0]) continue;
      if(genPt >= genBins[nGenBins]) continue;
      
      Double_t resATLAS = TMath::Sqrt(atlasPar0*atlasPar0 + atlasPar1*atlasPar1/genPt + atlasPar2*atlasPar2/(genPt*genPt));
      Double_t resCMSR4 = TMath::Sqrt(cmsR4Par0*cmsR4Par0 + cmsR4Par1*cmsR4Par1/genPt + cmsR4Par2*cmsR4Par2/(genPt*genPt));
      Double_t resCMSPPR4 = TMath::Sqrt(cmsPPR4Par0*cmsPPR4Par0 + cmsPPR4Par1*cmsPPR4Par1/genPt + cmsPPR4Par2*cmsPPR4Par2/(genPt*genPt));
      Double_t resCMSR10 = TMath::Sqrt(cmsR10Par0*cmsR10Par0 + cmsR10Par1*cmsR10Par1/genPt + cmsR10Par2*cmsR10Par2/(genPt*genPt));
      Double_t resCMSR10Smear = resCMSR10*1.1*1.1;
      
      Double_t recoATLAS = genPt*randGen_p->Gaus(1., resATLAS);
      Double_t recoCMSR4 = genPt*randGen_p->Gaus(1., resCMSR4);
      Double_t recoCMSPPR4 = genPt*randGen_p->Gaus(1., resCMSPPR4);
      Double_t recoCMSR10 = genPt*randGen_p->Gaus(1., resCMSR10);
      Double_t recoCMSR10Smear = genPt*randGen_p->Gaus(1., resCMSR10Smear);
      
      if(recoATLAS >= recoBins[0] && recoATLAS < recoBins[nRecoBins]){
	atlasResponse_h[pI]->Fill(recoATLAS, genPtToFill, fillWeight);
	atlasResReco_h[pI]->Fill(recoATLAS, fillWeight);
	atlasResGen_h[pI]->Fill(genPtToFill, fillWeight);
      }
      
      if(recoCMSR4 >= recoBins[0] && recoCMSR4 < recoBins[nRecoBins]){
	cmsR4Response_h[pI]->Fill(recoCMSR4, genPtToFill, fillWeight);
	cmsR4ResReco_h[pI]->Fill(recoCMSR4, fillWeight);
	cmsR4ResGen_h[pI]->Fill(genPtToFill, fillWeight);
      }

      if(recoCMSPPR4 >= recoBins[0] && recoCMSPPR4 < recoBins[nRecoBins]){
	cmsPPR4Response_h[pI]->Fill(recoCMSPPR4, genPtToFill, fillWeight);
	cmsPPR4ResReco_h[pI]->Fill(recoCMSPPR4, fillWeight);
	cmsPPR4ResGen_h[pI]->Fill(genPtToFill, fillWeight);
      }

      if(recoCMSR10 >= recoBins[0] && recoCMSR10 < recoBins[nRecoBins]){
	cmsR10Response_h[pI]->Fill(recoCMSR10, genPtToFill, fillWeight);
	cmsR10ResReco_h[pI]->Fill(recoCMSR10, fillWeight);
	cmsR10ResGen_h[pI]->Fill(genPtToFill, fillWeight);
      }

      if(recoCMSR10Smear >= recoBins[0] && recoCMSR10Smear < recoBins[nRecoBins]){
	cmsR10SmearResponse_h[pI]->Fill(recoCMSR10Smear, genPtToFill, fillWeight);
	cmsR10SmearResReco_h[pI]->Fill(recoCMSR10Smear, fillWeight);
	cmsR10SmearResGen_h[pI]->Fill(genPtToFill, fillWeight);
      }
    }
  }
  
  response1.stop();

  /*
  response2.start();
  constructResponse(atlasPar0, atlasPar1, atlasPar2, nFill, atlas2Response_h);
  constructResponse(cmsR4Par0, cmsR4Par1, cmsR4Par2, nFill, cms2Response_h);
  response2.stop();
  */
  std::cout << "TIMING: " << std::endl;
  std::cout << " First construction: " << response1.totalWall() << ", " << response1.totalCPU() << std::endl;
  std::cout << " Second construction: " << response2.totalWall() << ", " << response2.totalCPU() << std::endl;

  for(Int_t pI = 0; pI < nPowers; ++pI){
    for(ULong64_t gI = 0; gI < nFill; ++gI){
      Double_t prob = randGen_p->Uniform(TMath::Power(genBins[0]/genBins[nGenBins], 1), 1);
      Double_t genPt = genBins[0]/TMath::Power(prob, 1./1.);
      Double_t genPtToFill = genPt;
      Double_t weight2 = TMath::Power(genBins[0]/genPt, powers[pI]-2);
      fillWeight = weight2;
      //      fillWeight = 1.0;
      
      //    if(genPt < genBins[0]) genPtToFill = genBins[0]+1;
      //    if(genPt >= genBins[nGenBins]) genPtToFill = genBins[nGenBins]-1;
      
      if(genPt < genBins[0]) continue;
      if(genPt >= genBins[nGenBins]) continue;
      
      genSpectraPower_p[pI]->Fill(genPtToFill, fillWeight);
      
      Double_t resATLAS = TMath::Sqrt(atlasPar0*atlasPar0 + atlasPar1*atlasPar1/genPt + atlasPar2*atlasPar2/(genPt*genPt));
      Double_t resCMSR4 = TMath::Sqrt(cmsR4Par0*cmsR4Par0 + cmsR4Par1*cmsR4Par1/genPt + cmsR4Par2*cmsR4Par2/(genPt*genPt));
      Double_t resCMSPPR4 = TMath::Sqrt(cmsPPR4Par0*cmsPPR4Par0 + cmsPPR4Par1*cmsPPR4Par1/genPt + cmsPPR4Par2*cmsPPR4Par2/(genPt*genPt));
      Double_t resCMSR10 = TMath::Sqrt(cmsR10Par0*cmsR10Par0 + cmsR10Par1*cmsR10Par1/genPt + cmsR10Par2*cmsR10Par2/(genPt*genPt));
      Double_t resCMSR10Smear = resCMSR10*1.1*1.1;
      
      Double_t recoATLAS = genPt*randGen_p->Gaus(1., resATLAS);
      Double_t recoCMSR4 = genPt*randGen_p->Gaus(1., resCMSR4);
      Double_t recoCMSPPR4 = genPt*randGen_p->Gaus(1., resCMSPPR4);
      Double_t recoCMSR10 = genPt*randGen_p->Gaus(1., resCMSR10);
      Double_t recoCMSR10Smear = genPt*randGen_p->Gaus(1., resCMSR10Smear);

      
      if(recoATLAS >= recoBins[0] && recoATLAS < recoBins[nRecoBins]){
	recoATLASSpectraPower_p[pI]->Fill(recoATLAS, fillWeight);
	genSpectraATLASTruncPower_p[pI]->Fill(genPtToFill, fillWeight);
	recoGenSpectraATLASTruncPower_p[pI]->Fill(recoATLAS, genPtToFill, fillWeight);
      }
      
      if(recoCMSR4 >= recoBins[0] && recoCMSR4 < recoBins[nRecoBins]){
	recoCMSR4SpectraPower_p[pI]->Fill(recoCMSR4, fillWeight);
	genSpectraCMSR4TruncPower_p[pI]->Fill(genPtToFill, fillWeight);
	recoGenSpectraCMSR4TruncPower_p[pI]->Fill(recoCMSR4, genPtToFill, fillWeight);
      }      

      if(recoCMSPPR4 >= recoBins[0] && recoCMSPPR4 < recoBins[nRecoBins]){
	recoCMSPPR4SpectraPower_p[pI]->Fill(recoCMSPPR4, fillWeight);
	genSpectraCMSPPR4TruncPower_p[pI]->Fill(genPtToFill, fillWeight);
	recoGenSpectraCMSPPR4TruncPower_p[pI]->Fill(recoCMSPPR4, genPtToFill, fillWeight);
      }      

      if(recoCMSR10 >= recoBins[0] && recoCMSR10 < recoBins[nRecoBins]){
	recoCMSR10SpectraPower_p[pI]->Fill(recoCMSR10, fillWeight);
	genSpectraCMSR10TruncPower_p[pI]->Fill(genPtToFill, fillWeight);
	recoGenSpectraCMSR10TruncPower_p[pI]->Fill(recoCMSR10, genPtToFill, fillWeight);
      }      

      if(recoCMSR10Smear >= recoBins[0] && recoCMSR10Smear < recoBins[nRecoBins]){
	recoCMSR10SmearSpectraPower_p[pI]->Fill(recoCMSR10Smear, fillWeight);
	genSpectraCMSR10SmearTruncPower_p[pI]->Fill(genPtToFill, fillWeight);
	recoGenSpectraCMSR10SmearTruncPower_p[pI]->Fill(recoCMSR10Smear, genPtToFill, fillWeight);
      }      
    }
  }


  //  const Double_t errTarget = 5382.0000/20.;
 
  for(Int_t pI = 0; pI < nPowers; ++pI){
    std::cout << "Building distribution " << pI << "/" << nPowers << std::endl;
    for(Int_t rI = 0; rI < nRandomGen; ++rI){
      std::cout << " " << rI << "/" << nRandomGen << std::endl;

      for(Int_t tI = 0; tI < nTarget; ++tI){
	
	Double_t val = 0.0;
	Double_t err = 0.0;
	Double_t rat = 1.0;
	
	while(rat > 1./TMath::Sqrt(errTarget*TMath::Power(10, tI))){
	  Double_t prob = randGen_p->Uniform(TMath::Power(genBins[0]/genBins[nGenBins], 3-tI), 1);
	  Double_t genPt = genBins[0]/TMath::Power(prob, 1./(3.-tI));	

	  Double_t resATLAS = TMath::Sqrt(atlasPar0*atlasPar0 + atlasPar1*atlasPar1/genPt + atlasPar2*atlasPar2/(genPt*genPt));
	  Double_t recoATLAS = genPt*randGen_p->Gaus(1., resATLAS);
	  
	  if(recoATLAS >= recoBins[0] && recoATLAS < recoBins[nRecoBins]){
	    Double_t weight = TMath::Power(genBins[0]/genPt, powers[pI]-(4-tI));
	    genSpectraATLASTruncPower_Rand_p[pI][rI][tI]->Fill(genPt, weight);
	    recoGenSpectraATLASTruncPower_Rand_p[pI][rI][tI]->Fill(recoATLAS, genPt, weight);
	    
	    if(recoATLAS > targetPoint){
	      val += weight;
	      err = TMath::Sqrt(err*err + weight*weight);
	    }
	  }
	  
	  if(val <= TMath::Power(10, -20)) rat = 1.0;
	  else rat = err/val;
	}
	
	val = 0.0;
	err = 0.0;
	rat = 1.0;
	
	while(rat > 1./TMath::Sqrt(errTarget*TMath::Power(10, tI))){
	  Double_t prob = randGen_p->Uniform(TMath::Power(genBins[0]/genBins[nGenBins], 3-tI), 1);
	  Double_t genPt = genBins[0]/TMath::Power(prob, 1./(3.-tI));	
	  
	  Double_t resCMSR4 = TMath::Sqrt(cmsR4Par0*cmsR4Par0 + cmsR4Par1*cmsR4Par1/genPt + cmsR4Par2*cmsR4Par2/(genPt*genPt));
	  Double_t recoCMSR4 = genPt*randGen_p->Gaus(1., resCMSR4);
	  
	  if(recoCMSR4 >= recoBins[0] && recoCMSR4 < recoBins[nRecoBins]){
	    Double_t weight = TMath::Power(genBins[0]/genPt, powers[pI]-(4-tI));
	    genSpectraCMSR4TruncPower_Rand_p[pI][rI][tI]->Fill(genPt, weight);
	    recoGenSpectraCMSR4TruncPower_Rand_p[pI][rI][tI]->Fill(recoCMSR4, genPt, weight);
	    
	    if(recoCMSR4 > targetPoint){
	      val += weight;
	      err = TMath::Sqrt(err*err + weight*weight);
	    }
	  }
	  
	  if(val <= TMath::Power(10, -20)) rat = 1.0;
	  else rat = err/val;
	}
	
	val = 0.0;
	err = 0.0;
	rat = 1.0;
	
	while(rat > 1./TMath::Sqrt(errTarget*TMath::Power(10, tI))){
	  Double_t prob = randGen_p->Uniform(TMath::Power(genBins[0]/genBins[nGenBins], 3-tI), 1);
	  Double_t genPt = genBins[0]/TMath::Power(prob, 1./(3.-tI));	
	  
	  Double_t resCMSPPR4 = TMath::Sqrt(cmsPPR4Par0*cmsPPR4Par0 + cmsPPR4Par1*cmsPPR4Par1/genPt + cmsPPR4Par2*cmsPPR4Par2/(genPt*genPt));
	  Double_t recoCMSPPR4 = genPt*randGen_p->Gaus(1., resCMSPPR4);
	  
	  if(recoCMSPPR4 >= recoBins[0] && recoCMSPPR4 < recoBins[nRecoBins]){
	    Double_t weight = TMath::Power(genBins[0]/genPt, powers[pI]-(4-tI));
	    genSpectraCMSPPR4TruncPower_Rand_p[pI][rI][tI]->Fill(genPt, weight);
	    recoGenSpectraCMSPPR4TruncPower_Rand_p[pI][rI][tI]->Fill(recoCMSPPR4, genPt, weight);
	    
	    if(recoCMSPPR4 > targetPoint){
	      val += weight;
	      err = TMath::Sqrt(err*err + weight*weight);
	    }
	  }
	  
	  if(val <= TMath::Power(10, -20)) rat = 1.0;
	  else rat = err/val;
	}
	
	val = 0.0;
	err = 0.0;
	rat = 1.0;
	
	while(rat > 1./TMath::Sqrt(errTarget*TMath::Power(10, tI))){
	  Double_t prob = randGen_p->Uniform(TMath::Power(genBins[0]/genBins[nGenBins], 3-tI), 1);
	  Double_t genPt = genBins[0]/TMath::Power(prob, 1./(3.-tI));	
	  
	  Double_t resCMSR10 = TMath::Sqrt(cmsR10Par0*cmsR10Par0 + cmsR10Par1*cmsR10Par1/genPt + cmsR10Par2*cmsR10Par2/(genPt*genPt));
	  Double_t recoCMSR10 = genPt*randGen_p->Gaus(1., resCMSR10);
	  
	  if(recoCMSR10 >= recoBins[0] && recoCMSR10 < recoBins[nRecoBins]){
	    Double_t weight = TMath::Power(genBins[0]/genPt, powers[pI]-(4-tI));
	    genSpectraCMSR10TruncPower_Rand_p[pI][rI][tI]->Fill(genPt, weight);
	    recoGenSpectraCMSR10TruncPower_Rand_p[pI][rI][tI]->Fill(recoCMSR10, genPt, weight);
	    
	    if(recoCMSR10 > targetPoint){
	      val += weight;
	      err = TMath::Sqrt(err*err + weight*weight);
	    }
	  }
	  
	  if(val <= TMath::Power(10, -20)) rat = 1.0;
	  else rat = err/val;
	}
	

	val = 0.0;
	err = 0.0;
	rat = 1.0;
	
       	while(rat > 1./TMath::Sqrt(errTarget*TMath::Power(10, tI))){
	  Double_t prob = randGen_p->Uniform(TMath::Power(genBins[0]/genBins[nGenBins], 3-tI), 1);
	  Double_t genPt = genBins[0]/TMath::Power(prob, 1./(3.-tI));	
	  
	  Double_t resCMSR10Smear = 1.1*1.1*TMath::Sqrt(cmsR10Par0*cmsR10Par0 + cmsR10Par1*cmsR10Par1/genPt + cmsR10Par2*cmsR10Par2/(genPt*genPt));
	  Double_t recoCMSR10Smear = genPt*randGen_p->Gaus(1., resCMSR10Smear);

	  if(recoCMSR10Smear >= recoBins[0] && recoCMSR10Smear < recoBins[nRecoBins]){
	    Double_t weight = TMath::Power(genBins[0]/genPt, powers[pI]-(4-tI));
	    genSpectraCMSR10SmearTruncPower_Rand_p[pI][rI][tI]->Fill(genPt, weight);
	    recoGenSpectraCMSR10SmearTruncPower_Rand_p[pI][rI][tI]->Fill(recoCMSR10Smear, genPt, weight);
	    
	    if(recoCMSR10Smear > targetPoint){
	      val += weight;
	      err = TMath::Sqrt(err*err + weight*weight);
	    }
	  }
	  
	  if(val <= TMath::Power(10, -20)) rat = 1.0;
	  else rat = err/val;
	}
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
    cmsR4Response_h[pI]->Write("", TObject::kOverwrite);
    cmsR4ResReco_h[pI]->Write("", TObject::kOverwrite);
    cmsR4ResGen_h[pI]->Write("", TObject::kOverwrite);

    cmsPPR4Response_h[pI]->Write("", TObject::kOverwrite);
    cmsPPR4ResReco_h[pI]->Write("", TObject::kOverwrite);
    cmsPPR4ResGen_h[pI]->Write("", TObject::kOverwrite);

    cmsR10Response_h[pI]->Write("", TObject::kOverwrite);
    cmsR10ResReco_h[pI]->Write("", TObject::kOverwrite);
    cmsR10ResGen_h[pI]->Write("", TObject::kOverwrite);

    cmsR10SmearResponse_h[pI]->Write("", TObject::kOverwrite);
    cmsR10SmearResReco_h[pI]->Write("", TObject::kOverwrite);
    cmsR10SmearResGen_h[pI]->Write("", TObject::kOverwrite);
    
    genSpectraPower_p[pI]->Write("", TObject::kOverwrite);

    genSpectraATLASTruncPower_p[pI]->Write("", TObject::kOverwrite);
    recoATLASSpectraPower_p[pI]->Write("", TObject::kOverwrite);
    recoGenSpectraATLASTruncPower_p[pI]->Write("", TObject::kOverwrite);

    genSpectraCMSR4TruncPower_p[pI]->Write("", TObject::kOverwrite);
    recoCMSR4SpectraPower_p[pI]->Write("", TObject::kOverwrite);
    recoGenSpectraCMSR4TruncPower_p[pI]->Write("", TObject::kOverwrite);

    genSpectraCMSPPR4TruncPower_p[pI]->Write("", TObject::kOverwrite);
    recoCMSPPR4SpectraPower_p[pI]->Write("", TObject::kOverwrite);
    recoGenSpectraCMSPPR4TruncPower_p[pI]->Write("", TObject::kOverwrite);

    genSpectraCMSR10TruncPower_p[pI]->Write("", TObject::kOverwrite);
    recoCMSR10SpectraPower_p[pI]->Write("", TObject::kOverwrite);
    recoGenSpectraCMSR10TruncPower_p[pI]->Write("", TObject::kOverwrite);

    genSpectraCMSR10SmearTruncPower_p[pI]->Write("", TObject::kOverwrite);
    recoCMSR10SmearSpectraPower_p[pI]->Write("", TObject::kOverwrite);
    recoGenSpectraCMSR10SmearTruncPower_p[pI]->Write("", TObject::kOverwrite);

    for(Int_t rI = 0; rI < nRandomGen; ++rI){
      for(Int_t tI = 0; tI < nTarget; ++tI){
	genSpectraATLASTruncPower_Rand_p[pI][rI][tI]->Write("", TObject::kOverwrite);
	genSpectraCMSR4TruncPower_Rand_p[pI][rI][tI]->Write("", TObject::kOverwrite);
	genSpectraCMSPPR4TruncPower_Rand_p[pI][rI][tI]->Write("", TObject::kOverwrite);
	genSpectraCMSR10TruncPower_Rand_p[pI][rI][tI]->Write("", TObject::kOverwrite);
	genSpectraCMSR10SmearTruncPower_Rand_p[pI][rI][tI]->Write("", TObject::kOverwrite);
	
	recoGenSpectraATLASTruncPower_Rand_p[pI][rI][tI]->Write("", TObject::kOverwrite);
	recoGenSpectraCMSR4TruncPower_Rand_p[pI][rI][tI]->Write("", TObject::kOverwrite);
	recoGenSpectraCMSPPR4TruncPower_Rand_p[pI][rI][tI]->Write("", TObject::kOverwrite);
	recoGenSpectraCMSR10TruncPower_Rand_p[pI][rI][tI]->Write("", TObject::kOverwrite);
	recoGenSpectraCMSR10SmearTruncPower_Rand_p[pI][rI][tI]->Write("", TObject::kOverwrite);
      }
    }
  }
  

  for(Int_t pI = 0; pI < nPowers; ++pI){
    delete atlasResponse_h[pI];
    delete atlasResReco_h[pI];
    delete atlasResGen_h[pI];
    delete cmsR4Response_h[pI];
    delete cmsR4ResReco_h[pI];
    delete cmsR4ResGen_h[pI];

    delete cmsPPR4Response_h[pI];
    delete cmsPPR4ResReco_h[pI];
    delete cmsPPR4ResGen_h[pI];

    delete cmsR10Response_h[pI];
    delete cmsR10ResReco_h[pI];
    delete cmsR10ResGen_h[pI];

    delete cmsR10SmearResponse_h[pI];
    delete cmsR10SmearResReco_h[pI];
    delete cmsR10SmearResGen_h[pI];
    
    delete genSpectraPower_p[pI];

    delete genSpectraATLASTruncPower_p[pI];
    delete recoATLASSpectraPower_p[pI];
    delete recoGenSpectraATLASTruncPower_p[pI];

    delete genSpectraCMSR4TruncPower_p[pI];
    delete recoCMSR4SpectraPower_p[pI];
    delete recoGenSpectraCMSR4TruncPower_p[pI];

    delete genSpectraCMSPPR4TruncPower_p[pI];
    delete recoCMSPPR4SpectraPower_p[pI];
    delete recoGenSpectraCMSPPR4TruncPower_p[pI];

    delete genSpectraCMSR10TruncPower_p[pI];
    delete recoCMSR10SpectraPower_p[pI];
    delete recoGenSpectraCMSR10TruncPower_p[pI];

    delete genSpectraCMSR10SmearTruncPower_p[pI];
    delete recoCMSR10SmearSpectraPower_p[pI];
    delete recoGenSpectraCMSR10SmearTruncPower_p[pI];

    for(Int_t rI = 0; rI < nRandomGen; ++rI){
      for(Int_t tI = 0; tI < nTarget; ++tI){
	delete genSpectraATLASTruncPower_Rand_p[pI][rI][tI];
	delete genSpectraCMSR4TruncPower_Rand_p[pI][rI][tI];
	delete genSpectraCMSPPR4TruncPower_Rand_p[pI][rI][tI];
	delete genSpectraCMSR10TruncPower_Rand_p[pI][rI][tI];
	delete genSpectraCMSR10SmearTruncPower_Rand_p[pI][rI][tI];
	
	delete recoGenSpectraATLASTruncPower_Rand_p[pI][rI][tI];
	delete recoGenSpectraCMSR4TruncPower_Rand_p[pI][rI][tI];
	delete recoGenSpectraCMSPPR4TruncPower_Rand_p[pI][rI][tI];
	delete recoGenSpectraCMSR10TruncPower_Rand_p[pI][rI][tI];
	delete recoGenSpectraCMSR10SmearTruncPower_Rand_p[pI][rI][tI];
      }
    }
  }
  
  outFile_p->Close();
  delete outFile_p;

  delete randGen_p;
  
  total.stop();

  std::cout << "Total timing: " << std::endl;
  std::cout << " CPU: " << total.totalCPU() << std::endl;
  std::cout << " Wall: " << total.totalWall() << std::endl;

  return 0;
}


int main(int argc, char* argv[])
{
  if(argc != 3){
    std::cout << "Usage: ./bin/toyATLASCMS.exe <seed> <number>" << std::endl;
    return 1;
  }
  
  int retVal = 0;
  retVal += toyATLASCMS(std::stoi(argv[1]), std::stol(argv[2]));
  return retVal;
}
