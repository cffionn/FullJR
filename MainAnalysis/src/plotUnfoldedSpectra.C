//cpp dependencies
#include <iostream>
#include <string>
#include <vector>

//ROOT dependencies
#include "TFile.h"
#include "TH1D.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TDatime.h"
#include "TMath.h"
#include "TCollection.h"
#include "TKey.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TBox.h"
#include "TError.h"

//Local dependencies
#include "MainAnalysis/include/cutPropagator.h"
#include "MainAnalysis/include/doLocalDebug.h"
#include "MainAnalysis/include/texSlideCreator.h"

//Non-local dependencies
#include "Utility/include/checkMakeDir.h"
#include "Utility/include/cppWatch.h"
#include "Utility/include/doGlobalDebug.h"
#include "Utility/include/getLogBins.h"
#include "Utility/include/histDefUtility.h"
#include "Utility/include/kirchnerPalette.h"
#include "Utility/include/vanGoghPalette.h"
#include "Utility/include/lumiAndTAAUtil.h"
#include "Utility/include/returnRootFileContentsList.h"


std::vector<double> getSyst(TH1D* nominal_p, std::vector<TH1D*> syst_p, std::vector<std::string> systStr, Double_t minXVal, Double_t maxXVal, std::vector<std::string>* plotNames, bool doSmoothing)
{
  //Setting some style and color for TH1
  vanGoghPalette vg;
  const Int_t nColor = 5;
  const Int_t colors[nColor] = {vg.getColor(0), vg.getColor(1), vg.getColor(2), vg.getColor(3), vg.getColor(4)};

  const Int_t nStyle = 4;
  Int_t styles[nStyle] = {20, 21, 47, 34};

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  double smallDelta = 0.000001;
  std::vector<double> systVect;

  //Building parallel systematic histograms
  Int_t tempNBins = 0;
  for(Int_t bIX = 0; bIX < nominal_p->GetNbinsX(); ++bIX){
    Double_t binCenter = nominal_p->GetBinCenter(bIX+1);
    if(binCenter > maxXVal || binCenter < minXVal) continue;
    
    tempNBins++;
  }

  const Int_t nBins = tempNBins;
  Double_t bins[nBins+1];

  for(Int_t bIX = 0; bIX < nominal_p->GetNbinsX(); ++bIX){
    Double_t binCenter = nominal_p->GetBinCenter(bIX+1);
    if(binCenter < minXVal) continue;
    if(binCenter > maxXVal){
      bins[bIX] = nominal_p->GetBinLowEdge(bIX+1);
      break;
    }

    bins[bIX] = nominal_p->GetBinLowEdge(bIX+1);
  }

  const Int_t nSystHist = syst_p.size();
  TH1D* systHist_p[nSystHist];

  for(Int_t sI = 0; sI < nSystHist; ++sI){
    systHist_p[sI] = new TH1D(("systHist_" + std::to_string(sI)).c_str(), "", nBins, bins);

    Int_t binFillPos = 0;
    for(Int_t bIX = 0; bIX < syst_p.at(sI)->GetNbinsX(); ++bIX){
      Double_t binCenter = nominal_p->GetBinCenter(bIX+1);
      if(binCenter > maxXVal || binCenter < minXVal) continue;

      Double_t binVal = TMath::Abs(syst_p.at(sI)->GetBinContent(bIX+1) - nominal_p->GetBinContent(bIX+1))/(nominal_p->GetBinContent(bIX+1));
      Double_t binErr = nominal_p->GetBinContent(bIX+1)*nominal_p->GetBinContent(bIX+1)/(nominal_p->GetBinError(bIX+1)*nominal_p->GetBinError(bIX+1));
      binErr += syst_p.at(sI)->GetBinContent(bIX+1)*syst_p.at(sI)->GetBinContent(bIX+1)/(syst_p.at(sI)->GetBinError(bIX+1)*syst_p.at(sI)->GetBinError(bIX+1));
      binErr = binVal*TMath::Sqrt(binErr);

      systHist_p[sI]->SetBinContent(binFillPos+1, binVal);
      systHist_p[sI]->SetBinError(binFillPos+1, binErr);
      binFillPos++;
    }

    if(doSmoothing) systHist_p[sI]->Smooth();
    
    binFillPos = 0;

    for(Int_t bIX = 0; bIX < syst_p.at(sI)->GetNbinsX(); ++bIX){
      Double_t binCenter = nominal_p->GetBinCenter(bIX+1);
      if(binCenter > maxXVal || binCenter < minXVal) continue;

      Double_t binVal = systHist_p[sI]->GetBinContent(binFillPos+1)*nominal_p->GetBinContent(bIX+1);
      Double_t binErr = systHist_p[sI]->GetBinError(binFillPos+1)*nominal_p->GetBinContent(bIX+1);

      systHist_p[sI]->SetBinContent(binFillPos+1, binVal);
      systHist_p[sI]->SetBinError(binFillPos+1, binErr);
      binFillPos++;
    }

    delete systHist_p[sI];
  }


  for(Int_t bIX = 0; bIX < nominal_p->GetNbinsX(); ++bIX){
    Double_t binCenter = nominal_p->GetBinCenter(bIX+1);
    if(binCenter > maxXVal || binCenter < minXVal) continue;

    Double_t tempVal = 0.0;
    for(unsigned int sI = 0; sI < syst_p.size(); ++sI){
      tempVal = TMath::Sqrt(tempVal*tempVal + syst_p.at(sI)->GetBinContent(bIX+1)*syst_p.at(sI)->GetBinContent(bIX+1));

      if(binCenter > maxXVal || binCenter < minXVal){
	syst_p.at(sI)->SetBinContent(bIX+1, 0.0);
	syst_p.at(sI)->SetBinError(bIX+1, 0.0);
      }
      else{
	Double_t binVal = syst_p.at(sI)->GetBinContent(bIX+1)/nominal_p->GetBinContent(bIX+1);
	Double_t tempRelErr = syst_p.at(sI)->GetBinError(bIX+1)/nominal_p->GetBinContent(bIX+1);

	syst_p.at(sI)->SetBinContent(bIX+1, binVal);
	syst_p.at(sI)->SetBinError(bIX+1, tempRelErr*binVal);
	if(syst_p.at(sI)->GetBinContent(bIX+1) <= 0){
	  syst_p.at(sI)->SetBinContent(bIX+1, smallDelta);
	  syst_p.at(sI)->SetBinError(bIX+1, smallDelta);
	}
      }
    }

    systVect.push_back(tempVal);
  }

  double maxVal = -1;

  for(unsigned int sI = 0; sI < syst_p.size(); ++sI){
    if(syst_p.at(sI)->GetMaximum() > maxVal) maxVal = syst_p.at(sI)->GetMaximum();
  }

  const Int_t nMax = 2;
  const Double_t max[nMax] = {maxVal, 0.25};
  
  for(Int_t nI = 0; nI < nMax; ++nI){
    TLegend* leg_p = new TLegend(0.2, 0.52, 0.4, 0.99);
    leg_p->SetBorderSize(0);
    leg_p->SetFillColor(0);
    leg_p->SetFillStyle(0);
    leg_p->SetTextFont(43);
    leg_p->SetTextSize(16);
    
    
    TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);
    canv_p->SetTopMargin(0.01);
    canv_p->SetRightMargin(0.01);
    canv_p->SetLeftMargin(0.12);
    canv_p->SetBottomMargin(0.12);
    
    TH1D* tempHist_p = new TH1D("tempHist_p", ";Jet p_{T} (GeV);Fractional Systetmatics", 10, minXVal, maxXVal);
    tempHist_p->GetXaxis()->SetTitleOffset(1.6);
    
    canv_p->cd();
    centerTitles(tempHist_p);
    tempHist_p->SetMaximum(1.2*max[nI]);
    tempHist_p->DrawCopy("HIST E1 P");
    
    gStyle->SetOptStat(0);
    gPad->SetLogx();
    
    for(unsigned int sI = 0; sI < syst_p.size(); ++sI){
      syst_p.at(sI)->SetMarkerColor(colors[sI%nColor]);
      syst_p.at(sI)->SetMarkerStyle(styles[sI%nStyle]);
      syst_p.at(sI)->SetLineColor(colors[sI%nColor]);
      syst_p.at(sI)->SetMarkerSize(0.8);

      syst_p.at(sI)->DrawCopy("HIST E1 P SAME");
      syst_p.at(sI)->DrawCopy("HIST E1 SAME");
      
      leg_p->AddEntry(syst_p.at(sI), systStr.at(sI+1).c_str(), "P L");
    }
    
    leg_p->Draw("SAME");
    
    gPad->RedrawAxis();
    gPad->SetTicks(1,2);

    std::string maxValStr = "NormalMax";
    if(nI != 0) maxValStr = "ZOOM";
    
    std::string smoothStr = "Unsmoothed";
    if(doSmoothing) smoothStr = "Smoothed";

    std::string nominalName = nominal_p->GetName();
    const std::string dirName = "pdfDir/" + dateStr;
    checkMakeDir(dirName);
    const std::string saveName = "systErr_" + nominalName + "_" + maxValStr + "_" + smoothStr + "_" + dateStr +  ".pdf";

    quietSaveAs(canv_p, (dirName + "/" + saveName));
    plotNames->push_back(saveName);

    delete tempHist_p;
    delete canv_p;
    delete leg_p;
  }

  return systVect;
}


void drawSyst(TCanvas* canv_p, TH1D* nominal_p, std::vector<double> syst_, Double_t minXVal, Double_t maxXVal)
{
  canv_p->cd();

  TBox* tempBox_p = new TBox();
  tempBox_p->SetFillColorAlpha(nominal_p->GetMarkerColor(), .25);
  for(Int_t bIX = 0; bIX < nominal_p->GetNbinsX(); ++bIX){
    Double_t binCenter = nominal_p->GetBinCenter(bIX+1);
    Double_t binLowEdge = nominal_p->GetBinLowEdge(bIX+1);
    Double_t binHiEdge = nominal_p->GetBinLowEdge(bIX+2);
    Double_t binContent = nominal_p->GetBinContent(bIX+1);
    
    if(binCenter < minXVal) continue;
    if(binCenter > maxXVal) continue;
    
    tempBox_p->DrawBox(binLowEdge, binContent - syst_.at(bIX), binHiEdge, binContent + syst_.at(bIX));
  }

  return;
}


int plotUnfoldedSpectra(const std::string inFileNamePP, const std::string inFileNamePbPb)
{
  std::vector<std::string> slideTitles;
  std::vector<std::vector<std::string > > pdfPerSlide;

  double total = 0;
  cppWatch toPPFile;
  cppWatch ppFile;
  cppWatch pbpbFile;
  cppWatch comparisonOfFiles;

  toPPFile.start();

  std::cout << "Initializing files..." << std::endl;

  if(!checkFile(inFileNamePP) || !checkFile(inFileNamePbPb)){
    if(!checkFile(inFileNamePP)) std::cout << "inFileNamePP, \'" << inFileNamePP << "\', is invalid. return 1" << std::endl;
    if(!checkFile(inFileNamePbPb)) std::cout << "inFileNamePbPb, \'" << inFileNamePbPb << "\', is invalid. return 1" << std::endl;
    return 1;
  }

  kirchnerPalette kPalette;

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  const std::string dirStr = "pdfDir/" + dateStr;
  checkMakeDir("pdfDir");
  checkMakeDir(dirStr);

  toPPFile.stop();
  total += toPPFile.total();
  std::cout << "To PP File time: " << toPPFile.total() << "/" << total << std::endl;
  ppFile.start();

  TFile* inFilePP_p = new TFile(inFileNamePP.c_str(), "READ");
  cutPropagator cutPropPP;
  cutPropPP.Clean();
  cutPropPP.GetAllVarFromFile(inFilePP_p);

  std::vector<std::string> histTagPP = cutPropPP.GetHistTag();
  std::vector<int> histBestBayesPP = cutPropPP.GetHistBestBayes();

  std::map<std::string, int> histTagMapPP;
  for(unsigned int i = 0; i < histTagPP.size(); ++i){
    histTagMapPP[histTagPP.at(i)] = histBestBayesPP.at(i);
  }


  std::vector<std::string> jetPPList = returnRootFileContentsList(inFilePP_p, "TDirectoryFile", "JetAnalyzer", 1);

  std::cout << "JetPPList: " << std::endl;
  for(unsigned int jI = 0; jI < jetPPList.size(); ++jI){
    std::cout << " " << jI << "/" << jetPPList.size() << ": " << jetPPList.at(jI) << std::endl;
  }

  ppFile.stop();
  total += ppFile.total();
  std::cout << "PP File time: " << ppFile.total() << "/" << total << std::endl;
  pbpbFile.start();

  TFile* inFilePbPb_p = new TFile(inFileNamePbPb.c_str(), "READ");
  cutPropagator cutPropPbPb;
  cutPropPbPb.Clean();
  cutPropPbPb.GetAllVarFromFile(inFilePbPb_p);
  std::vector<std::string> jetPbPbList = returnRootFileContentsList(inFilePbPb_p, "TDirectoryFile", "JetAnalyzer", 1);

  std::vector<std::string> histTagPbPb = cutPropPbPb.GetHistTag();
  std::vector<int> histBestBayesPbPb = cutPropPbPb.GetHistBestBayes();
  std::map<std::string, int> histTagMapPbPb;
  for(unsigned int i = 0; i < histTagPbPb.size(); ++i){
    histTagMapPbPb[histTagPbPb.at(i)] = histBestBayesPbPb.at(i);
  }

  std::cout << "JetPbPbList: " << std::endl;
  for(unsigned int jI = 0; jI < jetPbPbList.size(); ++jI){
    std::cout << " " << jI << "/" << jetPbPbList.size() << ": " << jetPbPbList.at(jI) << std::endl;
  }

  pbpbFile.stop();
  total += pbpbFile.total();
  std::cout << "PBPB File time: " << pbpbFile.total() << "/" << total << std::endl;
  comparisonOfFiles.start();

  if(!cutPropPP.GetIsPP() || cutPropPbPb.GetIsPP() || !cutPropPP.CheckPropagatorsMatch(cutPropPbPb, true, false)){
    if(!cutPropPP.GetIsPP()) std::cout << "inFileNamePP \'" << inFileNamePP << "\' is not pp, return 1" << std::endl;
    if(cutPropPbPb.GetIsPP()) std::cout << "inFileNamePbPb \'" << inFileNamePbPb << "\' is not PbPb, return 1" << std::endl;
    if(!cutPropPP.CheckPropagatorsMatch(cutPropPbPb, true, false)) std::cout << "Cut propagators of pp and pbpb do not match. inspect. return 1" << std::endl;
  
    inFilePP_p->Close();
    delete inFilePP_p;
    
    inFilePbPb_p->Close();
    delete inFilePbPb_p;
    
    return 1;
  }

  comparisonOfFiles.stop();
  total += comparisonOfFiles.total();
  std::cout << "Comp File time: " << comparisonOfFiles.total() << "/" << total << std::endl;

  Int_t minForForLoop = 10000000;
  if(doLocalDebug || doGlobalDebug) minForForLoop = 1;

  std::cout << "Extracting cuts..." << std::endl;

  const Int_t nAdditionalSyst = 4;
  const std::string additionalSyst[nAdditionalSyst] = {"LumiUp", "LumiDown", "TAAUp", "TAADown"};

  const Int_t nJtPP = jetPPList.size();
  const Int_t nJtPbPb = jetPbPbList.size();

  const Int_t nSystOrig = cutPropPbPb.GetNSyst();
  //  const Int_t nSyst = TMath::Min(minForForLoop, nSystOrig + nAdditionalSyst);
  const Int_t nSyst = nSystOrig + nAdditionalSyst;
  std::vector<std::string> systStr = cutPropPbPb.GetSystStr();
  for(Int_t sI = 0; sI < nAdditionalSyst; ++sI){
    systStr.push_back(additionalSyst[sI]);
  }

  const Int_t nBayesCap = 10000;
  const Int_t nBayes = TMath::Min(nBayesCap, cutPropPbPb.GetNBayes());
  std::vector<int> bayesVal = cutPropPbPb.GetBayesVal();

  const Int_t nResponseMod = TMath::Min(minForForLoop, cutPropPbPb.GetNResponseMod());
  //  const Int_t nResponseMod = cutPropPbPb.GetNResponseMod();
  std::vector<double> responseMod = cutPropPbPb.GetResponseMod();

  const Int_t nCentBins = cutPropPbPb.GetNCentBins();
  std::vector<Int_t> centBinsLow = cutPropPbPb.GetCentBinsLow();
  std::vector<Int_t> centBinsHi = cutPropPbPb.GetCentBinsHi();
  std::vector<Double_t> centBinsScalingFact;
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    centBinsScalingFact.push_back(TMath::Power(10., cI+1.));
  }

  //  const Int_t nJtPtBins = cutPropPbPb.GetNJtPtBins();
  //  std::vector<Double_t> jtPtBinsTemp = cutPropPbPb.GetJtPtBins();

  //  const Int_t nJtAbsEtaBins = TMath::Min(minForForLoop, cutPropPbPb.GetNJtAbsEtaBins());
  const Int_t nJtAbsEtaBins = cutPropPbPb.GetNJtAbsEtaBins();
  std::vector<Double_t> jtAbsEtaBinsLow = cutPropPbPb.GetJtAbsEtaBinsLow();
  std::vector<Double_t> jtAbsEtaBinsHi = cutPropPbPb.GetJtAbsEtaBinsHi();

  //  const Int_t nID = cutPropPbPb.GetNID();
  const Int_t nID = TMath::Min(minForForLoop, cutPropPbPb.GetNID());
  std::vector<std::string> idStr = cutPropPbPb.GetIdStr();
  std::vector<double> jtPfCHMFCutLow = cutPropPbPb.GetJtPfCHMFCutLow();
  std::vector<double> jtPfCHMFCutHi = cutPropPbPb.GetJtPfCHMFCutHi();
  std::vector<double> jtPfMUMFCutLow = cutPropPbPb.GetJtPfMUMFCutLow();
  std::vector<double> jtPfMUMFCutHi = cutPropPbPb.GetJtPfMUMFCutHi();

  std::vector<int> jetPPMatchedToPbPb;
  for(Int_t jI = 0; jI < nJtPbPb; ++jI){
    const Int_t rValPbPb = getRVal(jetPbPbList.at(jI));

    bool isMatched = false;
    for(Int_t jI2 = 0; jI2 < nJtPP; ++jI2){
      const Int_t rValPP = getRVal(jetPPList.at(jI2));
      
      if(rValPbPb == rValPP){
	jetPPMatchedToPbPb.push_back(jI2);
	isMatched = true;
	break;
      }
    }

    if(!isMatched){
      std::cout << "Warning: " << jetPbPbList.at(jI) << " has no match." << std::endl;
      std::cout << " Options: ";
      for(Int_t jI2 = 0; jI2 < nJtPP; ++jI2){
	std::cout << jetPPList.at(jI2) << ", ";
      }
      std::cout << std::endl;
    }

  }

  if(jetPPMatchedToPbPb.size() != jetPbPbList.size()){
    std::cout << "Size of pp algos, " << jetPPMatchedToPbPb.size() << ",  matched to jetPbPbList of different size, " << jetPbPbList.size() << ". return 1" << std::endl;
    
    inFilePP_p->Close();
    delete inFilePP_p;
    
    inFilePbPb_p->Close();
    delete inFilePbPb_p;    

    return 1;
  }

  std::cout << "Processing the following pairs of jets: " << std::endl;
  for(Int_t jI = 0; jI < nJtPbPb; ++jI){
    std::cout << " " << jI << "/" << nJtPbPb << ": " << jetPbPbList.at(jI) << ", " << jetPPList.at(jetPPMatchedToPbPb.at(jI)) << std::endl;
  }

  if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  if(doLocalDebug || doGlobalDebug){
    std::cout << "nJtPbPb*nCentBins*nID*nResponseMod*nJtAbsEtaBins*nSyst*nBayes=Total" << std::endl;
    Int_t Total = nJtPbPb*nCentBins*nID*nResponseMod*nJtAbsEtaBins*nSyst*nBayes;
    std::cout << nJtPbPb << "*" << nCentBins << "*" << nID << "*" << nResponseMod << "*" << nJtAbsEtaBins << "*" << nSyst << "*" << nBayes << "=" << Total << std::endl;
  }


  std::string outFileName = inFileNamePP;
  while(outFileName.find("/") != std::string::npos){outFileName.replace(0, outFileName.find("/")+1, "");}
  if(outFileName.find(".txt") != std::string::npos) outFileName.replace(outFileName.find(".txt"), std::string(".txt").size(), "");
  else if(outFileName.find(".root") != std::string::npos) outFileName.replace(outFileName.find(".root"), std::string(".root").size(), "");

  std::string outFileName2 = inFileNamePbPb;
  while(outFileName2.find("/") != std::string::npos){outFileName2.replace(0, outFileName2.find("/")+1, "");}
  if(outFileName2.find(".txt") != std::string::npos) outFileName2.replace(outFileName2.find(".txt"), std::string(".txt").size(), "");
  else if(outFileName2.find(".root") != std::string::npos) outFileName2.replace(outFileName2.find(".root"), std::string(".root").size(), "");

  const Int_t sizeToTruncName = 40;
  while(outFileName.size() > sizeToTruncName){outFileName = outFileName.substr(0,outFileName.size()-1);}
  while(outFileName2.size() > sizeToTruncName){outFileName2 = outFileName2.substr(0,outFileName2.size()-1);}
  outFileName = "output/" + outFileName + "_" + outFileName2 + "_PlotUnfoldedSpectra_" + dateStr + ".root";


  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  //Following two lines necessary else you get absolutely clobbered on deletion timing. See threads here:
  //https://root-forum.cern.ch/t/closing-root-files-is-dead-slow-is-this-a-normal-thing/5273/16
  //https://root-forum.cern.ch/t/tfile-speed/17549/25
  //Bizarre
  outFile_p->SetBit(TFile::kDevNull);
  TH1::AddDirectory(kFALSE);
 
  TDirectory* dirPP_p[nJtPP];
  TDirectory* dirPbPb_p[nJtPbPb];
  
  for(Int_t tI = 0; tI < nJtPbPb; ++tI){
    std::string tempStr = jetPbPbList.at(tI);
    if(tempStr.find("/") != std::string::npos) tempStr.replace(tempStr.find("/"), tempStr.size() - tempStr.find("/"), "");
    dirPbPb_p[tI] = (TDirectory*)outFile_p->mkdir(tempStr.c_str());
    dirPbPb_p[tI]->cd();
  }

  TH1D* jtPtUnfolded_RecoTrunc_PbPb_h[nJtPbPb][nCentBins][nID][nResponseMod][nJtAbsEtaBins][nSyst][nBayes];
  

  for(Int_t tI = 0; tI < nJtPP; ++tI){
    std::string tempStr = jetPPList.at(tI);
    if(tempStr.find("/") != std::string::npos) tempStr.replace(tempStr.find("/"), tempStr.size() - tempStr.find("/"), "");
    dirPP_p[tI] = (TDirectory*)outFile_p->mkdir(tempStr.c_str());
    dirPP_p[tI]->cd();
  }

  TH1D* jtPtUnfolded_RecoTrunc_PP_h[nJtPP][nID][nResponseMod][nJtAbsEtaBins][nSyst][nBayes];

  if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  std::cout << "Grabbing histograms..." << std::endl;

  if(doLocalDebug || doGlobalDebug) std::cout << "idI,mI,aI,sI,bI,tI,cI" << std::endl;

  cppWatch histGrabTot;
  cppWatch histGrab;
  cppWatch histClone;

  histGrabTot.start();

  inFilePbPb_p->cd();
  for(unsigned int tI = 0; tI < jetPbPbList.size(); ++tI){
    cppWatch histGrabLocal;
    histGrabLocal.start();
    histGrab.start();
    std::cout << " Grabbing PbPb " << tI << "/" << jetPbPbList.size() << ": " << jetPbPbList.at(tI) << std::endl;

    TDirectory* dir_p = (TDirectory*)inFilePbPb_p->Get(jetPbPbList.at(tI).c_str());
    TIter nextPbPb(dir_p->GetListOfKeys());
    TKey* key = NULL;

    while( (key = (TKey*)nextPbPb()) ){
      const std::string name = key->GetName();
      const std::string className = key->GetClassName();

      if(className.find("TH1") == std::string::npos) continue;

      int idPos = -1;
      int modPos = -1;
      int absEtaPos = -1;
      int systPos = -1;
      int bayesPos = -1;
      int centPos = -1;

      for(unsigned int idI = 0; idI < idStr.size(); ++idI){
	if(name.find("_" + idStr.at(idI) + "_") != std::string::npos){
	  idPos = idI;
	  break;
	}
      }

      for(int mI = 0; mI < nResponseMod; ++mI){
	if(name.find("_ResponseMod" + prettyString(responseMod.at(mI), 2, true) + "_") != std::string::npos){
	  modPos = mI;
	  break;
	}
      }

      for(int aI = 0; aI < nJtAbsEtaBins; ++aI){
	const std::string jtAbsEtaStr = "AbsEta" + prettyString(jtAbsEtaBinsLow.at(aI), 1, true) + "to" + prettyString(jtAbsEtaBinsHi.at(aI), 1, true);	
	if(name.find(jtAbsEtaStr) != std::string::npos){
	  absEtaPos = aI;
	  break;
	}
      }

      for(Int_t sI = 1; sI < nSystOrig; ++sI){
	if(name.find(systStr.at(sI)) != std::string::npos){
	  systPos = sI;
	  break;
	}
      }
      if(systPos < 0) systPos = 0;

      for(Int_t bI = 0; bI < nBayes; ++bI){
	const std::string bayesStr = "_Bayes" + std::to_string(bayesVal.at(bI)) + "_";
	if(name.find(bayesStr) != std::string::npos){
	  bayesPos = bI;
	  break;
	}
      }

      for(Int_t cI = 0; cI < nCentBins; ++cI){		
	const std::string centStr = "Cent" + std::to_string(centBinsLow.at(cI)) + "to" + std::to_string(centBinsHi.at(cI));
	if(name.find(centStr) != std::string::npos){
	  centPos = cI;
	  break;
	}
      }


      if(idPos < 0) std::cout << "Warning: Cannot find idPos for \'" << name << "\'. Code will likely break" << std::endl;
      if(modPos < 0) std::cout << "Warning: Cannot find modPos for \'" << name << "\'. Code will likely break" << std::endl;
      if(absEtaPos < 0) std::cout << "Warning: Cannot find absEtaPos for \'" << name << "\'. Code will likely break" << std::endl;
      if(systPos < 0) std::cout << "Warning: Cannot find systPos for \'" << name << "\'. Code will likely break" << std::endl;
      if(bayesPos < 0) std::cout << "Warning: Cannot find bayesPos for \'" << name << "\'. Code will likely break" << std::endl;
      if(centPos < 0) std::cout << "Warning: Cannot find centPos for \'" << name << "\'. Code will likely break" << std::endl;

      jtPtUnfolded_RecoTrunc_PbPb_h[tI][centPos][idPos][modPos][absEtaPos][systPos][bayesPos] = (TH1D*)key->ReadObj();   
      /*

      if(systPos == 0){
	histGrab.stop();
	histGrabLocal.stop();
	histCloneLocal.start();
	histClone.start();

	const std::string centStr = "Cent" + std::to_string(centBinsLow.at(centPos)) + "to" + std::to_string(centBinsHi.at(centPos));
	const std::string bayesStr = "Bayes" + std::to_string(bayesVal.at(bayesPos));
	const std::string jtAbsEtaStr = "AbsEta" + prettyString(jtAbsEtaBinsLow.at(absEtaPos), 1, true) + "to" + prettyString(jtAbsEtaBinsHi.at(absEtaPos), 1, true);	
	const std::string resStr = "ResponseMod" + prettyString(responseMod.at(modPos), 2, true);

      	for(Int_t sI = nSystOrig; sI < nSyst; ++sI){
	  const std::string tempSystStr = systStr.at(sI) + "_";
	  const std::string newName = "jtPtUnfolded_RecoTrunc_" + jetPbPbList.at(tI) + "_" + centStr + "_" + idStr.at(idPos) + "_" + resStr + "_" + jtAbsEtaStr + "_" + tempSystStr + bayesStr + "_h";
	  
	  jtPtUnfolded_RecoTrunc_PbPb_h[tI][centPos][idPos][modPos][absEtaPos][sI][bayesPos] = (TH1D*) (((TH1D*)key->ReadObj())->Clone(newName.c_str()));   
	}

	histClone.stop();
	histCloneLocal.stop();

	histGrab.start();
	histGrabLocal.start();
      }
      */
    }    

    histGrabLocal.stop();
    histGrab.stop();

    cppWatch histCloneLocal;
    histCloneLocal.start();
    histClone.start();
    
    std::cout << "  Cloning..." << std::endl;

    for(Int_t cI = 0; cI < nCentBins; ++cI){
      const std::string centStr = "Cent" + std::to_string(centBinsLow.at(cI)) + "to" + std::to_string(centBinsHi.at(cI));

      for(unsigned int idI = 0; idI < idStr.size(); ++idI){
	for(Int_t mI = 0; mI < nResponseMod; ++mI){
	  const std::string resStr = "ResponseMod" + prettyString(responseMod.at(mI), 2, true);
	  
	  for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
	    const std::string jtAbsEtaStr = "AbsEta" + prettyString(jtAbsEtaBinsLow.at(aI), 1, true) + "to" + prettyString(jtAbsEtaBinsHi.at(aI), 1, true);	
	    
	    for(Int_t sI = nSystOrig; sI < nSyst; ++sI){
	      const std::string tempSystStr = systStr.at(sI) + "_";
	      
	      for(Int_t bI = 0; bI < nBayes; ++bI){
		const std::string bayesStr = "Bayes" + std::to_string(bayesVal.at(bI));
		
		const std::string newName = "jtPtUnfolded_RecoTrunc_" + jetPbPbList.at(tI) + "_" + centStr + "_" + idStr.at(idI) + "_" + resStr + "_" + jtAbsEtaStr + "_" + tempSystStr + bayesStr + "_h";
		
		jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idI][mI][aI][sI][bI] = (TH1D*)jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idI][mI][aI][0][bI]->Clone(newName.c_str());
	      }
	    }
	  }
	}
      }
    }

    histCloneLocal.stop();
    histClone.stop();

    std::cout << "  Local histgrab, clone times: "<< histGrabLocal.total() << ", " << histCloneLocal.total() << std::endl;
  }

  

  inFilePP_p->cd();
  for(unsigned int tI = 0; tI < jetPPList.size(); ++tI){
    cppWatch histGrabLocal;
    histGrabLocal.start();
    histGrab.start();
    std::cout << " Grabbing PP " << tI << "/" << jetPPList.size() << ": " << jetPPList.at(tI) << std::endl;

    TDirectory* dir_p = (TDirectory*)inFilePP_p->Get(jetPPList.at(tI).c_str());
    TIter nextPP(dir_p->GetListOfKeys());
    TKey* key = NULL;

    while( (key = (TKey*)nextPP()) ){
      const std::string name = key->GetName();
      const std::string className = key->GetClassName();

      if(className.find("TH1") == std::string::npos) continue;

      int idPos = -1;
      int modPos = -1;
      int absEtaPos = -1;
      int systPos = -1;
      int bayesPos = -1;

      for(unsigned int idI = 0; idI < idStr.size(); ++idI){
	if(name.find("_" + idStr.at(idI) + "_") != std::string::npos){
	  idPos = idI;
	  break;
	}
      }

      for(int mI = 0; mI < nResponseMod; ++mI){
	if(name.find("_ResponseMod" + prettyString(responseMod.at(mI), 2, true) + "_") != std::string::npos){
	  modPos = mI;
	  break;
	}
      }

      for(int aI = 0; aI < nJtAbsEtaBins; ++aI){
	const std::string jtAbsEtaStr = "AbsEta" + prettyString(jtAbsEtaBinsLow.at(aI), 1, true) + "to" + prettyString(jtAbsEtaBinsHi.at(aI), 1, true);	
	if(name.find(jtAbsEtaStr) != std::string::npos){
	  absEtaPos = aI;
	  break;
	}
      }

      for(Int_t sI = 1; sI < nSystOrig; ++sI){
	if(name.find(systStr.at(sI)) != std::string::npos){
	  systPos = sI;
	  break;
	}
      }
      if(systPos < 0) systPos = 0;

      for(Int_t bI = 0; bI < nBayes; ++bI){
	const std::string bayesStr = "_Bayes" + std::to_string(bayesVal.at(bI)) + "_";
	if(name.find(bayesStr) != std::string::npos){
	  bayesPos = bI;
	  break;
	}
      }

      if(idPos < 0) std::cout << "Warning: Cannot find idPos for \'" << name << "\'. Code will likely break" << std::endl;
      if(modPos < 0) std::cout << "Warning: Cannot find modPos for \'" << name << "\'. Code will likely break" << std::endl;
      if(absEtaPos < 0) std::cout << "Warning: Cannot find absEtaPos for \'" << name << "\'. Code will likely break" << std::endl;
      if(systPos < 0) std::cout << "Warning: Cannot find systPos for \'" << name << "\'. Code will likely break" << std::endl;
      if(bayesPos < 0) std::cout << "Warning: Cannot find bayesPos for \'" << name << "\'. Code will likely break" << std::endl;

      jtPtUnfolded_RecoTrunc_PP_h[tI][idPos][modPos][absEtaPos][systPos][bayesPos] = (TH1D*)key->ReadObj();   
    }

    histGrabLocal.stop();
    histGrab.stop();
    
    cppWatch histCloneLocal;
    histCloneLocal.start();
    histClone.start();

    std::cout << "  Cloning..." << std::endl;

    for(unsigned int idI = 0; idI < idStr.size(); ++idI){
      for(Int_t mI = 0; mI < nResponseMod; ++mI){
	const std::string resStr = "ResponseMod" + prettyString(responseMod.at(mI), 2, true);
	
	for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
	  const std::string jtAbsEtaStr = "AbsEta" + prettyString(jtAbsEtaBinsLow.at(aI), 1, true) + "to" + prettyString(jtAbsEtaBinsHi.at(aI), 1, true);	
	  
	  for(Int_t sI = nSystOrig; sI < nSyst; ++sI){
	    const std::string tempSystStr = systStr.at(sI) + "_";
	    
	    for(Int_t bI = 0; bI < nBayes; ++bI){
	      const std::string bayesStr = "Bayes" + std::to_string(bayesVal.at(bI));
	      
	      const std::string newName = "jtPtUnfolded_RecoTrunc_" + jetPPList.at(tI) + "_PP_" + idStr.at(idI) + "_" + resStr + "_" + jtAbsEtaStr + "_" + tempSystStr + bayesStr + "_h";
	      
	      jtPtUnfolded_RecoTrunc_PP_h[tI][idI][mI][aI][sI][bI] = (TH1D*)jtPtUnfolded_RecoTrunc_PP_h[tI][idI][mI][aI][0][bI]->Clone(newName.c_str());
	    }
	  }
	}
      }
    }

    histCloneLocal.stop();
    histClone.stop();

    std::cout << "  Local histgrab, clone times: " << histGrabLocal.total() << ", " << histCloneLocal.total() << std::endl;
  }

  
  histGrabTot.stop();
  std::cout << "Total histGrabTot: " << histGrabTot.total() << std::endl;
  std::cout << "  just grabbing: " << histGrab.total() << std::endl;
  std::cout << "  just cloning: " << histClone.total() << std::endl;


  if(doLocalDebug || doGlobalDebug){
    std::cout << "Started return" << std::endl;
    return 1;
  }


  const Double_t lumiFactor = getLumiFactor();
  const Double_t nMBEvents = getNMBEvents();

  std::cout << "Scaling histograms..." << std::endl;

  for(Int_t idI = 0; idI < nID; ++idI){
    for(Int_t mI = 0; mI < nResponseMod; ++mI){
      for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
	const Double_t etaBinWidth = TMath::Abs(jtAbsEtaBinsHi[aI] - jtAbsEtaBinsLow[aI]);

	for(Int_t sI = 0; sI < nSyst; ++sI){
	  for(Int_t bI = 0; bI < nBayes; ++bI){
	    for(Int_t tI = 0; tI < nJtPbPb; ++tI){
	      for(Int_t cI = 0; cI < nCentBins; ++cI){		
		const std::string centStr = "Cent" + std::to_string(centBinsLow.at(cI)) + "to" + std::to_string(centBinsHi.at(cI));
		const Double_t centBinWidth = TMath::Abs(centBinsHi[cI] - centBinsLow[cI])/100.;

		Double_t tempTAAFactor = getTAAScaleFactor(centStr);
		if(isStrSame(systStr[sI], "TAAUp")) tempTAAFactor += tempTAAFactor*getTAAScaleFactorUp(centStr);
		else if(isStrSame(systStr[sI], "TAADown")) tempTAAFactor -= tempTAAFactor*getTAAScaleFactorDown(centStr);

		Double_t totalPbPbFactor = tempTAAFactor*nMBEvents*2.*etaBinWidth*centBinWidth;

		jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idI][mI][aI][sI][bI]->Scale(1./totalPbPbFactor);
		divBinWidth(jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idI][mI][aI][sI][bI]);
	      }
	    }

	    for(Int_t tI = 0; tI < nJtPP; ++tI){
	      Double_t tempLumiFactor = lumiFactor;
	      if(isStrSame(systStr[sI], "LumiUp")) tempLumiFactor += tempLumiFactor*getLumiPercentError();
	      else if(isStrSame(systStr[sI], "LumiDown")) tempLumiFactor -= tempLumiFactor*getLumiPercentError();	      

	      Double_t totalPPFactor = tempLumiFactor*2.*etaBinWidth;
	      jtPtUnfolded_RecoTrunc_PP_h[tI][idI][mI][aI][sI][bI]->Scale(1./totalPPFactor);
	      divBinWidth(jtPtUnfolded_RecoTrunc_PP_h[tI][idI][mI][aI][sI][bI]);
	    }
	  }
	}
      }
    }
  }

  std::cout << "Writing histograms..." << std::endl;

  outFile_p->cd();
  for(Int_t tI = 0; tI < nJtPbPb; ++tI){
    dirPbPb_p[tI]->cd();

    for(Int_t cI = 0; cI < nCentBins; ++cI){		
      for(Int_t idI = 0; idI < nID; ++idI){
	for(Int_t mI = 0; mI < nResponseMod; ++mI){
	  for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
	    for(Int_t sI = 0; sI < nSyst; ++sI){
	      for(Int_t bI = 0; bI < nBayes; ++bI){
		std::string newName = jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idI][mI][aI][sI][bI]->GetName();
		newName.replace(newName.find("_h"), 2, "_Rescaled_h");
		jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idI][mI][aI][sI][bI]->Write(newName.c_str(), TObject::kOverwrite);
	      }
	    }
	  }
	}
      }
    }
  }

  for(Int_t tI = 0; tI < nJtPP; ++tI){
    dirPP_p[tI]->cd();
    
    for(Int_t idI = 0; idI < nID; ++idI){
      for(Int_t mI = 0; mI < nResponseMod; ++mI){
	for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
	  for(Int_t sI = 0; sI < nSyst; ++sI){
	    for(Int_t bI = 0; bI < nBayes; ++bI){
	      std::string newName = jtPtUnfolded_RecoTrunc_PP_h[tI][idI][mI][aI][sI][bI]->GetName();
	      newName.replace(newName.find("_h"), 2, "_Rescaled_h");
	      jtPtUnfolded_RecoTrunc_PP_h[tI][idI][mI][aI][sI][bI]->Write(newName.c_str(), TObject::kOverwrite);
	    }
	  }
	}
      }
    }
  }

  std::cout << "Plotting... " << std::endl;
  const std::string plotID = "LightMUAndCHID";
  const std::string plotAbsEtaStr = "AbsEta0p0to2p0";
  const Int_t plotBayesVal = 4;

  Int_t idPos = -1;
  Int_t absEtaPos = -1;
  Int_t bayesPos = -1;

  for(Int_t idI = 0; idI < nID; ++idI){
    if(isStrSame(idStr.at(idI), plotID)){
      idPos = idI;
      break;
    }
  }

  for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
    const std::string jtAbsEtaStr = "AbsEta" + prettyString(jtAbsEtaBinsLow.at(aI), 1, true) + "to" + prettyString(jtAbsEtaBinsHi.at(aI), 1, true);	
    
    if(isStrSame(jtAbsEtaStr, plotAbsEtaStr)){
      absEtaPos = aI;
      break;
    }
  }

  for(Int_t bI = 0; bI < nBayes; ++bI){
    if(bayesVal.at(bI) == plotBayesVal){
      bayesPos = bI;
      break;
    }
  }
  
  if(idPos == -1 || absEtaPos == -1 || bayesPos == -1){
    std::cout << "Plotting pos not found, skip plotting, return 1" << std::endl;

    std::cout << "Closing files..." << std::endl;
    
    outFile_p->Close();
    delete outFile_p;
  
    inFilePP_p->Close();
    delete inFilePP_p;

    inFilePbPb_p->Close();
    delete inFilePbPb_p;

    return 1;
  }

  const std::string idNameStr = idStr.at(idPos);
  const std::string plotBayesStr = "Bayes" + std::to_string(plotBayesVal);

  const Int_t nSystSmooth = 2;
  const std::string systSmooth[nSystSmooth] = {"Unsmoothed", "Smoothed"};
  const bool systSmoothBool[nSystSmooth] = {false, true};

  for(Int_t tI = 0; tI < nJtPbPb; ++tI){
    const Int_t rVal = getRVal(jetPbPbList.at(tI));
    const Int_t ppPos = jetPPMatchedToPbPb.at(tI);
    const std::string rValStr = std::to_string(rVal);

    for(Int_t mI = 0; mI < nResponseMod; ++mI){
      const std::string responseStr = prettyString(responseMod.at(mI), 2, true);

      for(Int_t usI = 0; usI < nSystSmooth; ++usI){
	TLatex* label_p = new TLatex();
	label_p->SetTextFont(43);
	label_p->SetTextSize(16);

	Int_t nYBins = 20;
	Double_t yBins[nYBins];
	
	TCanvas* spectCanv_p = new TCanvas("spectCanv_p", "", 450, 450);
	spectCanv_p->SetTopMargin(0.08);
	spectCanv_p->SetRightMargin(0.01);
	spectCanv_p->SetLeftMargin(0.14);
	spectCanv_p->SetBottomMargin(0.14);

	Double_t xMinVal = 200;
	Double_t xPointMinVal = xMinVal;
	if(getRVal(jetPbPbList.at(tI)) >= 8) xPointMinVal = 300;
	Double_t xMaxVal = 1200;
	Double_t xPointMaxVal = 1000;

	const std::string yAxisTitle = "#frac{d^{2}#sigma^{pp}}{dp_{T}d#eta} or #frac{1}{N_{evt}}#frac{1}{#LTTAA#GT}#frac{d^{2}N^{PbPb}}{dp_{T}d#eta} [#frac{nb}{GeV/c}]";
	TH1F* tempHist_p = new TH1F("tempHist_p", (";Jet p_{T} (GeV/c);" + yAxisTitle).c_str(), 10, xMinVal, xMaxVal);
	centerTitles(tempHist_p);

	tempHist_p->GetXaxis()->SetTitleOffset(1.8); 
	tempHist_p->GetYaxis()->SetTitleOffset(1.6);
	
	for(Int_t bIX = 0; bIX < jtPtUnfolded_RecoTrunc_PP_h[ppPos][idPos][mI][absEtaPos][0][bayesPos]->GetNbinsX(); ++bIX){
	  bool binIsBad = jtPtUnfolded_RecoTrunc_PP_h[ppPos][idPos][mI][absEtaPos][0][bayesPos]->GetBinCenter(bIX+1) < xPointMinVal;
	  binIsBad = binIsBad || jtPtUnfolded_RecoTrunc_PP_h[ppPos][idPos][mI][absEtaPos][0][bayesPos]->GetBinCenter(bIX+1) > xPointMaxVal;
	  if(binIsBad){
	    jtPtUnfolded_RecoTrunc_PP_h[ppPos][idPos][mI][absEtaPos][0][bayesPos]->SetBinContent(bIX+1, 0.0);
	    jtPtUnfolded_RecoTrunc_PP_h[ppPos][idPos][mI][absEtaPos][0][bayesPos]->SetBinError(bIX+1, 0.0);
	  }
	}

	Double_t minVal = getMinGTZero(jtPtUnfolded_RecoTrunc_PP_h[ppPos][idPos][mI][absEtaPos][0][bayesPos]);
	Double_t maxVal = getMax(jtPtUnfolded_RecoTrunc_PP_h[ppPos][idPos][mI][absEtaPos][0][bayesPos]);

	for(Int_t cI = 0; cI < nCentBins; ++cI){
	  for(Int_t sI = 0; sI < nSyst; ++sI){
	    scaleCentralAndErrorValues(jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idPos][mI][absEtaPos][sI][bayesPos], centBinsScalingFact.at(cI));
	  }

	  for(Int_t bIX = 0; bIX < jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idPos][mI][absEtaPos][0][bayesPos]->GetNbinsX(); ++bIX){
	    bool binIsBad = jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idPos][mI][absEtaPos][0][bayesPos]->GetBinCenter(bIX+1) < xPointMinVal;
	    binIsBad = binIsBad || jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idPos][mI][absEtaPos][0][bayesPos]->GetBinCenter(bIX+1) > xPointMaxVal;
	    if(binIsBad){
	      jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idPos][mI][absEtaPos][0][bayesPos]->SetBinContent(bIX+1, 0.0);
	      jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idPos][mI][absEtaPos][0][bayesPos]->SetBinError(bIX+1, 0.0);
	    }
	  }

	  
	  Double_t tempMinVal = getMinGTZero(jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idPos][mI][absEtaPos][0][bayesPos]);
	  Double_t tempMaxVal = getMax(jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idPos][mI][absEtaPos][0][bayesPos]);
	  if(minVal > tempMinVal) minVal = tempMinVal;
	  if(maxVal < tempMaxVal) maxVal = tempMaxVal;
	}
	
	maxVal *= 100.;
	minVal /= 100.;

	tempHist_p->SetMaximum(maxVal);
	tempHist_p->SetMinimum(minVal);

	tempHist_p->GetYaxis()->SetLabelFont(43);
	tempHist_p->GetYaxis()->SetLabelSize(11);
	
	std::cout << "Title offset: " << tempHist_p->GetYaxis()->GetTitleOffset() << std::endl;
	
	getLogBins(minVal, maxVal, nYBins, yBins);
	
	tempHist_p->DrawCopy("HIST E1 P");
	
	TLegend* leg_p = new TLegend(0.7, 0.65, 0.95, 0.88);
	leg_p->SetBorderSize(0);
	leg_p->SetFillColor(0);
	leg_p->SetFillStyle(0);
	leg_p->SetTextFont(43);
	leg_p->SetTextSize(16);


	jtPtUnfolded_RecoTrunc_PP_h[ppPos][idPos][mI][absEtaPos][0][bayesPos]->SetMarkerColor(kPalette.getColor(getColorPosFromCent("", true)));
	jtPtUnfolded_RecoTrunc_PP_h[ppPos][idPos][mI][absEtaPos][0][bayesPos]->SetLineColor(kPalette.getColor(getColorPosFromCent("", true)));
	jtPtUnfolded_RecoTrunc_PP_h[ppPos][idPos][mI][absEtaPos][0][bayesPos]->SetMarkerStyle(getStyleFromCent("", true));
	jtPtUnfolded_RecoTrunc_PP_h[ppPos][idPos][mI][absEtaPos][0][bayesPos]->SetMarkerSize(1);

	std::vector<TH1D*> systHistVectPP;
	systHistVectPP.reserve(nSyst-1);
	for(Int_t sI = 1; sI < nSyst; ++sI){
	  systHistVectPP.push_back((TH1D*)jtPtUnfolded_RecoTrunc_PP_h[ppPos][idPos][mI][absEtaPos][sI][bayesPos]->Clone(("systPP_" + systStr.at(sI)).c_str()));
	}      
	
	pdfPerSlide.push_back({});
	std::vector<double> systValVectPP = getSyst(jtPtUnfolded_RecoTrunc_PP_h[ppPos][idPos][mI][absEtaPos][0][bayesPos], systHistVectPP, systStr, xPointMinVal, xPointMaxVal, &(pdfPerSlide.at(pdfPerSlide.size()-1)), systSmoothBool[usI]);
    
	slideTitles.push_back("PP Systematic (" + jetPPList.at(ppPos) + ", " + responseStr +  ")");

	for(unsigned int sI = 0; sI < systHistVectPP.size(); ++sI){
	  delete systHistVectPP.at(sI);
	}
	drawSyst(spectCanv_p, jtPtUnfolded_RecoTrunc_PP_h[ppPos][idPos][mI][absEtaPos][0][bayesPos], systValVectPP, xPointMinVal, xPointMaxVal);


	for(Int_t cI = 0; cI < nCentBins; ++cI){
	  const std::string centStr = "Cent" + std::to_string(centBinsLow.at(cI)) + "to" + std::to_string(centBinsHi.at(cI));
	  const std::string centStr2 = std::to_string(centBinsLow.at(cI)) + "-" + std::to_string(centBinsHi.at(cI)) + "%";

	  jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idPos][mI][absEtaPos][0][bayesPos]->SetMarkerColor(kPalette.getColor(getColorPosFromCent(centStr, false)));
	  jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idPos][mI][absEtaPos][0][bayesPos]->SetLineColor(kPalette.getColor(getColorPosFromCent(centStr, false)));
	  jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idPos][mI][absEtaPos][0][bayesPos]->SetMarkerStyle(getStyleFromCent(centStr, false));
	  jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idPos][mI][absEtaPos][0][bayesPos]->SetMarkerSize(1);

	  std::vector<TH1D*> systHistVectPbPb;
	  systHistVectPbPb.reserve(nSyst-1);
	  for(Int_t sI = 1; sI < nSyst; ++sI){
	    systHistVectPbPb.push_back((TH1D*)jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idPos][mI][absEtaPos][sI][bayesPos]->Clone(("systPbPb_" + systStr.at(sI)).c_str()));
	  }      

	  pdfPerSlide.push_back({});
	  std::vector<double> systValVectPbPb = getSyst(jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idPos][mI][absEtaPos][0][bayesPos], systHistVectPbPb, systStr, xPointMinVal, xPointMaxVal, &(pdfPerSlide.at(pdfPerSlide.size()-1)), systSmoothBool[usI]);
	  slideTitles.push_back(centStr2 + " Systematic (" + jetPbPbList.at(tI) + ", " + responseStr +  ")");
	  
	  for(unsigned int sI = 0; sI < systHistVectPbPb.size(); ++sI){
	    delete systHistVectPbPb.at(sI);
	  }
	  drawSyst(spectCanv_p, jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idPos][mI][absEtaPos][0][bayesPos], systValVectPbPb, xPointMinVal, xPointMaxVal);	
	}


	jtPtUnfolded_RecoTrunc_PP_h[ppPos][idPos][mI][absEtaPos][0][bayesPos]->DrawCopy("HIST E1 P SAME");

	for(Int_t cI = 0; cI < nCentBins; ++cI){
	  jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idPos][mI][absEtaPos][0][bayesPos]->DrawCopy("HIST E1 P SAME");
	}


	for(Int_t cI = nCentBins-1; cI >= 0; --cI){
	  const std::string centLegStr = std::to_string(centBinsLow.at(cI)) + "-" + std::to_string(centBinsHi.at(cI)) + "% x 10^{" + std::to_string(cI+1) + "}";
	  jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idPos][mI][absEtaPos][0][bayesPos]->SetFillColorAlpha(jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idPos][mI][absEtaPos][0][bayesPos]->GetMarkerColor(), .25);      
	  leg_p->AddEntry(jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idPos][mI][absEtaPos][0][bayesPos], centLegStr.c_str(), "P L F");	
	}

	jtPtUnfolded_RecoTrunc_PP_h[ppPos][idPos][mI][absEtaPos][0][bayesPos]->SetFillColorAlpha(jtPtUnfolded_RecoTrunc_PP_h[ppPos][idPos][mI][absEtaPos][0][bayesPos]->GetMarkerColor(), .25);      
	leg_p->AddEntry(jtPtUnfolded_RecoTrunc_PP_h[ppPos][idPos][mI][absEtaPos][0][bayesPos], "pp", "P L F");
      
	gPad->SetLogy();
	gPad->SetLogx();
	gStyle->SetOptStat(0);

	drawWhiteBox(900, 1100, minVal/10, minVal*.9);

	label_p->DrawLatex(200, yBins[0]/4., "200");
	label_p->DrawLatex(400, yBins[0]/4., "400");
	label_p->DrawLatex(600, yBins[0]/4., "600");
	label_p->DrawLatex(800, yBins[0]/4., "800");
	label_p->DrawLatex(1000, yBins[0]/4., "1000");

	std::string rStr = "R=";
	if(getRVal(jetPbPbList.at(tI)) < 10) rStr = rStr + "0." + std::to_string(getRVal(jetPbPbList.at(tI)));
	else rStr = rStr + "1.0";
    
	label_p->DrawLatex(250, yBins[19], ("anti-k_{T} " + rStr).c_str());
	label_p->DrawLatex(250, yBins[18], "|#eta_{jets}| < 2");

	label_p->SetNDC();
      
	label_p->DrawLatex(0.1, 0.95, "#bf{CMS Preliminary}");
	label_p->DrawLatex(0.4, 0.95, "27.4 pb^{-1} pp + 404 #mub^{-1} PbPb (5.02 TeV)");

	leg_p->Draw("SAME");

	gPad->RedrawAxis();
	gPad->SetTicks(1,2);

	std::string saveName = "spectra_" + jetPbPbList.at(tI) + "_R" + rValStr + "_" + idNameStr + "_" + plotAbsEtaStr + "_" + plotBayesStr + "_" + responseStr + "_" + systSmooth[usI] + "_" + dateStr + ".pdf";
	slideTitles.push_back("Spectra (" + responseStr + ")");
	pdfPerSlide.push_back({});
	pdfPerSlide.at(pdfPerSlide.size()-1).push_back(saveName);

	saveName = "pdfDir/" + dateStr + "/" + saveName;

	quietSaveAs(spectCanv_p, saveName);
	
	delete spectCanv_p;
	delete label_p;
	delete leg_p;

	for(Int_t cI = 0; cI < nCentBins; ++cI){
	  for(Int_t sI = 0; sI < nSyst; ++sI){
	    scaleCentralAndErrorValues(jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idPos][mI][absEtaPos][sI][bayesPos], 1./centBinsScalingFact.at(cI));
	  }
	}
      }
    }
  }

  std::cout << "Closing files..." << std::endl;
  
  outFile_p->Close();
  delete outFile_p;
  
  inFilePP_p->Close();
  delete inFilePP_p;

  inFilePbPb_p->Close();
  delete inFilePbPb_p;

  texSlideCreator tex;
  tex.Clean();
  tex.Init(outFileName);
  tex.SetAuthor("Christopher McGinn");
  tex.SetSlideTitles(slideTitles);
  tex.SetSlidePdfs(pdfPerSlide);
  if(!(tex.CreateTexSlides())){
    std::cout << "Warning: .tex slide creation failed" << std::endl;
  }

  std::cout << "Job complete" << std::endl;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 3){
    std::cout << "Usage: ./bin/plotUnfoldedSpectra.exe <inFileNamePP> <inFileNamePbPb>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += plotUnfoldedSpectra(argv[1], argv[2]);
  return retVal;
}
