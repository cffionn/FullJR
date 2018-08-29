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
#include "Utility/include/getLinBins.h"
#include "Utility/include/getLogBins.h"
#include "Utility/include/histDefUtility.h"
#include "Utility/include/kirchnerPalette.h"
#include "Utility/include/vanGoghPalette.h"
#include "Utility/include/lumiAndTAAUtil.h"
#include "Utility/include/returnRootFileContentsList.h"

std::vector<double> getSyst(TH1D* nominal_p, std::vector<TH1D*> syst_p, std::vector<std::string> systStr, Double_t minXVal, Double_t maxXVal, std::vector<std::string>* plotNames, bool doSmoothing, std::vector<std::vector<std::string > > systToCombo)
{
  std::cout << "Starting getSyst: " << std::endl;
  std::cout << "SystHist input: " << std::endl;
  for(unsigned int i = 0; i < syst_p.size(); ++i){
    std::cout << syst_p.at(i)->GetName() << std::endl;
  }
  for(unsigned int sI = 0; sI < systStr.size(); ++sI){
    std::cout << systStr.at(sI) << std::endl;
  }

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
  std::vector<double> nomValVect;

  //Building parallel systematic histograms
  Int_t tempNBins = 0;
  for(Int_t bIX = 0; bIX < nominal_p->GetNbinsX(); ++bIX){
    Double_t binCenter = nominal_p->GetBinCenter(bIX+1);
    if(binCenter > maxXVal || binCenter < minXVal) continue;
    
    tempNBins++;
  }

  const Int_t nBins = tempNBins;
  Double_t bins[nBins+1];

  Int_t binFillPos = 0;
  for(Int_t bIX = 0; bIX < nominal_p->GetNbinsX(); ++bIX){
    Double_t binCenter = nominal_p->GetBinCenter(bIX+1);
    if(binCenter < minXVal) continue;
    if(binCenter > maxXVal){
      bins[binFillPos] = nominal_p->GetBinLowEdge(bIX+1);
      break;
    }

    bins[binFillPos] = nominal_p->GetBinLowEdge(bIX+1);
    binFillPos++;
  }

  std::cout << "CHECK BINS:" << std::endl;
  for(Int_t bIX = 0; bIX < nBins+1; ++bIX){
    std::cout << bins[bIX] << ", ";
  }
  std::cout << std::endl;

  const Int_t nSystHist = syst_p.size();
  bool systToBeCombined[nSystHist];
  std::vector<std::vector<int> > postToCombine;
  std::string systLegStr[nSystHist];
  bool systToBeSkipped[nSystHist];

  for(Int_t sI = 0; sI < nSystHist; ++sI){
    systToBeCombined[sI] = false;
    systLegStr[sI] = systStr.at(sI);
    systToBeSkipped[sI] = false;
    postToCombine.push_back({});

    for(unsigned int cI = 0; cI < systToCombo.size(); ++cI){
      if(isStrSame(systStr.at(sI), systToCombo.at(cI).at(0))){
	systToBeCombined[sI] = true;
	systLegStr[sI] = systToCombo.at(cI).at(systToCombo.at(cI).size()-1);

	for(Int_t sI2 = 0; sI2 < nSystHist; ++sI2){
	  bool doCombine = false;
	  for(unsigned int skI = 1; skI < systToCombo.at(cI).size()-1; ++skI){
	    if(isStrSame(systStr.at(sI2), systToCombo.at(cI).at(skI))){
	      doCombine = true;
	      break;
	    }
	  }

	  if(doCombine) postToCombine.at(postToCombine.size()-1).push_back(sI2);
	}

	break;
      }
      else{
	bool doBreak = false;
	for(unsigned int skI = 1; skI < systToCombo.at(cI).size()-1; ++skI){
	  if(isStrSame(systStr.at(sI), systToCombo.at(cI).at(skI))){
	    systToBeSkipped[sI] = true;
	    doBreak = true;
	    break;
	  }
	}
	if(doBreak) break;
      }
    }
  }

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

      syst_p.at(sI)->SetBinContent(bIX+1, binVal);
      syst_p.at(sI)->SetBinError(bIX+1, binErr);
      binFillPos++;
    }

    delete systHist_p[sI];
  }

  for(unsigned int sI = 0; sI < syst_p.size(); ++sI){
    if(!systToBeCombined[sI]) continue;

    for(unsigned int cI = 0; cI < postToCombine.at(sI).size(); ++cI){
      std::cout << "Combining " << systStr.at(sI) << " and " << systStr.at(postToCombine.at(sI).at(cI)) << std::endl;

      for(Int_t bIX = 0; bIX < syst_p.at(sI)->GetNbinsX(); ++bIX){
	Double_t binCenter = syst_p.at(sI)->GetBinCenter(bIX+1);
	if(binCenter > maxXVal || binCenter < minXVal) continue;

	Double_t binC1 = syst_p.at(sI)->GetBinContent(bIX+1);
	Double_t binC2 = syst_p.at(postToCombine.at(sI).at(cI))->GetBinContent(bIX+1);
	
	if(binC2 > binC1){
	  syst_p.at(sI)->SetBinContent(bIX+1, binC2);
	  syst_p.at(sI)->SetBinError(bIX+1, syst_p.at(sI)->GetBinError(bIX+1)*binC2/binC1);
	}
      }
    }
  }

  for(Int_t bIX = 0; bIX < nominal_p->GetNbinsX(); ++bIX){
    Double_t binCenter = nominal_p->GetBinCenter(bIX+1);
    if(binCenter > maxXVal || binCenter < minXVal){
      for(unsigned int sI = 0; sI < syst_p.size(); ++sI){
	syst_p.at(sI)->SetBinContent(bIX+1, 0.0);
	syst_p.at(sI)->SetBinError(bIX+1, 0.0);
      }
      continue;
    }

    Double_t tempVal = 0.0;
    for(unsigned int sI = 0; sI < syst_p.size(); ++sI){
      Double_t binVal = syst_p.at(sI)->GetBinContent(bIX+1)/nominal_p->GetBinContent(bIX+1);
      Double_t tempRelErr = syst_p.at(sI)->GetBinError(bIX+1)/nominal_p->GetBinContent(bIX+1);
      
      if(!systToBeSkipped[sI]) tempVal = TMath::Sqrt(tempVal*tempVal + syst_p.at(sI)->GetBinContent(bIX+1)*syst_p.at(sI)->GetBinContent(bIX+1));

      syst_p.at(sI)->SetBinContent(bIX+1, binVal);
      syst_p.at(sI)->SetBinError(bIX+1, tempRelErr);
      if(syst_p.at(sI)->GetBinContent(bIX+1) <= 0){
	syst_p.at(sI)->SetBinContent(bIX+1, smallDelta);
	syst_p.at(sI)->SetBinError(bIX+1, smallDelta);
      }
    }
  
    std::cout << "Pushing back: " << tempVal << ", " << nominal_p->GetBinLowEdge(bIX+1) << "-" << nominal_p->GetBinLowEdge(bIX+2) << std::endl;
    systVect.push_back(tempVal);
    nomValVect.push_back(nominal_p->GetBinContent(bIX+1));
  }

  double maxVal = -1;

  for(unsigned int sI = 0; sI < syst_p.size(); ++sI){
    if(systToBeSkipped[sI]) continue;
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
    
    TLatex* label_p = new TLatex();
    label_p->SetTextFont(43);
    label_p->SetTextSize(16);

    
    TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);
    canv_p->SetTopMargin(0.01);
    canv_p->SetRightMargin(0.001);
    canv_p->SetLeftMargin(0.12);
    canv_p->SetBottomMargin(0.12);
    
    TH1D* tempHist_p = new TH1D("tempHist_p", ";Jet p_{T} (GeV);Fractional Systetmatics", nBins, bins);
    tempHist_p->GetXaxis()->SetTitleOffset(1.6);
    tempHist_p->SetMarkerStyle(1);
    tempHist_p->SetMarkerSize(0.001);
    tempHist_p->SetMarkerColor(1);
    tempHist_p->SetLineColor(1);

    canv_p->cd();
    centerTitles(tempHist_p);
    tempHist_p->SetMaximum(1.5*max[nI]);
    tempHist_p->DrawCopy("HIST E1 P");
    
    gStyle->SetOptStat(0);
    gPad->SetLogx();
  
    double systMaxVal = -1;
    unsigned int systMaxPos = -1;
    double systSubMaxVal = -1;
    unsigned int systSubMaxPos = -1;

    for(unsigned int sI = 0; sI < syst_p.size(); ++sI){
      if(systToBeSkipped[sI]) continue;

      if(systMaxVal < syst_p.at(sI)->GetMaximum()){
	systSubMaxVal = systMaxVal;
	systSubMaxPos = systMaxPos;

	systMaxVal = syst_p.at(sI)->GetMaximum();
	systMaxPos = sI;
      }
      else if(systSubMaxVal < syst_p.at(sI)->GetMaximum()){
	systSubMaxVal = syst_p.at(sI)->GetMaximum();
	systSubMaxPos = sI;
      }
    }
  
    for(unsigned int sI = 0; sI < syst_p.size(); ++sI){
      if(systToBeSkipped[sI]) continue;

      Int_t binFillPos = 0;
      for(Int_t bIX = 0; bIX < syst_p.at(sI)->GetNbinsX(); ++bIX){
	Double_t binCenter = syst_p.at(sI)->GetBinCenter(bIX+1);
	if(binCenter > maxXVal || binCenter < minXVal) continue;

	Double_t binVal = tempHist_p->GetBinContent(binFillPos+1)*tempHist_p->GetBinContent(binFillPos+1);
	binVal += syst_p.at(sI)->GetBinContent(bIX+1)*syst_p.at(sI)->GetBinContent(bIX+1);
	binVal = TMath::Sqrt(binVal);

	tempHist_p->SetBinContent(binFillPos+1, binVal);
	tempHist_p->SetBinError(binFillPos+1, smallDelta);

	binFillPos++;
      }

      syst_p.at(sI)->SetMarkerColor(colors[sI%nColor]);
      syst_p.at(sI)->SetMarkerStyle(styles[sI%nStyle]);
      syst_p.at(sI)->SetLineColor(colors[sI%nColor]);
      syst_p.at(sI)->SetMarkerSize(0.8);
   
      syst_p.at(sI)->DrawCopy("HIST E1 P SAME");
      syst_p.at(sI)->DrawCopy("HIST E1 SAME");
      
      std::string legStr = systLegStr[sI];
      if(systMaxPos == sI) legStr = legStr + " (Dominant)";
      else if(systSubMaxPos == sI) legStr = legStr + " (Sub-dominant)";
      leg_p->AddEntry(syst_p.at(sI), legStr.c_str(), "P L");
    }

    tempHist_p->DrawCopy("HIST E1 SAME");
  
    if(nI == 0){
      std::cout << "Quick closure check: " << std::endl;
      for(Int_t bIX = 0; bIX < nBins; ++bIX){
	std::string binLowStr = prettyString(tempHist_p->GetBinLowEdge(bIX+1), 1, false);
	std::string binHiStr = prettyString(tempHist_p->GetBinLowEdge(bIX+2), 1, false);
	std::cout << " " << bIX << ", " << binLowStr << "-" << binHiStr << ": " << tempHist_p->GetBinContent(bIX+1) << ", " << systVect.at(bIX)/nomValVect.at(bIX) << std::endl;
      }
    }
    
    drawWhiteBox(900, 1100, -max[nI]*9./10., -max[nI]/100.);

    double yPos = 0 - max[nI]/20.;

    label_p->DrawLatex(200, yPos, "200");
    label_p->DrawLatex(400, yPos, "400");
    label_p->DrawLatex(600, yPos, "600");
    label_p->DrawLatex(800, yPos, "800");
	
    leg_p->Draw("SAME");

    std::string smoothStr = "Unsmoothed";
    if(doSmoothing) smoothStr = "Smoothed";

    std::string maxValStr = "NormalMax";
    if(nI != 0) maxValStr = "ZOOM";
    
    Int_t nLabelBins = 20;
    Double_t labelBins[nLabelBins+1];
    getLinBins(0, max[nI], nLabelBins, labelBins);

    label_p->DrawLatex(500, labelBins[19], maxValStr.c_str());
    label_p->DrawLatex(500, labelBins[17], smoothStr.c_str());

    gPad->RedrawAxis();
    gPad->SetTicks(1,2);

    std::string nominalName = nominal_p->GetName();
    std::string ppCentTag = "PP";
    if(nominalName.find("Cent") != std::string::npos){
      ppCentTag = nominalName.substr(nominalName.find("Cent"), nominalName.size());
      ppCentTag.replace(0, 4, "");
      ppCentTag.replace(ppCentTag.find("to"), 2, "-");
      ppCentTag.replace(ppCentTag.find("_"), ppCentTag.size(), "");
      ppCentTag = "PbPb " + ppCentTag + "%";
    }
    label_p->DrawLatex(500, labelBins[15], ppCentTag.c_str());

    const std::string dirName = "pdfDir/" + dateStr;
    checkMakeDir(dirName);
    const std::string saveName = "systErr_" + nominalName + "_" + maxValStr + "_" + smoothStr + "_" + dateStr +  ".pdf";

    quietSaveAs(canv_p, (dirName + "/" + saveName));
    plotNames->push_back(saveName);
  
    delete tempHist_p;
    delete canv_p;
    delete leg_p;
    delete label_p;
  }

  return systVect;
}


void drawSyst(TCanvas* canv_p, TH1D* nominal_p, std::vector<double> syst_, Double_t minXVal, Double_t maxXVal)
{
  canv_p->cd();

  TBox* tempBox_p = new TBox();
  tempBox_p->SetFillColorAlpha(nominal_p->GetMarkerColor(), .25);

  Int_t binFillPos = 0;
  for(Int_t bIX = 0; bIX < nominal_p->GetNbinsX(); ++bIX){
    Double_t binCenter = nominal_p->GetBinCenter(bIX+1);
    Double_t binLowEdge = nominal_p->GetBinLowEdge(bIX+1);
    Double_t binHiEdge = nominal_p->GetBinLowEdge(bIX+2);
    Double_t binContent = nominal_p->GetBinContent(bIX+1);
    
    if(binCenter < minXVal) continue;
    if(binCenter > maxXVal) continue;
    
    tempBox_p->DrawBox(binLowEdge, binContent - syst_.at(binFillPos), binHiEdge, binContent + syst_.at(binFillPos));
    ++binFillPos;
  }

  return;
}


void doRatioPlotting(std::vector<TH1D*> inHist_p, int divPos, const std::string xStr, const std::string yStr, const int xMin, const int xMax, const std::string labelStr, std::string saveName)
{
  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  TCanvas* canv_p = new TCanvas("convCanv_p", "", 450, 450);
  canv_p->SetTopMargin(0.01);
  canv_p->SetRightMargin(0.01);
  canv_p->SetLeftMargin(0.14);
  canv_p->SetBottomMargin(0.14);
  canv_p->cd();
	  
  TLatex* label_p = new TLatex();
  label_p->SetTextFont(43);
  label_p->SetTextSize(16);
	  
  TH1D* tempHist_p = new TH1D("tempHist_p", (";" + xStr + ";" + yStr).c_str(), 10, xMin, xMax);

  centerTitles(tempHist_p);
  tempHist_p->Draw("HIST");

  Double_t maxVal = -1;
  Double_t minVal = 1000000;
	  
  const Int_t nTempStyle = 3;
  const Int_t tempStyle[nTempStyle] = {20, 21, 34};

  kirchnerPalette kPalette;  
  const Int_t nTempColor = 4;
  const Int_t tempColor[nTempColor] = {kPalette.getColor(0), kPalette.getColor(1), kPalette.getColor(2), kPalette.getColor(3)};

  const Int_t nClones = inHist_p.size();
  std::vector<TH1D*> clones_p;
  clones_p.reserve(nClones);
  for(Int_t bI = 0; bI < nClones; ++bI){
    clones_p.push_back((TH1D*)inHist_p.at(bI)->Clone(("clone_" + std::to_string(bI)).c_str()));
  }
	
  for(unsigned int bI = 0; bI < clones_p.size(); ++bI){
    clones_p.at(bI)->SetMarkerStyle(tempStyle[bI%nTempStyle]);
    clones_p.at(bI)->SetMarkerSize(0.8);
    clones_p.at(bI)->SetMarkerColor(tempColor[bI%nTempColor]);
    clones_p.at(bI)->SetLineColor(tempColor[bI%nTempColor]);
    
    if(bI == (unsigned int)divPos) continue;
    
    clones_p.at(bI)->Divide(clones_p.at(divPos));
    
    for(Int_t bIX = 0; bIX < clones_p.at(bI)->GetNbinsX(); ++bIX){
      Double_t binCenter = clones_p.at(bI)->GetBinCenter(bIX+1);
      
      clones_p.at(bI)->SetBinError(bIX+1, 0.0);
      if(binCenter >= xMin) continue;
      if(binCenter <= xMax) continue;
      
      clones_p.at(bI)->SetBinContent(bIX+1, 0.0);
    }
    
    Double_t tempMaxVal = clones_p.at(bI)->GetMaximum();
    Double_t tempMinVal = getMinGTZero(clones_p.at(bI));
    if(tempMaxVal > maxVal) maxVal = tempMaxVal;
    if(tempMinVal < minVal) minVal = tempMinVal;
  }
  
  clones_p.at(divPos)->Divide(clones_p.at(divPos));
  
  for(Int_t bIX = 0; bIX < clones_p.at(divPos)->GetNbinsX(); ++bIX){
    Double_t binCenter = clones_p.at(divPos)->GetBinCenter(bIX+1);
    
    clones_p.at(divPos)->SetBinError(bIX+1, 0.0);
    if(binCenter >= 200.) continue;
    if(binCenter <= 1000.) continue;
    
    clones_p.at(divPos)->SetBinContent(bIX+1, 0.0);
  }

  const std::string maxValStr = "Max: " + prettyString(maxVal, 3, false); 
  const std::string minValStr = "Min: " + prettyString(minVal, 3, false);
 
  Double_t interval = maxVal - minVal;
  maxVal += interval/10.;
  minVal = TMath::Max(0., minVal - interval/10.);

  
  if(maxVal < 1.0001 && minVal > .9999){
    maxVal = 1.0001;
    minVal = 0.9999;
  }
  

  tempHist_p->SetMaximum(maxVal);
  tempHist_p->SetMinimum(minVal);
  
  for(Int_t bI = 0; bI < nClones; ++bI){
    clones_p.at(bI)->DrawCopy("HIST E1 P SAME");
  }
  
  const Int_t nTempBins = 20;
  Double_t tempBins[nTempBins+1];

  getLinBins(minVal, maxVal, nTempBins, tempBins);
  label_p->DrawLatex(250, tempBins[19], labelStr.c_str());
  label_p->DrawLatex(250, tempBins[18], maxValStr.c_str());
  label_p->DrawLatex(250, tempBins[17], minValStr.c_str());

  gStyle->SetOptStat(0);
  gPad->SetLogx();
	  
  saveName = "pdfDir/" + dateStr + "/" + saveName;
  quietSaveAs(canv_p, saveName);
  
  delete tempHist_p;
  delete canv_p;
  delete label_p;

  return;
}


int plotUnfoldedSpectra(const std::string inFileNamePP, const std::string inFileNamePbPb)
{
  std::vector<std::string> slideTitles;
  std::vector<std::vector<std::string > > pdfPerSlide;
  std::vector<std::string> slideTitlesMain;
  std::vector<std::vector<std::string > > pdfPerSlideMain;

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

  const Int_t nSystOrig = cutPropPbPb.GetNSyst() - 1; // Minus 1 because removing priorflat
  //  const Int_t nSyst = TMath::Min(minForForLoop, nSystOrig + nAdditionalSyst);
  const Int_t nSyst = nSystOrig + nAdditionalSyst; 
  std::vector<std::string> systStr = cutPropPbPb.GetSystStr();
  Int_t jerDataPos = -1;

  for(unsigned int sI = 0; sI < systStr.size(); ++sI){
    if(isStrSame(systStr.at(sI), "PriorFlat")){
      systStr.erase(systStr.begin()+sI);
      break;
    }
  }

  for(Int_t sI = 0; sI < nAdditionalSyst; ++sI){
    systStr.push_back(additionalSyst[sI]);
  }

  std::cout << "Processing full systematics: " << std::endl;
  
  for(Int_t sI = 0; sI < nSyst; ++sI){
    std::cout << " " << sI << "/" << nSyst << ": " << systStr.at(sI) << std::endl;
    if(isStrSame(systStr.at(sI), "JERData")) jerDataPos = sI - 1;
  }

  if(jerDataPos >= 0) std::cout << "JERDATAPOS: " << jerDataPos << std::endl;
  else std::cout << "WARNING: JERDATAPOS NOT FOUND" << std::endl;

  const Int_t nBayesCap = 10000;
  const Int_t nBayes = TMath::Min(nBayesCap, cutPropPbPb.GetNBayes());
  const Int_t nBigBayesSymm = cutPropPbPb.GetNBigBayesSymm();
  std::vector<int> bayesVal = cutPropPbPb.GetBayesVal();
  const Int_t nSuperBayes = cutPropPbPb.GetNSuperBayes();

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
  outFileName = "output/" + outFileName + "_" + outFileName2 + "_NSuperBayes" + std::to_string(nSuperBayes) + "_PlotUnfoldedSpectra_" + dateStr + ".root";


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
  Int_t bayesPosPbPb[nJtPbPb][nCentBins][nID][nResponseMod][nJtAbsEtaBins][nSyst];
  

  for(Int_t tI = 0; tI < nJtPP; ++tI){
    std::string tempStr = jetPPList.at(tI);
    if(tempStr.find("/") != std::string::npos) tempStr.replace(tempStr.find("/"), tempStr.size() - tempStr.find("/"), "");
    dirPP_p[tI] = (TDirectory*)outFile_p->mkdir(tempStr.c_str());
    dirPP_p[tI]->cd();
  }

  TH1D* jtPtUnfolded_RecoTrunc_PP_h[nJtPP][nID][nResponseMod][nJtAbsEtaBins][nSyst][nBayes];
  Int_t bayesPosPP[nJtPP][nID][nResponseMod][nJtAbsEtaBins][nSyst];


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
      if(systPos < 0 && name.find("PriorFlat") == std::string::npos) systPos = 0;

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
      if(systPos < 0){
	if(name.find("PriorFlat") == std::string::npos) std::cout << "Warning: Cannot find systPos for \'" << name << "\'. Code will likely break" << std::endl;
	else continue;
      }
      if(bayesPos < 0) std::cout << "Warning: Cannot find bayesPos for \'" << name << "\'. Code will likely break" << std::endl;
      if(centPos < 0) std::cout << "Warning: Cannot find centPos for \'" << name << "\'. Code will likely break" << std::endl;

      jtPtUnfolded_RecoTrunc_PbPb_h[tI][centPos][idPos][modPos][absEtaPos][systPos][bayesPos] = (TH1D*)key->ReadObj();   
    }

    //Building best bayes
    for(Int_t cI = 0; cI < nCentBins; ++cI){
      const std::string centStr = "Cent" + std::to_string(centBinsLow.at(cI)) + "to" + std::to_string(centBinsHi.at(cI));
      
      for(unsigned int idI = 0; idI < idStr.size(); ++idI){
	for(Int_t mI = 0; mI < nResponseMod; ++mI){
	  const std::string resStr = "ResponseMod" + prettyString(responseMod.at(mI), 2, true);
	  
	  for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
	    const std::string jtAbsEtaStr = "AbsEta" + prettyString(jtAbsEtaBinsLow.at(aI), 1, true) + "to" + prettyString(jtAbsEtaBinsHi.at(aI), 1, true);	
	    
	    for(Int_t sI = 0; sI < nSystOrig; ++sI){
	      const std::string tempHistTag = jetPbPbList.at(tI) + "_" + centStr + "_" + idStr.at(idI) + "_" + resStr + "_" + jtAbsEtaStr + "_" + systStr.at(sI);
	      if(histTagMapPbPb.count(tempHistTag) == 0){
		std::cout << "WARNING: tag \'" << tempHistTag << "\' not found. Setting to -1" << std::endl;
		bayesPosPbPb[tI][cI][idI][mI][aI][sI] = -1;
	      }
	      else{
		int tempSet = histTagMapPbPb[tempHistTag];
		std::cout << "tag \'" << tempHistTag << "\' found. Setting to " << tempSet << std::endl;
		bayesPosPbPb[tI][cI][idI][mI][aI][sI] = tempSet;
	      }	      
	    }
	  }
	}
      }
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

	      bayesPosPbPb[tI][cI][idI][mI][aI][sI] = bayesPosPbPb[tI][cI][idI][mI][aI][0];

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
      if(systPos < 0 && name.find("PriorFlat") == std::string::npos) systPos = 0;

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
      if(systPos < 0){
	if(name.find("PriorFlat") == std::string::npos) std::cout << "Warning: Cannot find systPos for \'" << name << "\'. Code will likely break" << std::endl;
	else continue;
      }
      if(bayesPos < 0) std::cout << "Warning: Cannot find bayesPos for \'" << name << "\'. Code will likely break" << std::endl;

      jtPtUnfolded_RecoTrunc_PP_h[tI][idPos][modPos][absEtaPos][systPos][bayesPos] = (TH1D*)key->ReadObj();   
    }

    //Building best bayes
    for(unsigned int idI = 0; idI < idStr.size(); ++idI){
      for(Int_t mI = 0; mI < nResponseMod; ++mI){
	const std::string resStr = "ResponseMod" + prettyString(responseMod.at(mI), 2, true);
	
	for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
	  const std::string jtAbsEtaStr = "AbsEta" + prettyString(jtAbsEtaBinsLow.at(aI), 1, true) + "to" + prettyString(jtAbsEtaBinsHi.at(aI), 1, true);	
	  
	  for(Int_t sI = 0; sI < nSystOrig; ++sI){
	    const std::string tempHistTag = jetPPList.at(tI) + "_PP_" + idStr.at(idI) + "_" + resStr + "_" + jtAbsEtaStr + "_" + systStr.at(sI);
	    if(histTagMapPP.count(tempHistTag) == 0){
	      std::cout << "WARNING: tag \'" << tempHistTag << "\' not found. Setting to -1" << std::endl;
	      bayesPosPP[tI][idI][mI][aI][sI] = -1;
	    }
	    else{
	      int tempSet = histTagMapPP[tempHistTag];
	      std::cout << "tag \'" << tempHistTag << "\' found. Setting to " << tempSet << std::endl;
	      bayesPosPP[tI][idI][mI][aI][sI] = tempSet;
	    }
	  }
	}
      }
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

	    bayesPosPP[tI][idI][mI][aI][sI] = bayesPosPP[tI][idI][mI][aI][0];
	    
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
  const Int_t plotBayesVal = 100;

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

  std::vector<std::vector<std::string> > systToCombo = {{"JECUpMC", "JECDownMC", "JECMC"}, {"JECUpData", "JECDownData", "JECData"}, {"JECUpUE", "JECDownUE", "JECUE"}, {"PriorUp1PowerPthat", "PriorDown1PowerPthat", "Prior1PowerPthat"}, {"LumiUp", "LumiDown", "Lumi"}, {"TAAUp", "TAADown", "TAA"}};

  const std::string idNameStr = idStr.at(idPos);
  const std::string plotBayesStr = "Bayes" + std::to_string(plotBayesVal);

  const Int_t nSystSmooth = 2;
  const std::string systSmooth[nSystSmooth] = {"Unsmoothed", "Smoothed"};
  const bool systSmoothBool[nSystSmooth] = {false, true};

  //Check convergence bands...
  std::cout << "Do convergence bands" << std::endl;
  for(Int_t tI = 0; tI < nJtPbPb; ++tI){
    const Int_t rVal = getRVal(jetPbPbList.at(tI));
    //    const Int_t ppPos = jetPPMatchedToPbPb.at(tI);
    const std::string rValStr = std::to_string(rVal);

    std::string rStr = "R=";
    if(getRVal(jetPbPbList.at(tI)) < 10) rStr = rStr + "0." + std::to_string(getRVal(jetPbPbList.at(tI)));
    else rStr = rStr + "1.0";

    const std::string jetStrPbPb = jetPbPbList.at(tI).substr(0, jetPbPbList.at(tI).find("JetAna"));;

    //	      jtPtUnfolded_RecoTrunc_PP_h[tI][idI][mI][aI][sI][bI] = (TH1D*)jtPtUnfolded_RecoTrunc_PP_h[tI][idI][mI][aI][0][bI]->Clone(newName.c_str());

    for(Int_t mI = 0; mI < nResponseMod; ++mI){
      for(Int_t cI = 0; cI < nCentBins; ++ cI){
	const std::string centStr = "Cent" + std::to_string(centBinsLow.at(cI)) + "to" + std::to_string(centBinsHi.at(cI));
	const std::string centStr2 = std::to_string(centBinsLow.at(cI)) + "-" + std::to_string(centBinsHi.at(cI)) + "%";

	
	slideTitles.push_back("Conv. " + jetStrPbPb + " (ModRes " + prettyString(responseMod[mI],2,false) + ", " + centStr2 + ")");
	pdfPerSlide.push_back({});

	for(Int_t sI = 0; sI < nSyst; ++sI){	  
	  const Int_t totBigBayes = 1 + 2*nBigBayesSymm;
	  std::vector<TH1D*> bigBayesClones_p;
	  bigBayesClones_p.reserve(totBigBayes);
	  for(Int_t bI = 0; bI < totBigBayes; ++bI){
	    Int_t binPos = nBayes - 1 -2*nBigBayesSymm + bI;
	    bigBayesClones_p.push_back((TH1D*)jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idPos][mI][absEtaPos][sI][binPos]->Clone(("bigBayesClone_" + std::to_string(bI)).c_str()));
	  }      

	  std::string saveName = "convBayes_NBigBayesSymm" + std::to_string(nBigBayesSymm) + "_" + jetPbPbList.at(tI) + "_R" + rValStr + "_" + centStr + "_" + idNameStr + "_ResponseMod" + prettyString(responseMod[mI], 2, true) + "_" + plotAbsEtaStr + "_" + systStr.at(sI) + "_" + dateStr + ".pdf";
	  pdfPerSlide.at(pdfPerSlide.size()-1).push_back(saveName);

	  doRatioPlotting(bigBayesClones_p, nBigBayesSymm, "Jet p_{T} (GeV/c)", "Ratio (Bayes Iter 100 #pm" + std::to_string(nBigBayesSymm) + ")", 200., 1000., systStr.at(sI), saveName);
	  for(unsigned int ucI = 0; ucI < bigBayesClones_p.size(); ++ucI){
	    delete bigBayesClones_p.at(ucI);
	    bigBayesClones_p.at(ucI) = NULL;
	  }
	}

	slideTitles.push_back("BestComp. " + jetStrPbPb + " (ModRes " + prettyString(responseMod[mI],2,false) + ", " + centStr2 + ")");
	pdfPerSlide.push_back({});

	for(Int_t sI = 0; sI < nSyst; ++sI){	  
	  std::vector<TH1D*> clones_p;

	  Int_t binPos = nBayes - 1 -nBigBayesSymm;	  
	  clones_p.push_back((TH1D*)jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idPos][mI][absEtaPos][sI][binPos]->Clone(("clone_" + std::to_string(binPos)).c_str()));
	  
	  if(bayesPosPbPb[tI][cI][idPos][mI][absEtaPos][sI] < 0) bayesPosPbPb[tI][cI][idPos][mI][absEtaPos][sI] = 25;

	  Int_t tempBayesPos = bayesPosPbPb[tI][cI][idPos][mI][absEtaPos][sI];

	  clones_p.push_back((TH1D*)jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idPos][mI][absEtaPos][sI][tempBayesPos]->Clone(("clone_" + std::to_string(tempBayesPos)).c_str()));

	  std::string saveName = "convBestBayes_NBigBayesSymm" + std::to_string(nBigBayesSymm) + "_" + jetPbPbList.at(tI) + "_R" + rValStr + "_" + centStr + "_" + idNameStr + "_ResponseMod" + prettyString(responseMod[mI], 2, true) + "_" + plotAbsEtaStr + "_" + systStr.at(sI) + "_" + dateStr + ".pdf";
	  pdfPerSlide.at(pdfPerSlide.size()-1).push_back(saveName);

	  std::string labelStr = systStr.at(sI) + ", Bayes" + std::to_string(tempBayesPos);
	  doRatioPlotting(clones_p, 0, "Jet p_{T} (GeV/c)", "Ratio (Bayes Iter 100 w/ " + std::to_string(tempBayesPos) + ")", 200., 1000., labelStr, saveName);

	  for(unsigned int ucI = 0; ucI < clones_p.size(); ++ucI){
	    delete clones_p.at(ucI);
	    clones_p.at(ucI) = NULL;
	  }

	}
      }
    }
  }


  //Actual spectra plots...
  for(Int_t tI = 0; tI < nJtPbPb; ++tI){
    const Int_t rVal = getRVal(jetPbPbList.at(tI));
    const Int_t ppPos = jetPPMatchedToPbPb.at(tI);
    const std::string rValStr = std::to_string(rVal);

    std::string rStr = "R=";
    if(getRVal(jetPbPbList.at(tI)) < 10) rStr = rStr + "0." + std::to_string(getRVal(jetPbPbList.at(tI)));
    else rStr = rStr + "1.0";

    Double_t xMinVal = 200;
    Double_t xPointMinVal = xMinVal;
    if(rVal >= 8) xPointMinVal = 300;
    Double_t xMaxVal = 1200;
    Double_t xPointMaxVal = 1000;

    TLegend* legRes_p = new TLegend(0.7, 0.65, 0.95, 0.88);
    legRes_p->SetBorderSize(0);
    legRes_p->SetFillColor(0);
    legRes_p->SetFillStyle(0);
    legRes_p->SetTextFont(43);
    legRes_p->SetTextSize(16);

    const Int_t nPad = 2;
    TCanvas* resModCompCanv_p[nSystSmooth];
    TPad* resModCompPad_p[nSystSmooth][nPad];
    TH1D* resModHist_p[nSystSmooth];
    TH1D* resModHistRat_p[nSystSmooth];

    const Double_t yPadFrac = 0.35;
    const Double_t yPadLow[nPad] = {yPadFrac, 0.0};
    const Double_t yPadHi[nPad] = {1.0, yPadFrac};
    const Double_t padLeftMarg = 0.14;

    for(Int_t usI = 0; usI < nSystSmooth; ++usI){
      resModCompCanv_p[usI] = new TCanvas(("resModCompCanv_" + systSmooth[usI]).c_str(), "", 450, 450);
      resModCompCanv_p[usI]->SetTopMargin(0.001);
      resModCompCanv_p[usI]->SetRightMargin(0.001);
      resModCompCanv_p[usI]->SetLeftMargin(0.001);
      resModCompCanv_p[usI]->SetBottomMargin(0.001);    

      for(Int_t pI = 0; pI < nPad; ++pI){
	resModCompCanv_p[usI]->cd();
	resModCompPad_p[usI][pI] = new TPad(("resModCompPad_" + systSmooth[usI] + "_" + std::to_string(pI)).c_str(), "", 0.0, yPadLow[pI], 1.0, yPadHi[pI]);
	resModCompPad_p[usI][pI]->SetRightMargin(0.001);
	resModCompPad_p[usI][pI]->SetLeftMargin(padLeftMarg);

	if(pI == 0) resModCompPad_p[usI][pI]->SetTopMargin(0.01);
	else resModCompPad_p[usI][pI]->SetTopMargin(0.001);

	if(pI == 0) resModCompPad_p[usI][pI]->SetBottomMargin(0.001);
	else resModCompPad_p[usI][pI]->SetBottomMargin(padLeftMarg/yPadFrac);

	resModCompPad_p[usI][pI]->Draw("SAME");
      }

      resModHist_p[usI] = new TH1D(("resModHist_" + systSmooth[usI]).c_str(), ";Jet p_{T} (GeV/c);Fractional Systematic", 10, xPointMinVal, xPointMaxVal);
      centerTitles(resModHist_p[usI]);

      resModHist_p[usI]->SetMaximum(0.0);
      resModHist_p[usI]->SetMinimum(-1);
      resModCompPad_p[usI][0]->cd();
      resModHist_p[usI]->GetYaxis()->SetTitleFont(43);
      resModHist_p[usI]->GetYaxis()->SetTitleSize(16);
      resModHist_p[usI]->GetYaxis()->SetLabelFont(43);
      resModHist_p[usI]->GetYaxis()->SetLabelSize(10);
      resModHist_p[usI]->Draw("HIST");


      resModHistRat_p[usI] = new TH1D(("resModHistRat_" + systSmooth[usI]).c_str(), ";Jet p_{T} (GeV/c);Old/New", 10, xPointMinVal, xPointMaxVal);
      centerTitles(resModHistRat_p[usI]);

      resModHistRat_p[usI]->SetMaximum(1.1);
      resModHistRat_p[usI]->SetMinimum(0.4);
      resModCompPad_p[usI][1]->cd();

      resModHistRat_p[usI]->GetYaxis()->SetTitleFont(43);
      resModHistRat_p[usI]->GetYaxis()->SetTitleSize(16);
      resModHistRat_p[usI]->GetYaxis()->SetLabelFont(43);
      resModHistRat_p[usI]->GetYaxis()->SetLabelSize(10);

      resModHistRat_p[usI]->GetXaxis()->SetTitleFont(43);
      resModHistRat_p[usI]->GetXaxis()->SetTitleSize(16);
      resModHistRat_p[usI]->GetXaxis()->SetLabelFont(43);
      resModHistRat_p[usI]->GetXaxis()->SetLabelSize(10);

      resModHistRat_p[usI]->GetXaxis()->SetTitleOffset(2.0);

      resModHistRat_p[usI]->Draw("HIST");
    }
  
    std::vector<TH1D*> prevJERData;

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

	Int_t slideSubVal = 0;
	if(usI == 0) pdfPerSlide.push_back({});
	else slideSubVal = nCentBins+1;

	std::cout << __LINE__ << std::endl;
	std::vector<std::string> tempSystStr(systStr.begin()+1, systStr.end());
	std::vector<double> systValVectPP = getSyst(jtPtUnfolded_RecoTrunc_PP_h[ppPos][idPos][mI][absEtaPos][0][bayesPos], systHistVectPP, tempSystStr, xPointMinVal, xPointMaxVal, &(pdfPerSlide.at(pdfPerSlide.size()- 1 - slideSubVal)), systSmoothBool[usI], systToCombo);
	std::cout << __LINE__ << std::endl;
    
	resModCompCanv_p[usI]->cd();
	resModCompPad_p[usI][0]->cd();
	Int_t resStyle = -1;
	Int_t resCol = kPalette.getColor(getColorPosFromCent("", true));
	std::string drawStr = "HIST P SAME";
	if(mI == 0) resStyle = getStyleFromCent("", true);
	else resStyle = getAltStyleFromCent("", true);

	systHistVectPP.at(jerDataPos)->SetMarkerSize(0.8);
	systHistVectPP.at(jerDataPos)->SetMarkerStyle(resStyle);
	systHistVectPP.at(jerDataPos)->SetMarkerColor(resCol);
	systHistVectPP.at(jerDataPos)->SetLineColor(resCol);

	systHistVectPP.at(jerDataPos)->DrawCopy(drawStr.c_str());

	if(systHistVectPP.at(jerDataPos)->GetMaximum() > resModHist_p[usI]->GetMaximum()) resModHist_p[usI]->SetMaximum(systHistVectPP.at(jerDataPos)->GetMaximum());

	if(mI == 0) prevJERData.push_back((TH1D*)systHistVectPP.at(jerDataPos)->Clone((std::string(systHistVectPP.at(jerDataPos)->GetName()) + "_CLONE").c_str()));
	else{
	  if(usI == 0) systHistVectPP.at(jerDataPos)->Divide(prevJERData.at(0));
	  else if(usI == 1) systHistVectPP.at(jerDataPos)->Divide(prevJERData.at(1 + nCentBins));

	  resModCompCanv_p[usI]->cd();
	  resModCompPad_p[usI][1]->cd();

	  systHistVectPP.at(jerDataPos)->DrawCopy("HIST P SAME");	  
	}
   
	std::cout << __LINE__ << std::endl;

	if(usI == 0) slideTitles.push_back("PP Systematic (" + jetPPList.at(ppPos) + ", " + responseStr +  ")");

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
	
	  Int_t centSlideSubVal = 0;
	  if(usI == 0) pdfPerSlide.push_back({});
	  else centSlideSubVal = slideSubVal - 1 - cI;

	  std::cout << __LINE__ << std::endl;

	  std::vector<double> systValVectPbPb = getSyst(jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idPos][mI][absEtaPos][0][bayesPos], systHistVectPbPb, tempSystStr, xPointMinVal, xPointMaxVal, &(pdfPerSlide.at(pdfPerSlide.size() - 1 - centSlideSubVal)), systSmoothBool[usI], systToCombo);
	  std::cout << __LINE__ << std::endl;
	  if(usI == 0) slideTitles.push_back(centStr2 + " Systematic (" + jetPbPbList.at(tI) + ", " + responseStr +  ")");

	  resModCompCanv_p[usI]->cd();
	  resModCompPad_p[usI][0]->cd();
	  resStyle = -1;
	  resCol = kPalette.getColor(getColorPosFromCent(centStr, false));
	  drawStr = "HIST P SAME";
	  if(mI == 0) resStyle = getStyleFromCent(centStr, false);
	  else resStyle = getAltStyleFromCent(centStr, false);
	  
	  systHistVectPbPb.at(jerDataPos)->SetMarkerSize(0.8);
	  systHistVectPbPb.at(jerDataPos)->SetMarkerStyle(resStyle);
	  systHistVectPbPb.at(jerDataPos)->SetMarkerColor(resCol);
	  systHistVectPbPb.at(jerDataPos)->SetLineColor(resCol);

	  std::cout << __LINE__ << std::endl;
	  
	  systHistVectPbPb.at(jerDataPos)->DrawCopy(drawStr.c_str());

	  if(systHistVectPbPb.at(jerDataPos)->GetMaximum() > resModHist_p[usI]->GetMaximum()) resModHist_p[usI]->SetMaximum(systHistVectPbPb.at(jerDataPos)->GetMaximum());

	  if(mI == 0) prevJERData.push_back((TH1D*)systHistVectPbPb.at(jerDataPos)->Clone((std::string(systHistVectPbPb.at(jerDataPos)->GetName()) + "_CLONE").c_str()));	  
	  else{
	    if(usI == 0) systHistVectPbPb.at(jerDataPos)->Divide(prevJERData.at(0 + cI + 1));
	    else if(usI == 1) systHistVectPbPb.at(jerDataPos)->Divide(prevJERData.at(1 + nCentBins + cI + 1));
	    
	    resModCompCanv_p[usI]->cd();
	    resModCompPad_p[usI][1]->cd();
	    
	    systHistVectPbPb.at(jerDataPos)->DrawCopy("HIST P SAME");
	  }
	

	  std::cout << __LINE__ << std::endl;

	  for(unsigned int sI = 0; sI < systHistVectPbPb.size(); ++sI){
	    delete systHistVectPbPb.at(sI);
	  }
	  drawSyst(spectCanv_p, jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idPos][mI][absEtaPos][0][bayesPos], systValVectPbPb, xPointMinVal, xPointMaxVal);	

	  std::cout << __LINE__ << std::endl;
	}

	jtPtUnfolded_RecoTrunc_PP_h[ppPos][idPos][mI][absEtaPos][0][bayesPos]->DrawCopy("HIST E1 P SAME");

	for(Int_t cI = 0; cI < nCentBins; ++cI){
	  jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idPos][mI][absEtaPos][0][bayesPos]->DrawCopy("HIST E1 P SAME");
	}
      
	
	std::cout << __LINE__ << std::endl;

	for(Int_t cI = nCentBins-1; cI >= 0; --cI){
	  const std::string centLegStr = std::to_string(centBinsLow.at(cI)) + "-" + std::to_string(centBinsHi.at(cI)) + "% x 10^{" + std::to_string(cI+1) + "}";
	  const std::string centLegStr2 = std::to_string(centBinsLow.at(cI)) + "-" + std::to_string(centBinsHi.at(cI)) + "%";

	  jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idPos][mI][absEtaPos][0][bayesPos]->SetFillColorAlpha(jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idPos][mI][absEtaPos][0][bayesPos]->GetMarkerColor(), .25);      
	  leg_p->AddEntry(jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idPos][mI][absEtaPos][0][bayesPos], centLegStr.c_str(), "P L F");	
	  if(usI == 0 && mI == 0) legRes_p->AddEntry(jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idPos][mI][absEtaPos][0][bayesPos], centLegStr2.c_str(), "P L");	
	}
	
	jtPtUnfolded_RecoTrunc_PP_h[ppPos][idPos][mI][absEtaPos][0][bayesPos]->SetFillColorAlpha(jtPtUnfolded_RecoTrunc_PP_h[ppPos][idPos][mI][absEtaPos][0][bayesPos]->GetMarkerColor(), .25);      
	leg_p->AddEntry(jtPtUnfolded_RecoTrunc_PP_h[ppPos][idPos][mI][absEtaPos][0][bayesPos], "pp", "P L F");
	if(usI == 0 && mI == 0) legRes_p->AddEntry(jtPtUnfolded_RecoTrunc_PP_h[ppPos][idPos][mI][absEtaPos][0][bayesPos], "pp", "P L");
      
	gPad->SetLogy();
	gPad->SetLogx();
	gStyle->SetOptStat(0);

	drawWhiteBox(900, 1100, minVal/10, minVal*.9);

	std::cout << __LINE__ << std::endl;

	label_p->DrawLatex(200, yBins[0]/4., "200");
	label_p->DrawLatex(400, yBins[0]/4., "400");
	label_p->DrawLatex(600, yBins[0]/4., "600");
	label_p->DrawLatex(800, yBins[0]/4., "800");
	label_p->DrawLatex(1000, yBins[0]/4., "1000");
          
	label_p->DrawLatex(250, yBins[19], ("anti-k_{T} " + rStr).c_str());
	label_p->DrawLatex(250, yBins[18], "|#eta_{jets}| < 2");

	label_p->SetNDC();
      
	label_p->DrawLatex(0.1, 0.95, "#bf{CMS Preliminary}");
	label_p->DrawLatex(0.4, 0.95, "27.4 pb^{-1} pp + 404 #mub^{-1} PbPb (5.02 TeV)");

	leg_p->Draw("SAME");

	std::cout << __LINE__ << std::endl;

	gPad->RedrawAxis();
	gPad->SetTicks(1,2);

	std::cout << __LINE__ << std::endl;
      
	std::string saveName = "spectra_" + jetPbPbList.at(tI) + "_R" + rValStr + "_" + idNameStr + "_" + plotAbsEtaStr + "_" + plotBayesStr + "_" + responseStr + "_" + systSmooth[usI] + "_" + dateStr + ".pdf";

	std::cout << __LINE__ << std::endl;

	if(usI == 0){
	  slideTitles.push_back("Spectra (" + responseStr + ")");
	  pdfPerSlide.push_back({});

	  slideTitlesMain.push_back("Spectra (" + responseStr + ")");
	  pdfPerSlideMain.push_back({});
	}

	std::cout << __LINE__ << std::endl;

	pdfPerSlide.at(pdfPerSlide.size() - 1).push_back(saveName);
	std::cout << __LINE__ << std::endl;

	pdfPerSlideMain.at(pdfPerSlideMain.size() - 1).push_back(saveName);

	std::cout << __LINE__ << std::endl;


	saveName = "pdfDir/" + dateStr + "/" + saveName;

	std::cout << __LINE__ << std::endl;

	quietSaveAs(spectCanv_p, saveName);

	std::cout << __LINE__ << std::endl;
	
	delete spectCanv_p;
	delete label_p;
	delete leg_p;

	std::cout << __LINE__ << std::endl;

	for(Int_t cI = 0; cI < nCentBins; ++cI){
	  for(Int_t sI = 0; sI < nSyst; ++sI){
	    scaleCentralAndErrorValues(jtPtUnfolded_RecoTrunc_PbPb_h[tI][cI][idPos][mI][absEtaPos][sI][bayesPos], 1./centBinsScalingFact.at(cI));
	  }
	}
      }
    }

    std::cout << __LINE__ << std::endl;

    for(Int_t usI = 0; usI < nSystSmooth; ++usI){
      TLatex* label_p = new TLatex();
      label_p->SetTextFont(43);
      label_p->SetTextSize(16);

      resModCompCanv_p[usI]->cd();
      resModCompPad_p[usI][0]->cd();
      Double_t maxVal = resModHist_p[usI]->GetMaximum()*1.2;
      resModHist_p[usI]->SetMaximum(maxVal);
      resModHist_p[usI]->SetMinimum(0.0);

      label_p->DrawLatex(500, .12, rStr.c_str());

      std::cout << __LINE__ << std::endl;
      legRes_p->Draw("SAME");
      std::cout << __LINE__ << std::endl;

      gPad->SetLogx();
      gPad->Modified();

      resModCompCanv_p[usI]->cd();
      resModCompPad_p[usI][1]->cd();
      gPad->SetLogx();

      drawWhiteBox(900, 1100, 0.2, .39);

      label_p->DrawLatex(200, .3, "200");
      label_p->DrawLatex(400, .3, "400");
      label_p->DrawLatex(600, .3, "600");
      label_p->DrawLatex(800, .3, "800");
    
      std::string saveName = "resErrComp_" + jetPbPbList.at(tI) + "_R" + rValStr + "_" + idNameStr + "_" + plotAbsEtaStr + "_" + plotBayesStr + "_" + systSmooth[usI] + "_"+ dateStr + ".pdf";

      if(usI == 0){
	slideTitles.push_back("Resolution Error Comp.");
	pdfPerSlide.push_back({});
      }
      
      pdfPerSlide.at(pdfPerSlide.size() - 1).push_back(saveName);
      
      saveName = "pdfDir/" + dateStr + "/" + saveName;

      quietSaveAs(resModCompCanv_p[usI], saveName);

      for(Int_t pI = 0; pI < nPad; ++pI){
	delete resModCompPad_p[usI][pI];
      }

      delete resModCompCanv_p[usI];
      delete resModHist_p[usI];
      std::cout << __LINE__ << std::endl;
      std::cout << __LINE__ << std::endl;
    }
    delete legRes_p;

    for(unsigned int i = 0; i < prevJERData.size(); ++i){
      delete prevJERData.at(i);
      prevJERData.at(i) = NULL;
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
  tex.InitTag("AllPlots");
  tex.SetAuthor("Christopher McGinn");
  tex.SetSlideTitles(slideTitles);
  tex.SetSlidePdfs(pdfPerSlide);
  if(!(tex.CreateTexSlides())){
    std::cout << "Warning: .tex slide creation failed" << std::endl;
  }


  texSlideCreator texMain;
  texMain.Clean();
  texMain.Init(outFileName);
  texMain.InitTag("MainPlots");
  texMain.SetAuthor("Christopher McGinn");
  texMain.SetSlideTitles(slideTitlesMain);
  texMain.SetSlidePdfs(pdfPerSlideMain);
  if(!(texMain.CreateTexSlides())){
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
