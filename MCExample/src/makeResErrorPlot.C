//cpp dependencisef
#include <iostream>
#include <string>

//ROOT dependencies
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TMath.h"
#include "TDatime.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"

//Non-local (Utility) dependencies
#include "Utility/include/plotUtilities.h"
#include "Utility/include/vanGoghPalette.h"

int makeResErrorPlot(const std::string inFileName, const Int_t rVal, const std::string paraStr, const std::string lowStr, const std::string hiStr, const bool doFake)
{
  if(rVal != 4 && rVal != 8){
    std::cout << "Given rVal \'" << rVal << "\' is invalid. please pick 4 or 8" << std::endl;
    return 1;
  }

  std::string lowStr2 = lowStr;
  std::string hiStr2 = hiStr;

  lowStr2.replace(lowStr2.find("p"), 1, ".");
  hiStr2.replace(hiStr2.find("p"), 1, ".");
  
  std::string rStr = "R0p" + std::to_string(rVal);

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");

  TH1D* gen_h = (TH1D*)inFile_p->Get(("genJtPt_" + rStr + paraStr + "_h").c_str());
  TH1D* unfold_h = (TH1D*)inFile_p->Get(("unfoldedJtPt_" + rStr + "_Bayes3" + paraStr + "_h").c_str());
  TH1D* unfoldLow_h = NULL;
  TH1D* unfoldHi_h = NULL;

  if(doFake){
    unfoldLow_h = (TH1D*)inFile_p->Get(("unfoldedJtPtFake_" + rStr + "_Bayes3" + paraStr + "_h").c_str());
    unfoldHi_h = (TH1D*)unfoldLow_h->Clone("hi_h");

    for(Int_t bIX = 0; bIX < unfoldLow_h->GetNbinsX(); ++bIX){
      if(unfoldLow_h->GetBinContent(bIX+1) > unfold_h->GetBinContent(bIX+1)){
	Double_t newVal = unfold_h->GetBinContent(bIX+1) - (unfoldLow_h->GetBinContent(bIX+1) - unfold_h->GetBinContent(bIX+1));
	Double_t errVal = newVal*unfoldLow_h->GetBinError(bIX+1)/unfold_h->GetBinContent(bIX+1);

	unfoldLow_h->SetBinContent(bIX+1, newVal);
	unfoldLow_h->SetBinError(bIX+1, errVal);
      }
      if(unfoldHi_h->GetBinContent(bIX+1) < unfold_h->GetBinContent(bIX+1)){
	Double_t newVal = unfold_h->GetBinContent(bIX+1) + (unfold_h->GetBinContent(bIX+1) - unfoldHi_h->GetBinContent(bIX+1));
	Double_t errVal = newVal*unfoldHi_h->GetBinError(bIX+1)/unfold_h->GetBinContent(bIX+1);

	unfoldHi_h->SetBinContent(bIX+1, newVal);
	unfoldHi_h->SetBinError(bIX+1, errVal);
      }
    }
  }
  else{
    unfoldLow_h = (TH1D*)inFile_p->Get(("unfoldedJtPtRes" + lowStr + "_" + rStr + "_Bayes3" + paraStr + "_h").c_str());
    unfoldHi_h = (TH1D*)inFile_p->Get(("unfoldedJtPtRes" + hiStr + "_" + rStr + "_Bayes3" + paraStr + "_h").c_str());
  }

  const Double_t yFrac = .4;
  const Double_t marginFrac = 0.15;
  const Double_t marginFracBottom = marginFrac/yFrac;

  vanGoghPalette col;
  const Int_t nStyles = 4;
  const Int_t styles[nStyles] = {20, 21, 47, 34};

  for(Int_t bI = 0; bI < gen_h->GetNbinsX(); ++bI){
    gen_h->SetBinContent(bI+1, gen_h->GetBinContent(bI+1)/(4.*gen_h->GetBinWidth(bI+1)));
    gen_h->SetBinError(bI+1, gen_h->GetBinError(bI+1)/(4.*gen_h->GetBinWidth(bI+1)));
  }

  for(Int_t bI = 0; bI < unfold_h->GetNbinsX(); ++bI){
    unfold_h->SetBinContent(bI+1, unfold_h->GetBinContent(bI+1)/(4.*unfold_h->GetBinWidth(bI+1)));
    unfold_h->SetBinError(bI+1, unfold_h->GetBinError(bI+1)/(4.*unfold_h->GetBinWidth(bI+1)));
  }

  for(Int_t bI = 0; bI < unfoldLow_h->GetNbinsX(); ++bI){
    unfoldLow_h->SetBinContent(bI+1, unfoldLow_h->GetBinContent(bI+1)/(4.*unfoldLow_h->GetBinWidth(bI+1)));
    unfoldLow_h->SetBinError(bI+1, unfoldLow_h->GetBinError(bI+1)/(4.*unfoldLow_h->GetBinWidth(bI+1)));
  }

  for(Int_t bI = 0; bI < unfoldHi_h->GetNbinsX(); ++bI){
    unfoldHi_h->SetBinContent(bI+1, unfoldHi_h->GetBinContent(bI+1)/(4.*unfoldHi_h->GetBinWidth(bI+1)));
    unfoldHi_h->SetBinError(bI+1, unfoldHi_h->GetBinError(bI+1)/(4.*unfoldHi_h->GetBinWidth(bI+1)));
  }

  prettyTH1(gen_h, 0.6, styles[0], 1);
  prettyTH1(unfold_h, 0.6, styles[1], col.getColor(0));
  prettyTH1(unfoldLow_h, 0.6, styles[2], col.getColor(1));
  prettyTH1(unfoldHi_h, 0.6, styles[3], col.getColor(2));

  gen_h->GetXaxis()->SetNdivisions(505);
  gen_h->GetYaxis()->SetNdivisions(505);
  unfold_h->GetXaxis()->SetNdivisions(505);
  unfold_h->GetYaxis()->SetNdivisions(505);

  gen_h->GetXaxis()->SetLabelFont(43);
  gen_h->GetYaxis()->SetLabelFont(43);
  gen_h->GetXaxis()->SetTitleFont(43);
  gen_h->GetYaxis()->SetTitleFont(43);

  gen_h->GetXaxis()->SetLabelSize(16);
  gen_h->GetYaxis()->SetLabelSize(14);
  gen_h->GetXaxis()->SetTitleSize(20);
  gen_h->GetYaxis()->SetTitleSize(18);

  unfold_h->GetXaxis()->SetLabelFont(43);
  unfold_h->GetYaxis()->SetLabelFont(43);
  unfold_h->GetXaxis()->SetTitleFont(43);
  unfold_h->GetYaxis()->SetTitleFont(43);

  unfold_h->GetXaxis()->SetLabelSize(16);
  unfold_h->GetYaxis()->SetLabelSize(14);
  unfold_h->GetXaxis()->SetTitleSize(20);
  unfold_h->GetYaxis()->SetTitleSize(18);

  gen_h->GetYaxis()->SetTitleOffset(gen_h->GetXaxis()->GetTitleOffset()*1.5);
  gen_h->GetXaxis()->SetTitleOffset(gen_h->GetYaxis()->GetTitleOffset());

  unfold_h->GetXaxis()->SetTitleOffset(gen_h->GetYaxis()->GetTitleOffset()*2.);
  unfold_h->GetYaxis()->SetTitleOffset(gen_h->GetYaxis()->GetTitleOffset());

  unfold_h->SetTitle("");
  unfold_h->GetYaxis()->SetTitle("#frac{Unfold}{Gen.}");
  gen_h->GetYaxis()->SetTitle("#frac{1}{N_{evt}} #frac{d^{2}N_{jet}}{dp_{T}d#eta}");


  TCanvas* canv_p = new TCanvas("canv_p", "", 500, 500);
  TLegend* leg_p = new TLegend(0.5, 0.6, 0.9, 0.9);
  TLatex* label_p = new TLatex();
  label_p->SetTextFont(43);
  label_p->SetTextSize(18);
  label_p->SetNDC();

  leg_p->SetTextFont(43);
  leg_p->SetTextSize(18);
  leg_p->SetBorderSize(0);
  leg_p->SetFillStyle(0);

  leg_p->AddEntry(gen_h, "Gen. PYTHIA", "P L");
  leg_p->AddEntry(unfold_h, "Unfolded", "P L");
  if(doFake){
    leg_p->AddEntry(unfoldLow_h, "Fake (Low)", "P L");
    leg_p->AddEntry(unfoldHi_h, "Fake (Hi)", "P L");
  }
  else{
    leg_p->AddEntry(unfoldLow_h, ("Unfolded (" + lowStr2 + "#sigma)").c_str(), "P L");
    leg_p->AddEntry(unfoldHi_h, ("Unfolded (" + hiStr2 + "#sigma)").c_str(), "P L");
  }

  canv_p->SetTopMargin(0.01);
  canv_p->SetRightMargin(0.001);
  canv_p->SetLeftMargin(0.01);
  canv_p->SetBottomMargin(0.01);

  TPad* pads[2];
  pads[0] = new TPad("pad0", "", 0.0, yFrac, 1.0, 1.0);
  pads[0]->SetRightMargin(0.001);
  pads[0]->SetBottomMargin(0.001);
  pads[0]->SetLeftMargin(marginFrac);
  pads[0]->SetTopMargin(marginFrac/2.);
  pads[0]->Draw("SAME");
  gStyle->SetOptStat(0);

  canv_p->cd();
  pads[1] = new TPad("pad1", "", 0.0, 0, 1.0, yFrac);
  pads[1]->SetTopMargin(0.001);
  pads[1]->SetRightMargin(0.001);
  pads[1]->SetBottomMargin(marginFracBottom);
  pads[1]->SetLeftMargin(marginFrac);
  pads[1]->Draw("SAME");
  gStyle->SetOptStat(0);

  canv_p->cd();
  pads[0]->cd();

  Double_t max = TMath::Max(TMath::Max(gen_h->GetMaximum(), unfold_h->GetMaximum()), TMath::Max(unfoldLow_h->GetMaximum(), unfoldHi_h->GetMaximum()));
  Double_t min = TMath::Min(TMath::Min(gen_h->GetMinimum(), unfold_h->GetMinimum()), TMath::Min(unfoldLow_h->GetMinimum(), unfoldHi_h->GetMinimum()));

  gen_h->SetMaximum(max*5.);
  gen_h->SetMinimum(min/5.);

  gen_h->DrawCopy("HIST E1 P");
  unfold_h->DrawCopy("HIST E1 P SAME");
  unfoldLow_h->DrawCopy("HIST E1 P SAME");
  unfoldHi_h->DrawCopy("HIST E1 P SAME");
  gPad->SetLogy();
  gPad->RedrawAxis();
  gPad->SetTicks(1, 2);

  leg_p->Draw("SAME");
  if(doFake){
    if(rVal == 4) label_p->DrawLatex(0.52, 0.54, "#sigma_{RC}=20.0");
    else if(rVal == 8) label_p->DrawLatex(0.52, 0.54, "#sigma_{RC}=40.0");
  }
  else{
    if(rVal == 4) label_p->DrawLatex(0.52, 0.54, "C,S,N = 0.06, 1.0, 20.0");
    else if(rVal == 8) label_p->DrawLatex(0.52, 0.54, "C,S,N = 0.06, 1.0, 40.0");
  }

  if(paraStr.size() == 0) label_p->DrawLatex(0.1, 0.94, "Generator Study, Spectra Matches Response");
  else label_p->DrawLatex(0.1, 0.94, "Generator Study, Spectra Independent of Response");

  canv_p->cd();
  pads[1]->cd();

  unfold_h->Divide(gen_h);
  unfoldLow_h->Divide(gen_h);
  unfoldHi_h->Divide(gen_h);

  unfold_h->SetMaximum(1.29);
  unfold_h->SetMinimum(0.71);

  unfold_h->DrawCopy("HIST E1 P");
  unfoldLow_h->Smooth();
  unfoldHi_h->Smooth();
  unfoldLow_h->DrawCopy("HIST E1 P SAME");
  unfoldHi_h->DrawCopy("HIST E1 P SAME");
  gPad->RedrawAxis();
  gPad->SetTicks(1, 2);

  std::string saveName = "pdfDir/resErr_" + rStr + paraStr + "_Err" + lowStr + "to" + hiStr + "_" + dateStr + ".pdf";
  if(doFake) saveName = "pdfDir/resErr_" + rStr + paraStr + "_Fake_" + dateStr + ".pdf";
  quietSaveAs(canv_p, saveName);

  delete label_p;
  delete leg_p;
  delete canv_p;

  inFile_p->Close();
  delete inFile_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/makeResErrorPlot.exe <inFileName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  
  retVal += makeResErrorPlot(argv[1], 4, "", "0p9", "1p1", false);
  retVal += makeResErrorPlot(argv[1], 8, "", "0p9", "1p1", false);
  retVal += makeResErrorPlot(argv[1], 4, "_Parallel", "0p9", "1p1", false);
  retVal += makeResErrorPlot(argv[1], 8, "_Parallel", "0p9", "1p1", false);

  retVal += makeResErrorPlot(argv[1], 4, "", "0p85", "1p15", false);
  retVal += makeResErrorPlot(argv[1], 8, "", "0p85", "1p15", false);
  retVal += makeResErrorPlot(argv[1], 4, "_Parallel", "0p85", "1p15", false);
  retVal += makeResErrorPlot(argv[1], 8, "_Parallel", "0p85", "1p15", false);
  
  retVal += makeResErrorPlot(argv[1], 4, "", "0p85", "1p15", true);
  retVal += makeResErrorPlot(argv[1], 8, "", "0p85", "1p15", true);
  retVal += makeResErrorPlot(argv[1], 4, "_Parallel", "0p85", "1p15", true);
  retVal += makeResErrorPlot(argv[1], 8, "_Parallel", "0p85", "1p15", true);

  return retVal;
}
