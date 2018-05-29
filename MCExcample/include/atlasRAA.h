#ifndef ATLASRAA_h
#define ATLASRAA_h

#include "TH1D.h"
#include "TColor.h"

class atlasRAA{
 public:
  //extracted with GIMP color picker from plot (url below)
  static const int rVal = 80;
  static const int gVal = 20;
  static const int bVal = 60;

  static const int markerStyle = 21;

  static const int nPoints = 15;
  //Extracted using https://automeris.io/WebPlotDigitizer/ with https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/HION-2017-10/fig_04a.png
  float bins[nPoints+1] = {100.10092096608406,
			   112.16007807513813,
			   126.03907442914226,
			   141.22300066829993,
			   158.2361343741133,
			   177.81670858322533,
			   199.82023686584935,
			   223.89257912481887,
			   250.86491625378744,
			   281.90763054025047,
			   316.79165561923617,
			   354.955543714172,
			   398.87871835517285,
			   500.7734844332488,
			   632.3756370274339,
			   999.63791851118};

  float yVals[nPoints] = {0.44293555384285344,
			  0.46934886418713295,
			  0.4902492511342238,
			  0.5078425514671413,
			  0.5210275153547279,
			  0.5408255400970944,
			  0.5628294016224764,
			  0.5782190899195119,
			  0.5671509729295947,
			  0.5693123147699372,
			  0.5957245127406503,
			  0.5956811301715441,
			  0.5989203619981408,
			  0.5646592561399051,
			  0.6262524929086185};

  //These are not the errors, but rather the maxVal within statistical error. Where error is not visibible, top of marker is taken. this will lead to underweighting of low pt points in fits, but is also the best that can be done without hepdata
  float yErrs[nPoints] = {0.4484473648664755,
			  0.4737583130060308,
			  0.495761062157846,
			  0.5133543624907634,
			  0.5265404387519167,
			  0.546338463494283,
			  0.5683412126460982,
			  0.5848332631478584,
			  0.5759709829409568,
			  0.5814382990219059,
			  0.6122599458115162,
			  0.6210354608802054,
			  0.6297865037304242,
			  0.6252891773997475,
			  0.7816855637747601};

  float yErrsSyst[nPoints] = {0.45579399141630883,
			      0.4815450643776822,
			      0.5030042918454934,
			      0.5201716738197423,
			      0.533047210300429,
			      0.5536480686695278,
			      0.5776824034334762,
			      0.5939914163090128,
			      0.5819742489270384,
			      0.5854077253218882,
			      0.6145922746781114,
			      0.6154506437768239,
			      0.6163090128755363,
			      0.5819742489270384,
			      0.6437768240343344};
  


  atlasRAA(){return;}
  void SetHistogram(TH1D** inHist_p, const std::string name);
  void StyleHistogram(TH1D* inHist_p);
};


void atlasRAA::SetHistogram(TH1D** inHist_p, const std::string name)
{
  (*inHist_p) = new TH1D(name.c_str(), "", nPoints, bins);

  for(Int_t pI = 0; pI < nPoints; ++pI){
    (*inHist_p)->SetBinContent(pI+1, yVals[pI]);
    (*inHist_p)->SetBinError(pI+1, yErrs[pI] - yVals[pI]);
  }

  TColor col;
  (*inHist_p)->SetMarkerColor((Int_t)col.GetColor(rVal, gVal, bVal));
  (*inHist_p)->SetLineColor((Int_t)col.GetColor(rVal, gVal, bVal));
  (*inHist_p)->SetMarkerStyle(markerStyle);
  (*inHist_p)->SetMarkerSize(1.);

  return;
}


void atlasRAA::StyleHistogram(TH1D* inHist_p)
{
  TColor col;
  inHist_p->SetMarkerColor((Int_t)col.GetColor(rVal, gVal, bVal));
  inHist_p->SetLineColor((Int_t)col.GetColor(rVal, gVal, bVal));
  inHist_p->SetMarkerStyle(markerStyle);
  inHist_p->SetMarkerSize(1.);

  return;
}

#endif
