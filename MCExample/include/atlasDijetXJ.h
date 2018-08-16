#ifndef ATLASDIJETXJ_h
#define ATLASDIJETXJ_h

#include "TH1D.h"
#include "TColor.h"

class atlasDijetXJ{
 public:
  //extracted with GIMP color picker from plot (url below)
  static const int rVal = 80;
  static const int gVal = 20;
  static const int bVal = 60;

  static const int markerStyle = 21;

  static const int nPoints = 10;
  //Extracted using https://automeris.io/WebPlotDigitizer/ with https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/HION-2012-11/figaux_05.png
  float bins[nPoints+1] = {0.31594533029612765,
			   0.35501138952164024,
			   0.3980637813211844, 
			   0.44669703872437355,
			   0.5009111617312073, 
			   0.5623006833712985, 
			   0.6308656036446468, 
			   0.7074031890660594, 
			   0.7943052391799544,
			   0.8915717539863326,
			   1.0000000000000002};


  float yVals[nPoints] = {0.05175413376974802,
			  0.2744856285615289,
			  0.7925476037773409,
			  1.755630257056712,
			  2.512439090017766,
			  2.3469500777731085,
			  1.7729403026885167,
			  1.4255138185490712,
			  1.3734546741269877,
			  1.3457364855400993};

  //These are not the errors, but rather the maxVal within statistical error. Where error is not visibible, top of marker is taken. this will lead to underweighting of low pt points in fits, but is also the best that can be done without hepdata
  float yErrs[nPoints] = {0.08816190076003849,
			  0.3189840104385522,
			  0.8370505930660301,
			  1.8041739463771007,
			  2.556932864483122,
			  2.3874077596184327,
			  1.813393377122174,
			  1.4619215855393621,
			  1.41391235597231,
			  1.3861895599737561};
			  
  float yErrsSyst[nPoints] = {0.11647905286359883,
			      0.48080091558484694,
			      1.140448651318457,
			      2.2087046907136694,
			      2.9169744417660026,
			      2.6018228763518154,
			      1.9549837450516403,
			      1.5711541013335695,
			      1.527180964386551,
			      1.5035034758313617};
  


  atlasDijetXJ(){return;}
  void SetHistogram(TH1D** inHist_p, const std::string name);
  void StyleHistogram(TH1D* inHist_p);
};


void atlasDijetXJ::SetHistogram(TH1D** inHist_p, const std::string name)
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


void atlasDijetXJ::StyleHistogram(TH1D* inHist_p)
{
  TColor col;
  inHist_p->SetMarkerColor((Int_t)col.GetColor(rVal, gVal, bVal));
  inHist_p->SetLineColor((Int_t)col.GetColor(rVal, gVal, bVal));
  inHist_p->SetMarkerStyle(markerStyle);
  inHist_p->SetMarkerSize(1.);

  return;
}

#endif
