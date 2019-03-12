#ifndef CANVNDCTOXY_H
#define CANVNDCTOXY_H

//cpp dependencies
#include <iostream>

//ROOT dependencies
#include "TPad.h"
#include "TH1D.h"

//Non-Local dependencies
#include "Utility/include/getLinBins.h"
#include "Utility/include/getLogBins.h"

class canvNDCToXY
{
 public:
  canvNDCToXY(double inCanvDimX, double inCanvDimY, double inCanvMargLeft, double inCanvMargRight, double inCanvMargTop, double inCanvMargBottom, double inXMin, double inXMax, double inYMin, double inYMax);
  canvNDCToXY(TPad* canv_p, TH1D* hist_p);

  double getXRelFromAbs(double xAbs, bool isLog);
  double getYRelFromAbs(double yAbs, bool isLog);

 private:
  double canvDimX_;
  double canvDimY_;

  double canvMargLeft_;
  double canvMargRight_;
  double canvMargTop_;
  double canvMargBottom_;

  double xMin_;
  double xMax_;
  double yMin_;
  double yMax_;

  static const int nBins_ = 100;
  double binsLinX_[nBins_+1];
  double binsLinY_[nBins_+1];
  double binsLogX_[nBins_+1];
  double binsLogY_[nBins_+1];

  double binsLinRelX_[nBins_+1];
  double binsLinRelY_[nBins_+1];

  void defineBins();
};


canvNDCToXY::canvNDCToXY(double inCanvDimX, double inCanvDimY, double inCanvMargLeft, double inCanvMargRight, double inCanvMargTop, double inCanvMargBottom, double inXMin, double inXMax, double inYMin, double inYMax)
{
  canvDimX_ = inCanvDimX;
  canvDimY_ = inCanvDimY;

  canvMargLeft_ = inCanvMargLeft;
  canvMargRight_ = inCanvMargRight;
  canvMargTop_ = inCanvMargTop;
  canvMargBottom_ = inCanvMargBottom;

  xMin_ = inXMin;
  xMax_ = inXMax;
  yMin_ = inYMin;
  yMax_ = inYMax;

  defineBins();

  return;
}

canvNDCToXY::canvNDCToXY(TPad* canv_p, TH1D* hist_p)
{
  canvDimX_ = canv_p->GetWw();
  canvDimY_ = canv_p->GetWh(); 
  canvMargLeft_ = canv_p->GetLeftMargin();
  canvMargRight_ = canv_p->GetRightMargin();
  canvMargTop_ = canv_p->GetTopMargin();
  canvMargBottom_ = canv_p->GetBottomMargin();
  
  xMin_ = hist_p->GetBinLowEdge(1);
  xMax_ = hist_p->GetBinLowEdge(hist_p->GetNbinsX()+1);
  yMin_ = hist_p->GetMinimum();
  yMax_ = hist_p->GetMaximum();

  defineBins();

  return;
}

void canvNDCToXY::defineBins()
{
  getLinBins(xMin_, xMax_, nBins_, binsLinX_);
  getLinBins(yMin_, yMax_, nBins_, binsLinY_);
  getLogBins(xMin_, xMax_, nBins_, binsLogX_);
  getLogBins(yMin_, yMax_, nBins_, binsLogY_);

  getLinBins(canvMargLeft_, 1.0 - canvMargRight_, nBins_, binsLinRelX_);
  getLinBins(canvMargBottom_, 1.0 - canvMargTop_, nBins_, binsLinRelY_);

  return;
}

double canvNDCToXY::getXRelFromAbs(double xAbs, bool isLog = false)
{
  double relPos = -1;
  int arrPos = -1;
  if(!isLog){
    for(int bI = 0; bI < nBins_; ++bI){
      if(xAbs >= binsLinX_[bI] && xAbs < binsLinX_[bI+1]){
	arrPos = bI;
      }
    }
  } 
  else{
    for(int bI = 0; bI < nBins_; ++bI){
      if(xAbs >= binsLogX_[bI] && xAbs < binsLogX_[bI+1]){
	arrPos = bI;
      }
    }
  } 

  relPos = binsLinRelX_[arrPos];
  return relPos;
}

double canvNDCToXY::getYRelFromAbs(double yAbs, bool isLog = false)
{
  double relPos = -1;
  int arrPos = -1;
  if(!isLog){
    for(int bI = 0; bI < nBins_; ++bI){
      if(yAbs >= binsLinY_[bI] && yAbs < binsLinY_[bI+1]){
	arrPos = bI;
      }
    }
  } 
  else{
    for(int bI = 0; bI < nBins_; ++bI){
      if(yAbs >= binsLogY_[bI] && yAbs < binsLogY_[bI+1]){
	arrPos = bI;
      }
    }
  } 

  relPos = binsLinRelY_[arrPos];
  return relPos;
}


#endif
