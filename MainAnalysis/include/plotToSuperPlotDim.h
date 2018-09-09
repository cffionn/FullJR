#ifndef PLOTTOSUPERPLOTDIM_H
#define PLOTTOSUPERPLOTDIM_H

class plotToSuperPlotDim{
 public:
  static const Int_t nPlots = 49;
  Int_t nPlotsX[nPlots];
  Int_t nPlotsY[nPlots];
  Double_t widthX[nPlots];
  
  plotToSuperPlotDim();
  int GetNPlotsX(const unsigned int nPlots);
  int GetNPlotsY(const unsigned int nPlots);
  double GetWidthX(const unsigned int nPlots);
};

plotToSuperPlotDim::plotToSuperPlotDim()
{
  nPlotsX[1-1] = 1;
  nPlotsY[1-1] = 1;
  widthX[1-1] = 0.8;

  nPlotsX[2-1] = 2;
  nPlotsY[2-1] = 1;
  widthX[2-1] = 0.45;

  nPlotsX[3-1] = 3;
  nPlotsY[3-1] = 1;
  widthX[3-1] = 0.325;

  nPlotsX[4-1] = 2;
  nPlotsY[4-1] = 2;
  widthX[4-1] = 0.45;


  nPlotsX[5-1] = 3;
  nPlotsY[5-1] = 2;
  widthX[5-1] = 0.325;
  nPlotsX[6-1] = 3;
  nPlotsY[6-1] = 2;
  widthX[6-1] = 0.325;

  nPlotsX[7-1] = 4;
  nPlotsY[7-1] = 2;
  widthX[7-1] = 0.24;
  nPlotsX[8-1] = 4;
  nPlotsY[8-1] = 2;
  widthX[8-1] = 0.24;

  nPlotsX[9-1] = 3;
  nPlotsY[9-1] = 3;
  widthX[9-1] = 0.3;

  nPlotsX[10-1] = 4;
  nPlotsY[10-1] = 3;
  widthX[10-1] = 0.24;
  nPlotsX[11-1] = 4;
  nPlotsY[11-1] = 3;
  widthX[11-1] = 0.24;
  nPlotsX[12-1] = 4;
  nPlotsY[12-1] = 3;
  widthX[12-1] = 0.24;

  nPlotsX[13-1] = 4;
  nPlotsY[13-1] = 4;
  widthX[13-1] = 0.24;
  nPlotsX[14-1] = 4;
  nPlotsY[14-1] = 4;
  widthX[14-1] = 0.24;
  nPlotsX[15-1] = 4;
  nPlotsY[15-1] = 4;
  widthX[15-1] = 0.24;
  nPlotsX[16-1] = 4;
  nPlotsY[16-1] = 4;
  widthX[16-1] = 0.24;

  nPlotsX[17-1] = 5;
  nPlotsY[17-1] = 4;
  widthX[17-1] = 0.19;
  nPlotsX[18-1] = 5;
  nPlotsY[18-1] = 4;
  widthX[18-1] = 0.19;
  nPlotsX[19-1] = 5;
  nPlotsY[19-1] = 4;
  widthX[19-1] = 0.19;
  nPlotsX[20-1] = 5;
  nPlotsY[20-1] = 4;
  widthX[20-1] = 0.19;


  //EDITING
  nPlotsX[21-1] = 6;
  nPlotsY[21-1] = 5;
  widthX[21-1] = 0.19;
  nPlotsX[22-1] = 6;
  nPlotsY[22-1] = 5;
  nPlotsX[23-1] = 6;
  nPlotsY[23-1] = 5;
  nPlotsX[24-1] = 6;
  nPlotsY[24-1] = 5;
  nPlotsX[25-1] = 6;
  nPlotsY[25-1] = 5;
  nPlotsX[26-1] = 6;
  nPlotsY[26-1] = 5;
  nPlotsX[27-1] = 6;
  nPlotsY[27-1] = 5;
  nPlotsX[28-1] = 6;
  nPlotsY[28-1] = 5;
  nPlotsX[29-1] = 6;
  nPlotsY[29-1] = 5;
  nPlotsX[30-1] = 6;
  nPlotsY[30-1] = 5;

  nPlotsX[31-1] = 6;
  nPlotsY[31-1] = 6;
  nPlotsX[32-1] = 6;
  nPlotsY[32-1] = 6;
  nPlotsX[33-1] = 6;
  nPlotsY[33-1] = 6;
  nPlotsX[34-1] = 6;
  nPlotsY[34-1] = 6;
  nPlotsX[35-1] = 6;
  nPlotsY[35-1] = 6;
  nPlotsX[36-1] = 6;
  nPlotsY[36-1] = 6;

  nPlotsX[37-1] = 7;
  nPlotsY[37-1] = 6;
  nPlotsX[38-1] = 7;
  nPlotsY[38-1] = 6;
  nPlotsX[39-1] = 7;
  nPlotsY[39-1] = 6;
  nPlotsX[40-1] = 7;
  nPlotsY[40-1] = 6;
  nPlotsX[41-1] = 7;
  nPlotsY[41-1] = 6;
  nPlotsX[42-1] = 7;
  nPlotsY[42-1] = 6;

  nPlotsX[42-1] = 7;
  nPlotsY[42-1] = 7;
  nPlotsX[43-1] = 7;
  nPlotsY[43-1] = 7;
  nPlotsX[44-1] = 7;
  nPlotsY[44-1] = 7;
  nPlotsX[45-1] = 7;
  nPlotsY[45-1] = 7;
  nPlotsX[46-1] = 7;
  nPlotsY[46-1] = 7;
  nPlotsX[47-1] = 7;
  nPlotsY[47-1] = 7;
  nPlotsX[48-1] = 7;
  nPlotsY[48-1] = 7;
  nPlotsX[49-1] = 7;
  nPlotsY[49-1] = 7;

  return;
}

int plotToSuperPlotDim::GetNPlotsX(const unsigned int nPlots)
{
  int tempNPlotsX = -1;
  if(nPlots <= 0) std::cout << "Warning plotToSuperPlotDum: Given nPlots \'" << nPlots << "\' is less than 1, return -1" << std::endl;
  else if(nPlots > 49) std::cout << "Warning plotToSuperPlotDum: Given nPlots \'" << nPlots << "\' is greater than 49, return -1" << std::endl;
  else tempNPlotsX = nPlotsX[nPlots-1];
  return tempNPlotsX;
}

int plotToSuperPlotDim::GetNPlotsY(const unsigned int nPlots)
{
  int tempNPlotsY = -1;
  if(nPlots <= 0) std::cout << "Warning plotToSuperPlotDum: Given nPlots \'" << nPlots << "\' is less than 1, return -1" << std::endl;
  else if(nPlots > 49) std::cout << "Warning plotToSuperPlotDum: Given nPlots \'" << nPlots << "\' is greater than 49, return -1" << std::endl;
  else tempNPlotsY = nPlotsY[nPlots-1];
  return tempNPlotsY;
}

double plotToSuperPlotDim::GetWidthX(const unsigned int nPlots)
{
  double tempWidthX = -1;
  if(nPlots <= 0) std::cout << "Warning plotToSuperPlotDum: Given nPlots \'" << nPlots << "\' is less than 1, return -1" << std::endl;
  else if(nPlots > 49) std::cout << "Warning plotToSuperPlotDum: Given nPlots \'" << nPlots << "\' is greater than 49, return -1" << std::endl;
  else tempWidthX = widthX[nPlots-1];
  return tempWidthX;
}


#endif
