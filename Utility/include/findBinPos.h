#ifndef FINDBINPOS_H
#define FINDBINPOS_H

#include <vector>

Int_t findBinPos(Float_t inVal, Int_t arrSize, Float_t arr[])
{
  Int_t pos = -1;

  for(Int_t i = 0; i < arrSize-1; ++i){
    if(inVal >= arr[i] && inVal < arr[i+1]){
      pos = i;
      break;
    }
  }
  if(pos == -1 && arr[arrSize-1] == inVal) pos = arrSize-1;

  return pos;
}

Int_t findBinPos(Double_t inVal, Int_t arrSize, Double_t arr[])
{
  Int_t pos = -1;

  for(Int_t i = 0; i < arrSize-1; ++i){
    if(inVal >= arr[i] && inVal < arr[i+1]){
      pos = i;
      break;
    }
  }
  if(pos == -1 && arr[arrSize-1] == inVal) pos = arrSize-1;

  return pos;
}

Int_t findBinPos(Float_t inVal, std::vector<float> arr)
{
  Int_t pos = -1;

  for(unsigned int i = 0; i < arr.size()-1; ++i){
    if(inVal >= arr[i] && inVal < arr[i+1]){
      pos = i;
      break;
    }
  }
  if(pos == -1 && arr[arr.size()-1] == inVal) pos = arr.size()-1;

  return pos;
}

Int_t findBinPos(Float_t inVal, std::vector<double> arr)
{
  Int_t pos = -1;

  for(unsigned int i = 0; i < arr.size()-1; ++i){
    if(inVal >= arr[i] && inVal < arr[i+1]){
      pos = i;
      break;
    }
  }
  if(pos == -1 && arr[arr.size()-1] == inVal) pos = arr.size()-1;

  return pos;
}

Int_t findBinPos(Float_t inVal, std::vector<double> arrLow, std::vector<double> arrHi)
{
  Int_t pos = -1;

  for(unsigned int i = 0; i < arrLow.size()-1; ++i){
    if(inVal >= arrLow[i] && inVal < arrHi[i]){
      pos = i;
      break;
    }
  }
  if(pos == -1 && arrHi[arrHi.size()-1] == inVal) pos = arrHi.size()-1;

  return pos;
}

#endif
