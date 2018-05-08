#ifndef READJSON_H
#define READJSON_H

#include <string>
#include <fstream>
#include <iostream>
#include <map>

#include "Utility/include/checkMakeDir.h"

class readJSON{
 public:
  readJSON(){};
  readJSON(const std::string inFileName);
  bool Init(const std::string inFileName);
  void Print();

  std::map<int, std::vector<int> > runRunToLumiMapLow;
  std::map<int, std::vector<int> > runRunToLumiMapHi;

  bool isGoodEvent(int run, int lumi);

  bool isInit = false;
};

readJSON::readJSON(const std::string inFileName)
{
  isInit = readJSON::Init(inFileName);
  return;
}

bool readJSON::Init(const std::string inFileName)
{
  if(!checkFile(inFileName)){
    std::cout << "READJSON ERROR: Given input file \'" << inFileName << "\' is not valid. return w/ initialization failed" << std::endl;
    return isInit;
  }

  std::ifstream file(inFileName.c_str());
  std::string tempStr;

  int pos = 0;
  while(std::getline(file, tempStr)){
    while(tempStr.find("\"") != std::string::npos){
      tempStr.replace(0, tempStr.find("\"")+1, "");
      int key = std::stoi(tempStr.substr(0, tempStr.find("\"")));
      tempStr.replace(0, tempStr.find("\"")+1, "");

      std::vector<int> tempLumiLow;
      std::vector<int> tempLumiHi;

      while(tempStr.substr(0,2).find("]]") == std::string::npos){
	if(std::string("0123456789").find(tempStr.substr(0,1)) == std::string::npos) tempStr.replace(0,1,"");
	else{
	  int posComma = tempStr.find(",");
	  int posBracket = tempStr.find("]");

	  if(posComma < posBracket && tempStr.find(",") != std::string::npos){
	    tempLumiLow.push_back(std::stoi(tempStr.substr(0, posComma)));
	    tempStr.replace(0, posComma+1,"");
	  }
	  else if(tempStr.find("]") != std::string::npos){
	    tempLumiHi.push_back(std::stoi(tempStr.substr(0, posBracket)));
	    tempStr.replace(0, posBracket,"");
	  }
	}

	if(tempStr.size() == 0 || tempStr.size() == 1) break;
      }
      
      runRunToLumiMapLow[key] = tempLumiLow;
      runRunToLumiMapHi[key] = tempLumiHi;
    }

    ++pos;
  }

  file.close();

  isInit = true;
  return isInit;
}


void readJSON::Print()
{
  for(std::map<int, std::vector<int> >::iterator it = runRunToLumiMapLow.begin(); it != runRunToLumiMapLow.end(); ++it){
    std::vector<int> lowVect = it->second;
    std::vector<int> hiVect = runRunToLumiMapHi[it->first];

    std::cout << "Run " << it->first << std::endl;
    for(unsigned int mI = 0; mI < lowVect.size(); ++mI){
      std::cout << " Lumi " << lowVect.at(mI) << "-" << hiVect.at(mI) << std::endl;
    }
  }

  return;
}



bool readJSON::isGoodEvent(int run, int lumi)
{
  if(!isInit){
    std::cout << "READJSON ERROR: Call of isGoodEvent without initialization. return false always" << std::endl;
    return false;
  }

  if(runRunToLumiMapLow.find(run) == runRunToLumiMapLow.end()) return false;

  std::vector<int> tempLow = runRunToLumiMapLow[run];
  std::vector<int> tempHi = runRunToLumiMapHi[run];
  
  bool isGoodLumi = false;
  for(unsigned int i = 0; i < tempLow.size(); ++i){
    if(lumi >= tempLow.at(i) && lumi <= tempHi.at(i)){
      isGoodLumi = true;
      break;
    }
  }

  return isGoodLumi;
}

#endif
