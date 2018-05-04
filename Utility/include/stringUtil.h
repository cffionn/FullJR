#ifndef STRINGUTIL_H
#define STRINGUTIL_H

#include <string>

std::string removeAllWhiteSpace(std::string inStr)
{
  while(inStr.find(" ") != std::string::npos){
    inStr.replace(inStr.find(" "), 1, "");
  }

  return inStr;
}


std::string returnAllCapsString(std::string inStr)
{
  const std::string lowStr = "abcdefghijklmnopqrstuvwxyz";
  const std::string hiStr = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";

  for(unsigned int lowIter = 0; lowIter < lowStr.size(); ++lowIter){
    while(inStr.find(lowStr.substr(lowIter, 1)) != std::string::npos){
      inStr.replace(inStr.find(lowStr.substr(lowIter, 1)), 1, hiStr.substr(lowIter, 1));
    }
  }

  return inStr;
}


bool isStrFromCharSet(const std::string inStr, const std::string charSet)
{
  for(unsigned int iter = 0; iter < inStr.size(); ++iter){
    if(charSet.find(inStr.substr(iter, 1)) == std::string::npos){
      return false;
    }
  }

  return true;
}


bool isStrAllAlpha(std::string inStr){return isStrFromCharSet(returnAllCapsString(inStr), "ABCDEFGHIJKLMNOPQRSTUVWXYZ");}
bool isStrInt(std::string inStr){return isStrFromCharSet(inStr, "-0123456789");}
bool isStrFloatOrDouble(std::string inStr){return isStrFromCharSet(inStr, ".-0123456789");}

bool isStrTrueOrFalse(std::string inStr)
{
  inStr = returnAllCapsString(inStr);
  if(!isStrAllAlpha(inStr)) return false;

  if(inStr.size() == 4 && inStr.find("TRUE") != std::string::npos) return true;
  if(inStr.size() == 5 && inStr.find("FALSE") != std::string::npos) return true;

  return false;
}


bool strToTrueOrFalse(std::string inStr)
{
  inStr = returnAllCapsString(inStr);

  if(inStr.size() == 4 && inStr.find("TRUE") != std::string::npos) return true;
  else if(inStr.size() == 5 && inStr.find("FALSE") != std::string::npos) return false;

  std::cout << "Call to strToTrueOrFalse is invalid; \'" << inStr << "\' neither true or false str. return false but you really ought to fix this" << std::endl;
  return false;
}

#endif
