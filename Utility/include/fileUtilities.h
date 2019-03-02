#ifndef FILEUTILITIES_H
#define FILEUTILITIES_H

//cpp dependencies
#include <iostream>
#include <string>

//Local Utilities
#include "Utility/include/checkMakeDir.h"

bool fileIsGood(const std::string name, const std::string ext)
{
  if(!checkFile(name) || name.find(ext) == std::string::npos){
    std::cout << "File \'" << name << "\' is invalid, ext \'" << ext << "\'. return false" << std::endl;
    return false;
  }
  return true;
}

#endif
