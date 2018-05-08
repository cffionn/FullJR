#ifndef VECTORSTRINGUTILITY
#define VECTORSTRINGUTILITY

#include <string>
#include <vector>

bool vectorContainsString(std::vector<std::string>* vect, const std::string inStr, const bool exactStr = false, const bool slashExactStr = false)
{
  bool hasString = false;

  for(unsigned int vI = 0; vI < vect->size(); ++vI){
    if(vect->at(vI).find(inStr) != std::string::npos){
      if(vect->at(vI).size() == inStr.size() || !exactStr || !slashExactStr){
	hasString = true;
	break;
      }
    }
  }

  return hasString;
}

void removeVectorDuplicates(std::vector<std::string>* vect)
{
  unsigned int pos = 0;
  while(pos < vect->size()){
    bool isGood = true;
    for(unsigned int p = pos+1; p < vect->size(); ++p){
      if(vect->at(p).size() == vect->at(pos).size() && vect->at(pos).find(vect->at(p)) != std::string::npos){
        vect->erase(vect->begin() + p);
        isGood = false;
        break;
      }
    }
    if(isGood) ++pos;
  }

  return;
}

#endif
