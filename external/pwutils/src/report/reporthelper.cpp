
#include <vector>
#include "pwutils/report/reporthelper.h"

namespace pw{

fs::path createDirectory(const std::string& dir_name, bool overwrite)
{
    fs::path dir_path = fs::current_path() / dir_name;
    createDirectory(dir_path,overwrite);
    return dir_path;
}

void createDirectory(const std::filesystem::path& dir_path,bool overwrite)
{
	if(overwrite)
	    clearDirectory(dir_path);
	if(!fs::exists(dir_path)){
        fs::create_directory(dir_path);
	}
}

void clearDirectory(const std::filesystem::path& dir_path)
{
	if(fs::exists(dir_path)){
        fs::remove_all(dir_path);
	}
}

fs::path filePath(const fs::path& dir_path,
		const std::string& nm,int repNum,const std::string& extension) 
{
	fs::path local_path(nm+"_"+std::to_string(repNum)+"." + extension);
	fs::path full_path = dir_path / local_path;
	return full_path;
}

fs::path filePath(const fs::path& dir_path,
		const std::string& nm,const std::string& extension) 
{
	fs::path local_path(nm+"."+extension);
	fs::path full_path = dir_path / local_path;
	return full_path;
}



void AdjustPhase(std::vector<double>& in,int sz)
{
  const double MAX_DIFF = PI/2.0;
  int count = 0;
  std::vector<int> indxVec;
  std::vector<int> countVec;
  int stindx = sz/2;
  while(stindx > 0){
    stindx--;
    double diff = in[stindx] - in[stindx-1];
    if(fabs(diff) > MAX_DIFF){
      indxVec.push_back(stindx-1);
      if(diff < 0)
        count--;
      else
        count++;
      countVec.push_back(count);
    }
  }
  for(unsigned int i = 0; i < indxVec.size(); i++){
    int indx1 = indxVec[i];
    int indx2;
    count = countVec[i];
    if(i != (indxVec.size()-1))
      indx2 = indxVec[i+1]; 
    else
      indx2 = -1;
    for(int j = indx1; j > indx2; j--){
      in[j] += 2*PI*count;
    }
  }
  indxVec.clear();
  countVec.clear();
  int endindx = sz/2;
  count = 0;
  while(endindx < sz-1){
    endindx++;
    double diff = in[endindx] - in[endindx-1];
    if(fabs(diff) > MAX_DIFF){
      indxVec.push_back(endindx);
      if(diff > 0)
        count--;
      else
        count++;
      countVec.push_back(count);
    }
  }
  for(unsigned int i = 0; i < indxVec.size(); i++){
    int indx1 = indxVec[i];
    int indx2;
    count = countVec[i];
    if(i != (indxVec.size()-1))
      indx2 = indxVec[i+1]; 
    else
      indx2 = sz;
    for(int j = indx1; j < indx2; j++){
      in[j] += 2*PI*count;
    }
  }
}
   
bool CheckSignChange(std::vector<double>& in,int indx1,int indx2)
{
  if(in[indx1] > 0){
    if(in[indx2] < 0)
      return true;
    else 
      return false;
  }
  else{
    if(in[indx2] < 0)
      return false;
    else 
      return true;
  }
}

void getPhaseLimits(std::vector<double>& in,int& stindx,int& endindx,int sz)
{
  const double MAX_DIFF = PI/2.0;
  stindx = sz/2;
  endindx = sz/2;
  double diff = fabs(in[stindx] - in[stindx-1]);
  while(stindx > 0 && diff < MAX_DIFF){
    stindx--;
    diff = fabs(in[stindx] - in[stindx-1]);
  }
  diff = fabs(in[endindx] - in[endindx+1]);
  while(endindx < sz && diff < MAX_DIFF){
    endindx++;
    diff = fabs(in[endindx] - in[endindx+1]);
  }
  if( (endindx - sz/6) > (stindx + sz/6) ){
    stindx += sz/6;
    endindx -= sz/6;
  }
  else if( (endindx - sz/8) > (stindx + sz/8) ){
    stindx += sz/8;
    endindx -= sz/8;
  }
  else if( (endindx - sz/12) > (stindx + sz/12) ){
    stindx += sz/12;
    endindx -= sz/12;
  }
  else if( (endindx - 2) > (stindx + 2)){
    stindx+=2;
    endindx-=2;
  }
}

int ComputeOutN_2D(int nD,int stride)
{
  double val = (double) nD/stride;
  val = ceil(val);
  int n_out = (int) val;
  return n_out;
}

}

