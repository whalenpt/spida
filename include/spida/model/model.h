
#ifndef MODEL_H_
#define MODEL_H_

#include "spida/helper/constants.h"
#include <vector>
#include <pwutils/pwthreads.h>

namespace spida{

//enum class 

enum class Dimension{D1,D2,D3};

class ModelCV
{
  public:
      ModelCV(int c_nthreads = 1) : thmgt(c_nthreads) {}
      virtual ~ModelCV() {};
      virtual void nonLinResponse(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) = 0;
      virtual const std::vector<dcmplx>& linOp() = 0;
      unsigned int numThreads() {return thmgt.getNumThreads();}
      void setNumThreads(unsigned int val) {thmgt.setNumThreads(val);}
      pw::ThreadManager& threadManager() {return thmgt;}
      virtual Dimension dimension() = 0;
  private:
      pw::ThreadManager thmgt;

};

}

#endif




