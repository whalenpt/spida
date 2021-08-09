
#ifndef MODEL_H_
#define MODEL_H_

#include "spida/constants.h"
#include <vector>
#include <pwutils/pwthreads.h>

namespace spida{

class ModelDC
{
  public:
      ModelDC(int c_nthreads = 1) : thmgt(c_nthreads) {}
      virtual ~ModelDC() {};
      virtual void nonLinResponse(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) = 0;
      virtual const std::vector<dcmplx>& linOperator() = 0;
      unsigned int numThreads() {return thmgt.getNumThreads();}
      void setNumThreads(unsigned int val) {thmgt.setNumThreads(val);}
      pw::ThreadManager& threadManager() {return thmgt;}
      virtual int numDimensions() const = 0;
      virtual int specDimSize() const = 0;
  private:
      pw::ThreadManager thmgt;

};

}

#endif




