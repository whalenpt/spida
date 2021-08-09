
#ifndef SPIDA_X_H_
#define SPIDA_X_H_

#include <vector>

namespace spida{

class GridX;
class TransformX;

// Interface class
class SpidaX
{
    public:
        SpidaX(const GridX& grid) {}

        // No copy constructors or assignments please 
        SpidaX() = delete;
        SpidaX(const SpidaX&)=delete;
        SpidaX& operator=(const SpidaX&)=delete;

        virtual ~SpidaX() {};
        virtual void X_To_SX(const std::vector<double>& in,std::vector<double>& out) = 0; 
        virtual void SX_To_X(const std::vector<double>& in,std::vector<double>& out) = 0;
        virtual void dX(const std::vector<double>& in,std::vector<double>& out,int n) = 0; 
        virtual void dX(const std::vector<double>& in,std::vector<double>& out) = 0; 
        virtual const GridX& getGridX() = 0;
        virtual const TransformX& getTransformX() = 0;
    private:
        virtual void dSX(const std::vector<double>& in,std::vector<double>& out,int n) = 0; 
        virtual void dSX(const std::vector<double>& in,std::vector<double>& out) = 0; 

        virtual void dSX(std::vector<double>& a,int n) = 0; 
        virtual void dSX(std::vector<double>& a) = 0; 
};

}

#endif


