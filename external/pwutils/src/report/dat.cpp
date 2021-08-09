
#include "pwutils/report/dat.hpp"
#include "pwutils/report/reporthelper.h"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <ctime>
#include <iostream>


namespace dat{

using pw::dcmplx;
using pw::REPORT_PADING;

void streamToDat(std::ofstream& os,const pw::metadataMap& str_map) 
{
    pw::metadataMap::const_iterator it;
    for(it = str_map.begin(); it!= str_map.end(); it++){
        os << "#" << (*it).first + ": " << (*it).second << std::endl;
    }
}



/*
void ReportComplexData2D::reportData(std::ofstream& os) const
{
    if(this->getPhase())
		writePhaseDat2D(os,this->getX(),this->getY(),this->getZ(),\
		        this->precision());
    else if(this->getPower()){
		writePowerDat2D(os,this->getX(),this->getY(),this->getZ(),\
		        this->precision());
	} else
        writeDat2D(os,this->getX(),this->getY(),this->getZ(),this->precision());
}


void ReportRealData2D::reportData(std::ofstream& os) const
{
    writeDat2D(os,this->getX(),this->getY(),this->getZ(),this->precision());
}

void ReportRealTrackerMax::reportTracker(std::ofstream& os,double t) const
{
    writeDatRow1D(os,t,*std::max_element(m_v.begin(),m_v.end()),this->precision());
}

void ReportComplexTrackerMax::reportTracker(std::ofstream& os,double t) const
{
    double max_val = 0;
    for(auto val : m_v)
        max_val = std::max(abs(val),max_val);
    writeDatRow1D(os,t,max_val,this->precision());
}
*/


/*
void writeDat2D(std::ofstream& os,const std::vector<double>& x,const std::vector<double>& y,
        const std::vector<double>& z,int precision) 
{
    assert (x.size()*y.size() == z.size());
    writeRowVec(os,x,precision);
    writeRowVec(os,y,precision);
    int nx = x.size(); 
    int ny = y.size();
	for(unsigned int i = 0; i < nx; i++){
        for(unsigned int j = 0; j < ny; j++){
            os << std::scientific << std::setprecision(precision) \
                << std::setw(precision+REPORT_PADING) << z[nx*i+j]; 
        }
        os << std::endl;
	}
}

void writeDat2D(std::ofstream& os,const std::vector<double>& x,const std::vector<double>& y,
        const std::vector<dcmplx>& z,int precision) 
{
    assert (x.size()*y.size() == z.size());
    writeRowVec(os,x,precision);
    writeRowVec(os,y,precision);
    int nx = x.size(); 
    int ny = y.size();
	for(unsigned int i = 0; i < nx; i++){
        for(unsigned int j = 0; j < ny; j++){
            os << std::scientific << std::setprecision(precision) \
               << std::setw(precision+REPORT_PADING) << z[nx*i+j].real();
            os << std::scientific << std::setprecision(precision) \
               << std::setw(precision+REPORT_PADING) << z[nx*i+j].imag();
        }
        os << std::endl;
	}
}

void writePowerDat2D(std::ofstream& os,const std::vector<double>& x,const std::vector<double>& y,
        const std::vector<dcmplx>& z,int precision) 
{
    assert (x.size()*y.size() == z.size());
    writeRowVec(os,x,precision);
    writeRowVec(os,y,precision);
    int nx = x.size(); 
    int ny = y.size();
	for(unsigned int i = 0; i < nx; i++){
        for(unsigned int j = 0; j < ny; j++){
            os << std::scientific << std::setprecision(precision) \
               << std::setw(precision+REPORT_PADING) << pow(abs(z[nx*i+j]),2);
        }
        os << std::endl;
	}
}

void writePhaseDat2D(std::ofstream& os,const std::vector<double>& x,const std::vector<double>& y,
        const std::vector<dcmplx>& z,int precision) 
{
    assert (x.size()*y.size() == z.size());
    std::vector<double> phaseVec(z.size());
    for(unsigned int i=0; i < z.size(); i++)
        phaseVec[i] = arg(z[i]);
    pw::AdjustPhase(phaseVec,z.size());
    writeDat2D(os,x,y,phaseVec,precision);
}
*/

}

