
#pragma once

#include <complex> 
#include <string> 
#include <vector> 
#include <map> 
#include <fstream>
#include <cassert>
#include "pwutils/report/basedata.hpp"
#include "pwutils/report/basetrack.hpp"
#include "pwutils/report/reporthelper.h"
#include "pwutils/pwmath.hpp"

namespace dat{

void streamToDat(std::ofstream& os,const pw::metadataMap& str_map);

template<typename T1,typename T2>
void writeDatRow1D(std::ofstream& os,T1 x,T2 y)
{
    os << x << " " << y << std::endl;
}

template<typename T1,typename T2>
void writeDat1D(std::ofstream& os,const std::vector<T1>& x,const std::vector<T2>& y) 
{
    assert (x.size() == y.size());
	for(unsigned int i = 0; i < x.size(); i++){
		os << x[i] << " " << y[i] << std::endl;
	}
}

template<typename T1,typename T2>
void writeDat1D(std::ofstream& os,const std::vector<T1>& x,\
    const std::vector<std::complex<T2>>& y)
{
    assert (x.size() == y.size());
	for(unsigned int i = 0; i < x.size(); i++)
		os << x[i] << " " <<  y[i].real() << " " <<  y[i].imag() << std::endl;
}

template<typename T1,typename T2>
void writeDat2D_xy(std::ofstream& os,const std::vector<T1>& x,unsigned int strideX,\
        const std::vector<T2>& y,unsigned int strideY)
{
    unsigned int nx = pw::intceil(x.size(),strideX);
    unsigned int ny = pw::intceil(y.size(),strideY);
    os << nx << " " << ny << std::endl;
    for(unsigned int i = 0; i < x.size(); i+=strideX)
        os << x[i] << std::endl; 
    for(unsigned int j = 0; j < y.size(); j+=strideY)
        os << y[j] << std::endl;
}

template<typename T1,typename T2,typename T3>
void writeDat2D(std::ofstream& os,const std::vector<T1>& x,unsigned int strideX,\
        const std::vector<T2>& y,unsigned int strideY,\
		const std::vector<T3>& z) 
{
    assert (x.size()*y.size() == z.size());
    writeDat2D_xy(os,x,strideX,y,strideY);
    for(unsigned int i = 0; i < x.size(); i+=strideX){
    	for(unsigned int j = 0; j < y.size(); j+=strideY){
       	    os << z[i*y.size()+j] << " ";
    	}
    	os << std::endl;
    }
}

// dat 2D format (real,imaginary,real,imaginary)
template<typename T1,typename T2,typename T3>
void writeDat2D(std::ofstream& os,const std::vector<T1>& x,unsigned int strideX,\
        const std::vector<T2>& y,unsigned int strideY,\
    const std::vector<std::complex<T3>>& z)
{
    assert (x.size()*y.size() == z.size());
    writeDat2D_xy(os,x,strideX,y,strideY);
    for(unsigned int i = 0; i < x.size(); i+=strideX){
    	for(unsigned int j = 0; j < y.size(); j+=strideY){
            os << z[i*y.size()+j].real() << " " << z[i*y.size()+j].imag() << " ";
    	}
    	os << std::endl;
    }
}

template<typename T1,typename T2,typename T3>
void writePowerDat2D(std::ofstream& os,const std::vector<T1>& x,unsigned int strideX,\
        const std::vector<T2>& y,unsigned int strideY,\
        const std::vector<std::complex<T3>>& z)
{
    assert (x.size()*y.size() == z.size());

    // 2D DAT format prints the integer sizes of the x array and y array on the same line
    writeDat2D_xy(os,x,strideX,y,strideY);
    for(unsigned int i = 0; i < x.size(); i+=strideX){
        for(unsigned int j = 0; j < y.size(); j+=strideY){
            os << pow(abs(z[i*y.size()+j]),2) << " ";
        }
    	os << std::endl;
    }
}

template<typename T1,typename T2,typename T3>
void writePhaseDat2D(std::ofstream& os,const std::vector<T1>& x,unsigned int strideX,\
        const std::vector<T2>& y,unsigned int strideY,\
        const std::vector<std::complex<T3>>& z)
{
    assert (x.size()*y.size() == z.size());
    writeDat2D_xy(os,x,strideX,y,strideY);
    for(unsigned int i = 0; i < x.size(); i+=strideX){
        for(unsigned int j = 0; j < y.size(); j+=strideY){
            os << arg(z[i]) << " ";
        }
    	os << std::endl;
    }
}

template<typename T1,typename T2>
void writePowerDat1D(std::ofstream& os,const std::vector<T1>& x,\
    const std::vector<std::complex<T2>>& y)
{
  assert (x.size() == y.size());
	for(unsigned int i = 0; i < x.size(); i++){
		os << x[i] << " " << pow(abs(y[i]),2) << std::endl;
	}
}

template<typename T1,typename T2>
void writePhaseDat1D(std::ofstream& os,const std::vector<T1>& x,
        const std::vector<std::complex<T2>>& y)
{
    assert (x.size() == y.size());
    std::vector<double> phaseVec(y.size());
    for(unsigned int i=0; i < y.size(); i++)
        phaseVec[i] = arg(y[i]);
    pw::AdjustPhase(phaseVec,phaseVec.size());
    writeDat1D(os,x,phaseVec);
}


template<typename T>
void writeRowVec(std::ofstream& os,const std::vector<T>& x) {
	for(unsigned int i = 0; i < x.size(); i++){
		os << x[i] << " "; 
    }
	os << std::endl;
}

template<typename T>
void writeColVec(std::ofstream& os,const std::vector<T>& x)
{
	for(unsigned int i = 0; i < x.size(); i++){
		os << x[i] << std::endl; 
    }
}

template<typename T>
void writeColVec(std::ofstream& os,const std::vector<std::complex<T>>& x)
{
	for(unsigned int i = 0; i < x.size(); i++){
		os << x[i].real() << " " <<  x[i].imag() << std::endl;
	}
}

template<class T1,class T2>
class ReportData1D : public pw::ReportDataBase1D<T1,T2>
{
    public:
        ReportData1D(const std::string& name,
            const std::vector<T1>& x, 
            const std::vector<T2>& y) :
                pw::ReportDataBase1D<T1,T2>(name,x,y) {
                    pw::ReportBase::setFileExtension("dat");
                    pw::ReportBase::setFileSignature(pw::FileSignature::DAT);
                }
        ~ReportData1D() {};
    private:
		void reportMetadata(std::ofstream& os) const {streamToDat(os,this->getMetadata());}
		void reportData(std::ofstream& os) const; 
		void reportImplement(std::ofstream& os) const {
            if(pw::ReportBase::metadataOn())
                reportMetadata(os);
            reportData(os);
        }
};

template<class T1,class T2>
void ReportData1D<T1,T2>::reportData(std::ofstream& os) const
{
    writeDat1D(os,this->getX(),this->getY());
}

template<class T1,class T2>
class ReportComplexData1D : public pw::ReportComplexDataBase1D<T1,T2>
{
	public :
        ReportComplexData1D(const std::string& name,
            const std::vector<T1>& x,
            const std::vector<std::complex<T2>>& y) :
                pw::ReportComplexDataBase1D<T1,T2>(name,x,y) {
                    pw::ReportBase::setFileExtension("dat");
                    pw::ReportBase::setFileSignature(pw::FileSignature::DAT);
                }
		~ReportComplexData1D() {}
    private:
		void reportMetadata(std::ofstream& os) const {streamToDat(os,this->getMetadata());}
		void reportData(std::ofstream& os) const; 
		void reportImplement(std::ofstream& os) const {
            if(pw::ReportBase::metadataOn())
                reportMetadata(os);
            reportData(os);
        }

};

template<class T1,class T2>
void ReportComplexData1D<T1,T2>::reportData(std::ofstream& os) const
{
    if(this->getPhase())
		writePhaseDat1D(os,this->getX(),this->getY());
    else if(this->getPower()){
		writePowerDat1D(os,this->getX(),this->getY());
	} else
        writeDat1D(os,this->getX(),this->getY());
}

template<class T1,class T2,class T3>
class ReportData2D : public pw::ReportDataBase2D<T1,T2,T3>
{
    public:
        ReportData2D(const std::string& name,
            const std::vector<T1>& x, 
            const std::vector<T2>& y, 
            const std::vector<T3>& z) :
                pw::ReportDataBase2D<T1,T2,T3>(name,x,y,z) {
                pw::ReportBase::setFileExtension("dat");
                pw::ReportBase::setFileSignature(pw::FileSignature::DAT);
            }
        ~ReportData2D() {};
    private:
        void reportMetadata(std::ofstream& os) const {streamToDat(os,this->getMetadata());}
        void reportData(std::ofstream& os) const; 
		void reportImplement(std::ofstream& os) const {
            if(pw::ReportBase::metadataOn())
                reportMetadata(os);
            reportData(os);
        }

};

template<class T1,class T2,class T3>
void ReportData2D<T1,T2,T3>::reportData(std::ofstream& os) const
{
    writeDat2D(os,this->getX(),this->getStrideX(),\
            this->getY(),this->getStrideY(),this->getZ());
}


template<class T1,class T2,class T3>
class ReportComplexData2D : public pw::ReportComplexDataBase2D<T1,T2,T3>
{
    public :
        ReportComplexData2D(const std::string& name,
            const std::vector<T1>& x,
            const std::vector<T2>& y,
            const std::vector<std::complex<T3>>& z) :
            pw::ReportComplexDataBase2D<T1,T2,T3>(name,x,y,z) {
                 pw::ReportBase::setFileExtension("dat");
                 pw::ReportBase::setFileSignature(pw::FileSignature::DAT);
            }
        ~ReportComplexData2D() {}
    private:
    	void reportMetadata(std::ofstream& os) const {streamToDat(os,this->getMetadata());}
    	void reportData(std::ofstream& os) const; 
        void reportImplement(std::ofstream& os) const {
            if(pw::ReportBase::metadataOn())
                reportMetadata(os);
            reportData(os);
        }
};

template<class T1,class T2,class T3>
void ReportComplexData2D<T1,T2,T3>::reportData(std::ofstream& os) const
{
    if(this->getPhase()){
    		writePhaseDat2D(os,this->getX(),this->getStrideX(),\
    		        this->getY(),this->getStrideY(),this->getZ());
    }
    else if(this->getPower()){
    		writePowerDat2D(os,this->getX(),this->getStrideX(),\
    		        this->getY(),this->getStrideY(),this->getZ());
	} else{
        writeDat2D(os,this->getX(),this->getStrideX(),\
                this->getY(),this->getStrideY(),this->getZ());
    }
}


template<class T>
class TrackData : public pw::TrackDataBase<T>
{
    public:
        TrackData(const std::string& name,
            pw::TrackType ttype,
            const std::vector<T>& data) :
                pw::TrackDataBase<T>(name,ttype,data) {
                    pw::ReportBase::setFileExtension("dat");
                    pw::ReportBase::setFileSignature(pw::FileSignature::DAT);
                }
    private:
		void reportMetadata(std::ofstream& os) const {streamToDat(os,this->getMetadata());}
		void reportData(std::ofstream& os) const; 
        void reportImplement(std::ofstream& os) const {
            if(pw::ReportBase::metadataOn())
                reportMetadata(os);
            reportData(os);
        }
};

template<class T>
void TrackData<T>::reportData(std::ofstream& os) const 
{
    writeDat1D(os,this->getX(),this->getY());
}

template<class T>
class TrackComplexData : public pw::TrackComplexDataBase<T>
{
    public:
        TrackComplexData(const std::string& name,
            pw::TrackType ttype,
            const std::vector<std::complex<T>>& data,
            pw::ComplexOp cmplxop = pw::ComplexOp::None) : 
                pw::TrackComplexDataBase<T>(name,ttype,data) {
                    pw::ReportBase::setFileExtension("dat");
                    pw::ReportBase::setFileSignature(pw::FileSignature::DAT);
                }
    private:
		void reportMetadata(std::ofstream& os) const {streamToDat(os,this->getMetadata());}
		void reportData(std::ofstream& os) const; 
        void reportImplement(std::ofstream& os) const {
            if(pw::ReportBase::metadataOn())
                reportMetadata(os);
            reportData(os);
        }
};

template<class T>
void TrackComplexData<T>::reportData(std::ofstream& os) const 
{

	if(pw::TrackComplexDataBase<T>::getComplexOp() == pw::ComplexOp::None)
        writeDat1D(os,this->getX(),this->getY());
    else
        writeDat1D(os,this->getX(),this->getOpY());
}




}



