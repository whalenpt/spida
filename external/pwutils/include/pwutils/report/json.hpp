// json.hpp
#pragma once

#include "pwutils/report/basedata.hpp"
#include "pwutils/report/basetrack.hpp"
#include "pwutils/report/reporthelper.h"
#include "pwutils/pwmath.hpp"
#include <complex> 
#include <string> 
#include <vector> 
#include <memory>
#include <fstream>
#include <iostream>
#include <map> 

namespace json{

void streamToJSON(std::ofstream& os,const pw::metadataMap& str_map);
void writeJSONLabel(std::ofstream& os,const std::string& label,const std::string& indent="\t");
void writeJSONValue(std::ofstream& os,const std::string& value,
		const std::string& indent="",bool end_value=false);
void mapToJSON(std::ofstream& os,const pw::metadataMap& str_map,bool end_value=false);

template<typename T>
void writeJSONVector(std::ofstream& os,const std::string& label,const std::vector<T>& v,\
        unsigned int stride=1,const std::string& indent="\t",bool end_value=false)
{
    writeJSONLabel(os,label,indent);
    os << "[";
    unsigned int i;
    for(i = 0; i < v.size()-stride; i+=stride){
    	os << v[i] << ", ";
    }
    std::string end_string = end_value ? "]" : "],";
    os << v[i] << end_string << std::endl;
}

template<typename T>
void writeJSONVector(std::ofstream& os,const std::string& label,const std::vector<std::complex<T>>& v,\
        unsigned int stride,const std::string& indent,bool end_value)
{
    writeJSONLabel(os,label,indent);
    os << "[";
    unsigned int i;
    for(i = 0; i < v.size()-stride; i+=stride){
    	os << v[i].real() << ", " << v[i].imag() << ", ";
    }
    std::string end_string = end_value ? "]" : "],";
    os << v[i].real() << ", " << v[i].imag() << end_string << std::endl;

}

template<typename T>
void writeJSONPowerVector(std::ofstream& os,const std::string& label,const std::vector<std::complex<T>>& v,\
        unsigned int stride=1,const std::string& indent="\t",bool end_value = false)
{
    writeJSONLabel(os,label,indent);
	os << "[";
    unsigned int i;
	for(i=0; i < v.size()-stride; i+=stride){
        os << pow(abs(v[i]),2) << ", ";
	}
    std::string end_string = end_value ? "]" : "],";
	os << pow(abs(v[i]),2) << end_string << std::endl;
}



template<typename T>
void writeJSONPhaseVector(std::ofstream& os,const std::string& label,const std::vector<std::complex<T>>& v,\
        unsigned int stride=1,const std::string& indent="\t",bool end_value=false)
{
		std::vector<T> phaseVec(v.size(),0.0);
        unsigned int count = 0;
        unsigned int i;
		for(i = 0; i < v.size(); i+=stride){
            phaseVec[count] = arg(v[i]);
            count++;
        }
        phaseVec[count] = arg(v[i]);
        count++;
        phaseVec.resize(count);
//		pw::AdjustPhase(phaseVec,count);
		writeJSONVector<T>(os,label,phaseVec,1,indent,end_value);
}

template<typename T1,typename T2>
void writeJSON2D_xy(std::ofstream& os,\
        const std::vector<T1>& x,unsigned int strideX,\
        const std::vector<T2>& y,unsigned int strideY,\
        const std::string& indent="\t")
{
    writeJSONVector(os,pw::XLABEL,x,strideX,indent);
    writeJSONVector(os,pw::YLABEL,y,strideY,indent);
}

template<typename T1,typename T2,typename T3>
void writeJSON2D(std::ofstream& os,\
        const std::vector<T1>& x,unsigned int strideX,\
        const std::vector<T2>& y,unsigned int strideY,\
        const std::vector<T3>& z,\
        const std::string& indent="\t")
{
    assert (x.size()*y.size() == z.size());
    writeJSON2D_xy(os,x,strideX,y,strideY,indent);
    writeJSONLabel(os,pw::ZLABEL,indent);
	os << "[";
    unsigned int i = 0;
	for(i = 0; i < x.size()-strideX; i+=strideX){
	    for(unsigned int j=0; j < y.size(); j+=strideY){
            os << z[i*y.size()+j] << ", ";
        }
	}
    unsigned int j;
	for(j = 0; j < y.size()-strideY; j++)
        os << z[i*y.size()+j] << ", ";
    os << z[i*y.size()+j];
    os << "]" << std::endl;
}

template<typename T1,typename T2,typename T3>
void writeJSON2D(std::ofstream& os,\
        const std::vector<T1>& x,unsigned int strideX,\
        const std::vector<T2>& y,unsigned int strideY,\
        const std::vector<std::complex<T3>>& z,\
        const std::string& indent="\t")
{
    assert (x.size()*y.size() == z.size());
    writeJSON2D_xy(os,x,strideX,y,strideY,indent);
    writeJSONLabel(os,pw::ZLABEL,indent);
	os << "[";
    unsigned int i = 0;
	for(i = 0; i < x.size()-strideX; i+=strideX){
	    for(unsigned int j=0; j < y.size(); j+=strideY){
            os << z[i*y.size()+j].real() << ", " << z[i*y.size()+j].imag() << ", ";
        }
	}
    unsigned int j;
	for(j = 0; j < y.size()-strideY; j++)
        os << z[i*y.size()+j].real() << ", " << z[i*y.size()+j].imag() << ", ";
    os << z[i*y.size()+j].real() << ", " << z[i*y.size()+j].imag();
    os << "]" << std::endl;
}

template<typename T1,typename T2,typename T3>
void writeJSONPower2D(std::ofstream& os,\
        const std::vector<T1>& x,unsigned int strideX,\
        const std::vector<T2>& y,unsigned int strideY,\
        const std::vector<std::complex<T3>>& z,\
        const std::string& indent="\t")
{
    assert (x.size()*y.size() == z.size());
    writeJSON2D_xy(os,x,strideX,y,strideY,indent);
    writeJSONLabel(os,pw::ZLABEL,indent);
	os << "[";
    unsigned int i = 0;
	for(i = 0; i < x.size()-strideX; i+=strideX){
	    for(unsigned int j=0; j < y.size(); j+=strideY){
            os << pow(abs(z[i*y.size()+j]),2) << ", ";
        }
	}
    unsigned int j;
	for(j = 0; j < y.size()-strideY; j++)
        os << pow(abs(z[i*y.size()+j]),2) << ", ";
    os << pow(abs(z[i*y.size()+j]),2);
    os << "]" << std::endl;
}

template<typename T1,typename T2,typename T3>
void writeJSONPhase2D(std::ofstream& os,\
        const std::vector<T1>& x,unsigned int strideX,\
        const std::vector<T2>& y,unsigned int strideY,\
        const std::vector<std::complex<T3>>& z,\
        const std::string& indent="\t")
{
    assert (x.size()*y.size() == z.size());
    writeJSON2D_xy(os,x,strideX,y,strideY,indent);
    writeJSONLabel(os,pw::ZLABEL,indent);
	os << "[";
    unsigned int i = 0;
	for(i = 0; i < x.size()-strideX; i+=strideX){
	    for(unsigned int j=0; j < y.size(); j+=strideY){
            os << arg(z[i*y.size()+j]) << ", ";
        }
	}
    unsigned int j;
	for(j = 0; j < y.size()-strideY; j++)
        os << arg(z[i*y.size()+j]) << ", ";
    os << arg(z[i*y.size()+j],2);
    os << "]" << std::endl;
}

template<class T1,class T2>
class ReportData1D : public pw::ReportDataBase1D<T1,T2>
{
    public:
        ReportData1D(const std::string& name,
            const std::vector<T1>& x, 
            const std::vector<T2>& y) :
                pw::ReportDataBase1D<T1,T2>(name,x,y) {
                    pw::ReportBase::setFileExtension("json");
                    pw::ReportBase::setFileSignature(pw::FileSignature::JSON);
                }
        ~ReportData1D() {};
    private:
		void reportMetadata(std::ofstream& os) const {streamToJSON(os,this->getMetadata());}
		void reportData(std::ofstream& os) const; 
        void reportImplement(std::ofstream& os) const {
            os << "{" << std::endl;
            if(pw::ReportBase::metadataOn())
                reportMetadata(os);
            reportData(os);
            os << "}" << std::endl;
        }

};

template<class T1,class T2>
void ReportData1D<T1,T2>::reportData(std::ofstream& os) const 
{
	writeJSONVector(os,pw::XLABEL,this->getX(),1,"\t",false);
    writeJSONVector(os,pw::YLABEL,this->getY(),1,"\t",true);
}

template<class T1,class T2>
class ReportComplexData1D : public pw::ReportComplexDataBase1D<T1,T2>
{
	public :
        ReportComplexData1D(const std::string& name,
            const std::vector<T1>& x,
            const std::vector<std::complex<T2>>& y) :
                pw::ReportComplexDataBase1D<T1,T2>(name,x,y) {
                    pw::ReportBase::setFileExtension("json");
                    pw::ReportBase::setFileSignature(pw::FileSignature::JSON);
                }
		~ReportComplexData1D() {}
    private:
		void reportMetadata(std::ofstream& os) const {streamToJSON(os,this->getMetadata());}
		void reportData(std::ofstream& os) const; 
        void reportImplement(std::ofstream& os) const {
            os << "{" << std::endl;
            if(pw::ReportBase::metadataOn())
                reportMetadata(os);
            reportData(os);
            os << "}" << std::endl;
        }
};

template<class T1,class T2>
void ReportComplexData1D<T1,T2>::reportData(std::ofstream& os) const
{
    if(!this->getPhase() && !this->getPower()){
        writeJSONLabel(os,"dtype");
        writeJSONValue(os,"complex");
    }
    writeJSONVector(os,pw::XLABEL,this->getX(),1,"\t",false);
	if(this->getPhase())
		writeJSONPhaseVector(os,pw::YLABEL,this->getY(),1,"\t",true);
    else if(this->getPower()){
	    writeJSONPowerVector(os,pw::YLABEL,this->getY(),1,"\t",true);
	} else{
        writeJSONVector(os,pw::YLABEL,this->getY(),1,"\t",true);
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
                    pw::ReportBase::setFileExtension("json");
                    pw::ReportBase::setFileSignature(pw::FileSignature::JSON);
                }
        ~TrackData() {};
    private:
		void reportMetadata(std::ofstream& os) const {streamToJSON(os,this->getMetadata());}
		void reportData(std::ofstream& os) const; 
        void reportImplement(std::ofstream& os) const {
            os << "{" << std::endl;
            if(pw::ReportBase::metadataOn())
                reportMetadata(os);
            reportData(os);
            os << "}" << std::endl;
        }

};

template<class T>
void TrackData<T>::reportData(std::ofstream& os) const 
{
	writeJSONVector(os,pw::XLABEL,this->getX(),"\t",false);
    writeJSONVector(os,pw::YLABEL,this->getY(),"\t",true);
}

template<class T>
class TrackComplexData : public pw::TrackComplexDataBase<T>
{
    public:
        TrackComplexData(const std::string& name,
            pw::TrackType ttype,
            const std::vector<std::complex<T>>& data, 
            pw::ComplexOp cmplxop = pw::ComplexOp::None) : 
                pw::TrackComplexDataBase<T>(name,ttype,data,cmplxop) {
                    pw::ReportBase::setFileExtension("json");
                    pw::ReportBase::setFileSignature(pw::FileSignature::JSON);
                }
    private:
		void reportMetadata(std::ofstream& os) const {streamToJSON(os,this->getMetadata());}
		void reportData(std::ofstream& os) const; 
        void reportImplement(std::ofstream& os) const {
            os << "{" << std::endl;
            if(pw::ReportBase::metadataOn())
                reportMetadata(os);
            reportData(os);
            os << "}" << std::endl;
        }

};

template<class T>
void TrackComplexData<T>::reportData(std::ofstream& os) const 
{
	writeJSONVector(os,pw::XLABEL,this->getX(),1,"\t",false);
	if(pw::TrackComplexDataBase<T>::getComplexOp() == pw::ComplexOp::None)
        writeJSONVector(os,pw::YLABEL,this->getY(),1,"\t",true);
    else
        writeJSONVector(os,pw::YLABEL,this->getOpY(),1,"\t",true);
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
                    pw::ReportBase::setFileExtension("json");
                    pw::ReportBase::setFileSignature(pw::FileSignature::JSON);
                }
        ~ReportData2D() {};
    private:
		void reportMetadata(std::ofstream& os) const {streamToJSON(os,this->getMetadata());}
		void reportData(std::ofstream& os) const; 
        void reportImplement(std::ofstream& os) const {
            os << "{" << std::endl;
            if(pw::ReportBase::metadataOn())
                reportMetadata(os);
            reportData(os);
            os << "}" << std::endl;
        }

};

template<class T1,class T2,class T3>
void ReportData2D<T1,T2,T3>::reportData(std::ofstream& os) const 
{
    writeJSON2D(os,this->getX(),this->getStrideX(),this->getY(),this->getStrideY(),\
            this->getZ(),"\t");
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
                    pw::ReportBase::setFileExtension("json");
                    pw::ReportBase::setFileSignature(pw::FileSignature::JSON);
                }
		~ReportComplexData2D() {}
    private:
		void reportMetadata(std::ofstream& os) const {streamToJSON(os,this->getMetadata());}
		void reportData(std::ofstream& os) const; 
        void reportImplement(std::ofstream& os) const {
            os << "{" << std::endl;
            if(pw::ReportBase::metadataOn())
                reportMetadata(os);
            reportData(os);
            os << "}" << std::endl;
        }
};

template<class T1,class T2,class T3>
void ReportComplexData2D<T1,T2,T3>::reportData(std::ofstream& os) const
{
    if(!this->getPhase() && !this->getPower()){
        writeJSONLabel(os,"dtype");
        writeJSONValue(os,"complex");
    }

    if(this->getPhase()){
        writeJSONPower2D(os,this->getX(),this->getStrideX(),\
            this->getY(),this->getStrideY(),this->getZ(),"\t");
    }
    else if(this->getPower()){
        writeJSONPower2D(os,this->getX(),this->getStrideX(),this->getY(),this->getStrideY(),\
            this->getZ(),"\t");
	} else{
        writeJSON2D(os,this->getX(),this->getStrideX(),this->getY(),this->getStrideY(),\
            this->getZ(),"\t");
    }
}




}



