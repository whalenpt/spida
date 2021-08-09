
#ifndef DAT_HPP_
#define DAT_HPP_ 

#include <complex> 
#include <string> 
#include <vector> 
#include <map> 
#include <fstream>
#include <cassert>
#include "pwutils/report/basedata.hpp"
#include "pwutils/report/basetrack.hpp"
#include "pwutils/report/reporthelper.h"

namespace dat{

using pw::dcmplx;

void streamToDat(std::ofstream& os,const pw::metadataMap& str_map);

template<typename T1,typename T2>
void writeDatRow1D(std::ofstream& os,T1 x,T2 y,int precision=pw::REPORT_PRECISION)
{
    os << std::scientific << std::setprecision(precision) << std::setw(precision+pw::REPORT_PADING) << x \
       << std::scientific << std::setprecision(precision) << std::setw(precision+pw::REPORT_PADING) << y << std::endl;
}

template<typename T1,typename T2>
void writeDat1D(std::ofstream& os,const std::vector<T1>& x,const std::vector<T2>& y,int precision=pw::REPORT_PRECISION) 
{
    assert (x.size() == y.size());
	for(unsigned int i = 0; i < x.size(); i++){
		os << std::scientific << std::setprecision(precision) << std::setw(precision+pw::REPORT_PADING) << x[i] \
		   << std::scientific << std::setprecision(precision) << std::setw(precision+pw::REPORT_PADING) << y[i] << std::endl;
	}
}

template<typename T1>
void writeDat1D(std::ofstream& os,const std::vector<T1>& x,\
    const std::vector<dcmplx>& y,int precision=pw::REPORT_PRECISION)
{
    assert (x.size() == y.size());
	for(unsigned int i = 0; i < x.size(); i++){
		os << std::scientific << std::setprecision(precision) << std::setw(precision+pw::REPORT_PADING) << x[i] \
		   << std::scientific << std::setprecision(precision) << std::setw(precision+pw::REPORT_PADING) << y[i].real() \
		   << std::scientific << std::setprecision(precision) << std::setw(precision+pw::REPORT_PADING) << y[i].imag()
		   << std::endl;
	}
}

template<typename T1>
void writePowerDat1D(std::ofstream& os,const std::vector<T1>& x,\
    const std::vector<dcmplx>& y,int precision=pw::REPORT_PRECISION)
{
    assert (x.size() == y.size());
	for(unsigned int i = 0; i < x.size(); i++){
		os << std::scientific << std::setprecision(precision) << std::setw(precision+pw::REPORT_PADING) << x[i] \
		   << std::scientific << std::setprecision(precision) << std::setw(precision+pw::REPORT_PADING) \
		   << pow(abs(y[i]),2) << std::endl;
	}
}

template<typename T1>
void writePhaseDat1D(std::ofstream& os,const std::vector<T1>& x,
        const std::vector<dcmplx>& y,int precision)
{
    assert (x.size() == y.size());
    std::vector<double> phaseVec(y.size());
    for(unsigned int i=0; i < y.size(); i++)
        phaseVec[i] = arg(y[i]);
    pw::AdjustPhase(phaseVec,phaseVec.size());
    writeDat1D(os,x,phaseVec,precision);
}


/*
void writeDat2D(std::ofstream& os,const std::vector<double>& x,const std::vector<double>& y,
    const std::vector<double>& z,int precision=pw::REPORT_PRECISION); 
void writeDat2D(std::ofstream& os,const std::vector<double>& x,const std::vector<double>& y,
    const std::vector<dcmplx>& z,int precision=pw::REPORT_PRECISION); 
void writePhaseDat2D(std::ofstream& os,const std::vector<double>& x,const std::vector<double>& y,
    const std::vector<dcmplx>& z,int precision=pw::REPORT_PRECISION); 
void writePowerDat2D(std::ofstream& os,const std::vector<double>& x,const std::vector<double>& y,
    const std::vector<dcmplx>& z,int precision=pw::REPORT_PRECISION); 
*/

template<typename T1>
void writeRowVec(std::ofstream& os,const std::vector<T1>& x,int precision=pw::REPORT_PRECISION) {
	for(unsigned int i = 0; i < x.size(); i++){
		os << std::scientific << std::setprecision(precision) << std::setw(precision+pw::REPORT_PADING) << x[i]; 
    }
	os << std::endl;
}

template<typename T1>
void writeColVec(std::ofstream& os,const std::vector<T1>& x,int precision=pw::REPORT_PRECISION)
{
	for(unsigned int i = 0; i < x.size(); i++){
		os << std::scientific << std::setprecision(precision) \
		   << std::setw(precision+pw::REPORT_PADING) << x[i] << std::endl; 
    }
}

template<>
void writeColVec(std::ofstream& os,const std::vector<dcmplx>& x,int precision)
{
	for(unsigned int i = 0; i < x.size(); i++){
		os << std::scientific << std::setprecision(precision) << std::setw(precision+pw::REPORT_PADING) << x[i].real() \
		   << std::scientific << std::setprecision(precision) << std::setw(precision+pw::REPORT_PADING) << x[i].imag() \
		   << std::endl;
	}
}

template<class T1,class T2>
class ReportData1D : public pw::ReportDataBase1D<T1,T2>
{
    public:
        ReportData1D(const std::string& name,
            const std::vector<T1>& x, 
            const std::vector<T2>& y, 
            const std::string& x_label = "x",
            const std::string& y_label = "y") : 
                pw::ReportDataBase1D<T1,T2>(name,x,y,x_label,y_label) {
                    pw::ReportBase::setFileExtension("dat");}
        ~ReportData1D() {};
    private:
		void reportMetadata(std::ofstream& os) const {streamToDat(os,this->metadata());}
		void reportData(std::ofstream& os) const; 
};

template<class T1,class T2>
void ReportData1D<T1,T2>::reportData(std::ofstream& os) const
{
    writeDat1D(os,this->getX(),this->getY(),this->precision());
}

template<class T1>
class ReportComplexData1D : public pw::ReportComplexDataBase1D<T1>
{
	public :
        ReportComplexData1D(const std::string& name,
            const std::vector<T1>& x,
            const std::vector<dcmplx>& y,
            std::string x_label="x",
            std::string y_label="y") : 
                pw::ReportComplexDataBase1D<T1>(name,x,y,x_label,y_label) {
                    pw::ReportBase::setFileExtension("dat");}
		~ReportComplexData1D() {}
    private:
		void reportMetadata(std::ofstream& os) const {streamToDat(os,this->metadata());}
		void reportData(std::ofstream& os) const; 
};

template<class T1>
void ReportComplexData1D<T1>::reportData(std::ofstream& os) const
{
    if(this->getPhase())
		writePhaseDat1D(os,this->getX(),this->getY(),this->precision());
    else if(this->getPower()){
		writePowerDat1D(os,this->getX(),this->getY(),this->precision());
	} else
        writeDat1D(os,this->getX(),this->getY(),this->precision());
}

template<class T>
class TrackData : public pw::TrackDataBase<T>
{
    public:
        TrackData(const std::string& name,
            pw::TrackType ttype,
            const std::vector<T>& data, 
            const std::string& x_label = "x",
            const std::string& y_label = "y") : 
                pw::TrackDataBase<T>(name,ttype,data,x_label,y_label) {
                    pw::ReportBase::setFileExtension("dat");}
    private:
		void reportMetadata(std::ofstream& os) const {streamToDat(os,this->metadata());}
		void reportData(std::ofstream& os) const; 
};

template<class T>
void TrackData<T>::reportData(std::ofstream& os) const 
{
    writeDat1D(os,this->getX(),this->getY(),this->precision());
}

class TrackComplexData : public pw::TrackComplexDataBase
{
    public:
        TrackComplexData(const std::string& name,
            pw::TrackType ttype,
            const std::vector<dcmplx>& data, 
            const std::string& x_label = "x",
            const std::string& y_label = "y",
            pw::ComplexOp cmplxop = pw::ComplexOp::None) : 
                pw::TrackComplexDataBase(name,ttype,data,x_label,y_label) {
                    pw::ReportBase::setFileExtension("dat");}
    private:
		void reportMetadata(std::ofstream& os) const {streamToDat(os,this->metadata());}
		void reportData(std::ofstream& os) const; 
};

void TrackComplexData::reportData(std::ofstream& os) const 
{
	if(getComplexOp() == pw::ComplexOp::None)
        writeDat1D(os,this->getX(),this->getY(),this->precision());
    else
        writeDat1D(os,this->getX(),this->getOpY(),this->precision());
}




/*
class ReportRealData2D : public pw::ReportBaseRealData2D
{
    public:
        ReportRealData2D(const std::string& name,
            const std::vector<double>& x,
            const std::vector<double>& y,
            const std::vector<double>& z,
            std::string x_label = "x",
            std::string y_label = "y",
            std::string z_label = "z") : pw::ReportBaseRealData2D(
                name,x,y,z,x_label,y_label,z_label) {
                    pw::ReportBase::setFileExtension("dat"); }
        ~ReportRealData2D() {}
    private:
		void reportMetadata(std::ofstream& os) const {streamToDat(os,this->metadata());}
		void reportData(std::ofstream& os) const; 
};

class ReportComplexData2D : public pw::ReportBaseComplexData2D
{
    public:
        ReportComplexData2D(const std::string& name,
            const std::vector<double>& x,
            const std::vector<double>& y,
            const std::vector<dcmplx>& z,
            std::string x_label = "x",
            std::string y_label = "y",
            std::string z_label = "z") : pw::ReportBaseComplexData2D(
                name,x,y,z,x_label,y_label,z_label) {
                    pw::ReportBase::setFileExtension("dat");}
		~ReportComplexData2D() {}
    private:
		void reportMetadata(std::ofstream& os) const {streamToDat(os,this->metadata());}
		void reportData(std::ofstream& os) const; 
};

class ReportTracker : public pw::ReportBaseTracker {
    public:
        ReportTracker(const std::string& nm,const std::string& tlabel="x") : 
            pw::ReportBaseTracker(nm,tlabel) {
                pw::ReportBase::setFileExtension("dat");}

        virtual ~ReportTracker() {}
		void report(std::filesystem::path& filePath,double t){
            std::ofstream os(filePath.string(),std::ios::app);
		    if(ReportBase::metadataOn()){
		        reportMetadata(os);
                ReportBase::setReportMetadata(false);
            }
		    reportTracker(os,t);
        }
    private:
		void reportMetadata(std::ofstream& os) const {streamToDat(os,this->metadata());}
		virtual void reportTracker(std::ofstream& os,double t) = 0;
};

class ReportRealTrackerMax : public ReportTracker
{
    public:
        ReportRealTrackerMax(const std::string& nm,const std::vector<double>& v,
                const std::string& tlabel = "x", const std::string& ylabel="y") :
            ReportTracker(nm,tlabel), m_v(v),m_ylabel(ylabel) {}
        ~ReportRealTrackerMax() {}
    private:
        const std::vector<double>& m_v;
        std::string m_ylabel;
		void reportTracker(std::ofstream& os,double t) const;
};

class ReportComplexTrackerMax : public ReportTracker
{
    public:
        ReportComplexTrackerMax(const std::string& nm,const std::vector<dcmplx>& v,
                const std::string& tlabel = "x", const std::string& ylabel="y") :
            ReportTracker(nm,tlabel), m_v(v),m_ylabel(ylabel) {}
        ~ReportComplexTrackerMax() {}
    private:
        const std::vector<dcmplx>& m_v;
        std::string m_ylabel;
		void reportTracker(std::ofstream& os,double t) const;
};

*/






}

#endif


