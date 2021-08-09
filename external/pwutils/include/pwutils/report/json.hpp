
#ifndef JSON_HPP_
#define JSON_HPP_ 

#include "pwutils/report/basedata.hpp"
#include "pwutils/report/basetrack.hpp"
#include "pwutils/report/reporthelper.h"
#include <complex> 
#include <string> 
#include <vector> 
#include <memory>
#include <fstream>
#include <map> 

namespace json{

using pw::dcmplx;         

void streamToJSON(std::ofstream& os,const pw::metadataMap& str_map);
void writeJSONLabel(std::ofstream& os,const std::string& label,const std::string& indent="\t");
void writeJSONValue(std::ofstream& os,const std::string& value,
		const std::string& indent="",bool end_value=false);

template<typename T>
void writeJSONVector(std::ofstream& os,const std::string& label,const std::vector<T>& v,
		const std::string& indent="\t",bool end_value=false,int precision=pw::REPORT_PRECISION)
{
	writeJSONLabel(os,label,indent);
	os << "[";
	for(unsigned int i = 0; i < v.size()-1; i++){
		os << std::scientific << std::setprecision(precision) << v[i] << ", ";
	}
	std::string end_string = end_value ?  "]" : "],"; 
	os << v.back() << end_string << std::endl;
}

// Complex value specialization of writeJSONVector
template<>
void writeJSONVector(std::ofstream& os,const std::string& label,const std::vector<dcmplx>& v,
		const std::string& indent,bool end_value,int precision)
{
	std::string two_indent = indent + indent;
    writeJSONLabel(os,label,indent);
	os << "{" << std::endl;
    writeJSONLabel(os,"dtype",two_indent);
    writeJSONValue(os,"complex128");
	writeJSONLabel(os,"real",two_indent);
	os << "[";
	for(unsigned int i = 0; i < v.size()-1; i++){
		os << std::scientific << std::setprecision(precision) << v[i].real() << ", ";
	}
	os << v.back().real() << "]," << std::endl;
	writeJSONLabel(os,"imag",two_indent);
	os << "[";
	for(unsigned int i = 0; i < v.size()-1; i++){
		os << std::scientific << std::setprecision(precision) << v[i].imag() << ", ";
	}
	os << v.back().imag() << "]" << std::endl;
	std::string end_string = end_value ?  "}" : "},"; 
	os << indent << end_string << std::endl;
}


void writeJSONPowerVector(std::ofstream& os,const std::string& label,const std::vector<dcmplx>& v,
		const std::string& indent="\t",bool end_value=false,int precision=pw::REPORT_PRECISION);
void writeJSONPhaseVector(std::ofstream& os,const std::string& label,const std::vector<dcmplx>& v,
		const std::string& indent="\t",bool end_value=false,int precision=pw::REPORT_PRECISION);
void mapToJSON(std::ofstream& os,const pw::metadataMap& str_map,bool end_value=false);


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
                    pw::ReportBase::setFileExtension("json");}
        ~ReportData1D() {};
        void report(std::ofstream& os) const {
            os << "{" << std::endl;
            if(pw::ReportBase::metadataOn())
                reportMetadata(os);
            reportData(os);
            os << "}" << std::endl;
        }
    private:
		void reportMetadata(std::ofstream& os) const {streamToJSON(os,this->metadata());}
		void reportData(std::ofstream& os) const; 
};

template<class T1,class T2>
void ReportData1D<T1,T2>::reportData(std::ofstream& os) const 
{
	writeJSONVector(os,this->getLabelX(),this->getX(),"\t",false,this->precision());
    writeJSONVector(os,this->getLabelY(),this->getY(),"\t",true,this->precision());
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
                    pw::ReportBase::setFileExtension("json");}
		~ReportComplexData1D() {}
        void report(std::ofstream& os) const {
            os << "{" << std::endl;
            if(pw::ReportBase::metadataOn())
                reportMetadata(os);
            reportData(os);
            os << "}" << std::endl;
        }
    private:
		void reportMetadata(std::ofstream& os) const {streamToJSON(os,this->metadata());}
		void reportData(std::ofstream& os) const; 
};

template<class T1>
void ReportComplexData1D<T1>::reportData(std::ofstream& os) const
{
    writeJSONVector(os,this->getLabelX(),this->getX(),"\t",false,this->precision());
	if(this->getPhase())
		writeJSONPhaseVector(os,this->getLabelY(),this->getY(),"\t",true,this->precision());
    else if(this->getPower()){
	    writeJSONPowerVector(os,this->getLabelY(),this->getY(),"\t",true,this->precision());
	} else
	   writeJSONVector(os,this->getLabelY(),this->getY(),"\t",true,this->precision());
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
                    pw::ReportBase::setFileExtension("json");}
        ~TrackData() {};
        void report(std::ofstream& os) const {
            os << "{" << std::endl;
            if(pw::ReportBase::metadataOn())
                reportMetadata(os);
            reportData(os);
            os << "}" << std::endl;
        }
    private:
		void reportMetadata(std::ofstream& os) const {streamToJSON(os,this->metadata());}
		void reportData(std::ofstream& os) const; 
};

template<class T>
void TrackData<T>::reportData(std::ofstream& os) const 
{
	writeJSONVector(os,this->getLabelX(),this->getX(),"\t",false,this->precision());
    writeJSONVector(os,this->getLabelY(),this->getY(),"\t",true,this->precision());
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
                pw::TrackComplexDataBase(name,ttype,data,x_label,y_label,cmplxop) {
                    pw::ReportBase::setFileExtension("json");}
        void report(std::ofstream& os) const {
            os << "{" << std::endl;
            if(pw::ReportBase::metadataOn())
                reportMetadata(os);
            reportData(os);
            os << "}" << std::endl;
        }
    private:
		void reportMetadata(std::ofstream& os) const {streamToJSON(os,this->metadata());}
		void reportData(std::ofstream& os) const; 
};

void TrackComplexData::reportData(std::ofstream& os) const 
{
	writeJSONVector(os,this->getLabelX(),this->getX(),"\t",false,this->precision());
	if(getComplexOp() == pw::ComplexOp::None)
        writeJSONVector(os,this->getLabelY(),this->getY(),"\t",true,this->precision());
    else
        writeJSONVector(os,this->getLabelY(),this->getOpY(),"\t",true,this->precision());
}





/*
class ReportRealData2D : public pw::VBReportRealData2D
{
    public:
        ReportRealData2D(const std::string& name,
            const std::vector<double>& x,
            const std::vector<double>& y,
            const std::vector<double>& z,
            std::string x_label = "x",
            std::string y_label = "y",
            std::string z_label = "z") : pw::VBReportRealData2D(
                name,x,y,z,x_label,y_label,z_label) {
                    pw::VBReport::setFileExtension("json"); }
        ~ReportRealData2D() {}
        void report(std::ofstream& os) const {
            os << "{" << std::endl;
            if(VBReportData::metadataOn())
                reportMetadata(os);
            reportData(os);
            os << "}" << std::endl;
        }

    private:
		void reportMetadata(std::ofstream& os) const {mapToJSON(os,this->metadata(),false);}
        void reportData(std::ofstream& os) const;
};

class ReportComplexData2D : public pw::VBReportComplexData2D
{
    public:
        ReportComplexData2D(const std::string& name,
            const std::vector<double>& x,
            const std::vector<double>& y,
            const std::vector<dcmplx>& z,
            std::string x_label = "x",
            std::string y_label = "y",
            std::string z_label = "z") : pw::VBReportComplexData2D(
                name,x,y,z,x_label,y_label,z_label) {
                    pw::VBReport::setFileExtension("json");}
		~ReportComplexData2D() {}
        void report(std::ofstream& os) const {
            os << "{" << std::endl;
            if(VBReportData::metadataOn())
                reportMetadata(os);
            reportData(os);
            os << "}" << std::endl;
        }
    private:
		void reportMetadata(std::ofstream& os) const {mapToJSON(os,this->metadata(),false);}
        void reportData(std::ofstream& os) const;

};

class ReportTracker : public pw::VBReportTracker {
    public:
        ReportTracker(const std::string& nm, const std::string& tlabel="x") : 
            pw::VBReportTracker(nm,tlabel), m_t() {
                pw::VBReport::setFileExtension("json");}
        virtual ~ReportTracker() {}
        std::vector<double>& getT() {return m_t;}
		void report(std::filesystem::path& filePath,double t){
            std::ofstream os(filePath.string(),std::ios::trunc);
            os << "{" << std::endl;
		    if(VBReport::metadataOn()){
                reportMetadata(os);
            }
		    reportTracker(os,t);
            os << "}" << std::endl;
        }
    private:
        std::vector<double> m_t;
		void reportMetadata(std::ofstream& os) const {mapToJSON(os,this->metadata(),false);}
		virtual void reportTracker(std::ofstream& os,double t) = 0;
};

class ReportRealTrackerMax : public ReportTracker
{
    public:
        ReportRealTrackerMax(const std::string& nm,const std::vector<double>& v,
                const std::string& tlabel = "x",const std::string& ylabel="y") :
            ReportTracker(nm,tlabel), m_v(v),m_ylabel(ylabel) {}
        ~ReportRealTrackerMax() {}
        std::string getLabelY() const {return m_ylabel;}
    private:
        const std::vector<double>& m_v;
        std::string m_ylabel;
        std::vector<double> m_y;
		void reportTracker(std::ofstream& os,double t);
};

class ReportComplexTrackerMax : public ReportTracker
{
    public:
        ReportComplexTrackerMax(const std::string& nm,const std::vector<dcmplx>& v,
                const std::string& tlabel = "x", const std::string& ylabel="y") :
            ReportTracker(nm,tlabel), m_v(v),m_ylabel(ylabel) {}
        ~ReportComplexTrackerMax() {}
        std::string getLabelY() const {return m_ylabel;}
    private:
        const std::vector<dcmplx>& m_v;
        std::string m_ylabel;
        std::vector<double> m_y;
		void reportTracker(std::ofstream& os,double t);
};

*/





}

#endif


