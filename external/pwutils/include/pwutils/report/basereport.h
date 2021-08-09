#ifndef BASEREPORT_H
#define BASEREPORT_H

#include <string>
#include <fstream>
#include <filesystem>
#include <memory>
#include "pwutils/report/reporthelper.h"

namespace pw{

enum class TrackType { Max,
    Min
};

enum class ComplexOp {None,
    Power
};

class ReportBase{
	public:
		ReportBase(const std::string& nm) 
		  : m_name(nm),m_report_metadata(true),
		  m_metadata_map(),
		  m_extension(),
		  m_precision(REPORT_PRECISION) {}
		virtual ~ReportBase() {};
		virtual void report(std::ofstream& os) const {
		    if(metadataOn())
		        reportMetadata(os);
		    reportData(os);
        }
		friend std::ofstream& operator<<(std::ofstream& os,const ReportBase& def){ def.report(os);
			return os; }
		friend std::ofstream& operator<<(std::ofstream& os,const ReportBase* def){ def->report(os);
			return os; }
		friend std::ofstream& operator<<(std::ofstream& os,const std::unique_ptr<ReportBase> def){ \
		    def->report(os); return os; }

		void setItem(const std::string& key,double val);
		void setItem(const std::string&,const std::string&); 
		void removeItem(const std::string&); 
		void setReportMetadata(bool val) {m_report_metadata = val;}
		void setFileExtension(const std::string& extension) {
		    m_extension=extension;}
		void setPrecision(int precision) {m_precision = precision;}
		int precision() const {return m_precision;}
		std::string name() const {return m_name;}
		bool metadataOn() const {return m_report_metadata;}
		const metadataMap& metadata() const {return m_metadata_map;} 
		std::string fileExtension() const {return m_extension;}
        std::filesystem::path filePath(const std::filesystem::path& dir_path) const
            { return pw::filePath(dir_path,m_name,m_extension);}
        std::filesystem::path filePath(const std::filesystem::path& dir_path,int repNum) const
            { return pw::filePath(dir_path,m_name,repNum,m_extension);}
	private:
		const std::string m_name;
		bool m_report_metadata;
        metadataMap m_metadata_map; 
        std::string m_extension;
		int m_precision;
		virtual void reportMetadata(std::ofstream& os) const = 0;
		virtual void reportData(std::ofstream& os) const = 0;
};

// Need a non-templated base class for holding all ReportData1D instances in an STL container
// without specifying a data type (which is dictated by the template subclass)
class ReportData1D : public ReportBase
{
	public:
		ReportData1D(const std::string& nm) 
		  : ReportBase(nm) {}
		virtual ~ReportData1D() {}
	private:
		virtual void reportData(std::ofstream& os) const = 0;
};


// Need a non-templated base class for holding all TrackData instances in an STL container
// without specifying a data type (which is dictated by the template subclass)
class TrackData : public ReportBase
{
	public:
		TrackData(const std::string& nm,TrackType ttype) 
		  : ReportBase(nm), m_ttype(ttype) {}
		virtual ~TrackData() {}
		virtual void updateTracker(double t) = 0; // assume time t is a double value (or convert to)
		const TrackType getTrackType() {return m_ttype;}
		void setTrackType(TrackType ttype) {m_ttype = ttype;}
    private:
        TrackType m_ttype;
		virtual void reportData(std::ofstream& os) const = 0;
};


/*
class ReportBaseData2D : public ReportBaseData
{
	public:
        ReportBaseData2D(const std::string& name) : ReportBaseData(name) {}
        virtual ~ReportBaseData2D() {};
		std::string getLabelX() const {return m_xlabel;}
		std::string getLabelY() const {return m_ylabel;}
		std::string getLabelZ() const {return m_zlabel;}
		void setLabelX(const std::string& xlabel) {m_xlabel=xlabel;}
		void setLabelY(const std::string& ylabel) {m_ylabel=ylabel;}
		void setLabelZ(const std::string& zlabel) {m_zlabel=zlabel;}
	private:
		std::string m_xlabel;
		std::string m_ylabel;
		std::string m_zlabel;
		virtual void reportData(std::ofstream& os) const = 0;
};


class ReportBaseRealData2D : public ReportBaseData2D
{
    public:
        ReportBaseRealData2D(const std::string& name,
            const std::vector<double>& x,
            const std::vector<double>& y,
            const std::vector<double>& z,
            std::string x_label = "x",
            std::string y_label = "y",
            std::string z_label = "z") : ReportBaseData2D(name), m_x(x), m_y(y),
                m_z(z) {setLabelX(x_label); setLabelY(y_label); setLabelZ(z_label);}
        virtual ~ReportBaseRealData2D() {}
		const std::vector<double>& getX() const {return m_x;}
		const std::vector<double>& getY() const {return m_y;}
		const std::vector<double>& getZ() const {return m_z;}
	private:
        const std::vector<double>& m_x;
        const std::vector<double>& m_y;
        const std::vector<double>& m_z;
		virtual void reportData(std::ofstream& os) const = 0;
};

class ReportBaseComplexData2D : public ReportBaseData2D
{
    public:
        ReportBaseComplexData2D(const std::string& name,
            const std::vector<double>& x,
            const std::vector<double>& y,
            const std::vector<dcmplx>& z,
            std::string x_label = "x",
            std::string y_label = "y",
            std::string z_label = "z") : ReportBaseData2D(name), m_x(x), m_y(y),
                m_z(z),m_power(false),m_phase(false) {
                    setLabelX(x_label); setLabelY(y_label); setLabelZ(z_label);}
        virtual ~ReportBaseComplexData2D() {}
		void setPower(bool val) {m_power= val;}  
		void setPhase(bool val) {m_phase = val;}  
		bool getPower() const {return m_power;}
		bool getPhase() const {return m_phase;}
		const std::vector<double>& getX() const {return m_x;}
		const std::vector<double>& getY() const {return m_y;}
		const std::vector<dcmplx>& getZ() const {return m_z;}
	private:
        const std::vector<double>& m_x;
        const std::vector<double>& m_y;
        const std::vector<dcmplx>& m_z;
  		bool m_power;
  		bool m_phase;
		virtual void reportData(std::ofstream& os) const = 0;
};

class ReportBaseTracker : public ReportBase
{
	public:
		ReportBaseTracker(const std::string& nm,std::string tlabel="x") : 
		    ReportBase(nm), m_tlabel(tlabel) {}
		virtual ~ReportBaseTracker() {}
		std::string getLabelT() const {return m_tlabel;}
		void setLabelT(const std::string& tlabel) {m_tlabel = tlabel;}
		virtual void report(std::filesystem::path& filePath,double t) = 0; 
	private:
		virtual void reportMetadata(std::ofstream& os) const = 0;
		virtual void reportTracker(std::ofstream& os,double t) = 0;
        std::string m_tlabel;
};

*/



}




#endif
