// basereport.h
#pragma once

#include <string>
#include <fstream>
#include <filesystem>
#include <memory>
#include <cassert>
#include "pwutils/pwdefs.h"
#include "pwutils/report/reporthelper.h"

namespace pw{

enum class TrackType { Max, Min };
enum class ComplexOp {None, Power };

class ReportBase{
	public:
		ReportBase(const std::string& nm) 
		  : m_name(nm),m_report_metadata(true),
		  m_metadata_map(),
		  m_extension(),
	      m_dirpath("outfolder") {}
		virtual ~ReportBase() {};

		void report(std::ofstream& os) const;
		void report(std::ofstream& os,unsigned rep_num) const;
		friend std::ofstream& operator<<(std::ofstream& os,const ReportBase& def);
		friend std::ofstream& operator<<(std::ofstream& os,const ReportBase* def);
		friend std::ofstream& operator<<(std::ofstream& os,const std::unique_ptr<ReportBase> def);

        void setFileSignature(FileSignature file_sig);
        void setDataSignature(DataSignature data_sig);
        void setOperatorSignature(OperatorSignature op_sig);
        FileSignature fileSignature() const;
        DataSignature dataSignature() const;
        OperatorSignature operatorSignature() const;

		void setName(const std::string& nm) {m_name = nm;}
		void setItem(const std::string& key,double val);
		void setItem(const std::string& key,float val);
		void setItem(const std::string& key,int val);
		void setItem(const std::string& key,const std::string&); 
		void removeItem(const std::string&); 
		void setReportMetadata(bool val) {m_report_metadata = val;}
		void setFileExtension(const std::string& extension) {
		    m_extension=extension;}
		void setDirPath(const std::filesystem::path& dirpath) {
            pw::createDirectory(dirpath,false);
		    m_dirpath = dirpath;}

		std::string getName() const {return m_name;}
		const metadataMap& getMetadata() const {return m_metadata_map;} 
		std::string getFileExtension() const {return m_extension;}

        std::filesystem::path dirpath() const {return m_dirpath;}
        std::filesystem::path path() const {
            return pw::filePath(m_dirpath,m_name,m_extension);}
        std::filesystem::path path(int rep_num) const {
            return pw::filePath(m_dirpath,m_name,rep_num,m_extension);}
		bool metadataOn() const {return m_report_metadata;}

	private:
		std::string m_name;
		bool m_report_metadata;
        metadataMap m_metadata_map; 
        std::string m_extension;
        virtual void reportImplement(std::ofstream& os) const = 0;
		virtual void reportMetadata(std::ofstream& os) const = 0;
		virtual void reportData(std::ofstream& os) const = 0;
        std::filesystem::path m_dirpath;
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

// Need a non-templated base class for holding all ReportData1D instances in an STL container
// without specifying a data type (which is dictated by the template subclass)
class ReportData2D : public ReportBase
{
	public:
		ReportData2D(const std::string& nm) :
            ReportBase(nm),
            m_strideX(1), 
            m_strideY(1) {}
		virtual ~ReportData2D() {}
		unsigned int getStrideX() const {return m_strideX;}
		unsigned int getStrideY() const {return m_strideY;}
		void setStrideX(unsigned int strideX) { assert (strideX >= 1);
		    m_strideX  = strideX;}
		void setStrideY(unsigned int strideY) { assert (strideY >= 1);
		    m_strideY  = strideY;}
	private:
		virtual void reportData(std::ofstream& os) const = 0;
		unsigned int m_strideX;
		unsigned int m_strideY;
};



}




