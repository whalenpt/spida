// basedata.hpp
#pragma once

#include <string>
#include <fstream>
#include <cassert>
#include "pwutils/report/basereport.h"

namespace pw{

template<class T1,class T2>
class ReportDataBase1D : public ReportData1D
{
	public:
        ReportDataBase1D(const std::string& name,
            const std::vector<T1>& x, 
            const std::vector<T2>& y) :
                ReportData1D(name), m_x(x),m_y(y) {}
        virtual ~ReportDataBase1D() {};
        const std::vector<T1>& getX() const {return m_x;}
        const std::vector<T2>& getY() const {return m_y;}
	private:
        const std::vector<T1>& m_x;
        const std::vector<T2>& m_y;
		virtual void reportData(std::ofstream& os) const = 0;
};

template<class T1,class T2>
class ReportComplexDataBase1D : public ReportDataBase1D<T1,std::complex<T2>>
{
	public:
        ReportComplexDataBase1D(const std::string& name,
            const std::vector<T1>& x,const std::vector<std::complex<T2>>& y) : 
                ReportDataBase1D<T1,std::complex<T2>>(name,x,y),
                    m_power(false),m_phase(false) {}
        virtual ~ReportComplexDataBase1D() {};
		void setPower(bool val) {m_power= val;}  
		void setPhase(bool val) {m_phase = val;}  
		bool getPower() const {return m_power;}
		bool getPhase() const {return m_phase;}
	private:
  		bool m_power;
  		bool m_phase;
		virtual void reportData(std::ofstream& os) const = 0;
};


template<class T1,class T2,class T3>
class ReportDataBase2D : public ReportData2D
{
    public:
        ReportDataBase2D(const std::string& name,
        const std::vector<T1>& x, 
        const std::vector<T2>& y, 
        const std::vector<T3>& z) :
            ReportData2D(name), m_x(x),m_y(y),m_z(z) {}
        virtual ~ReportDataBase2D() {};
        const std::vector<T1>& getX() const {return m_x;}
        const std::vector<T2>& getY() const {return m_y;}
        const std::vector<T3>& getZ() const {return m_z;}
    private:
        const std::vector<T1>& m_x;
        const std::vector<T2>& m_y;
        const std::vector<T3>& m_z;
        virtual void reportData(std::ofstream& os) const = 0;
}; 

template<class T1,class T2,class T3>
class ReportComplexDataBase2D : public ReportDataBase2D<T1,T2,std::complex<T3>>
{
    public:
        ReportComplexDataBase2D(const std::string& name,
            const std::vector<T1>& x,
            const std::vector<T2>& y,
	        const std::vector<std::complex<T3>>& z) :
                ReportDataBase2D<T1,T2,std::complex<T3>>(name,x,y,z),
                    m_power(false),m_phase(false) {}
        virtual ~ReportComplexDataBase2D() {};
        void setPower(bool val) {m_power= val;}  
        void setPhase(bool val) {m_phase = val;}  
        bool getPower() const {return m_power;}
        bool getPhase() const {return m_phase;}
    private:
        bool m_power;
        bool m_phase;
        virtual void reportData(std::ofstream& os) const = 0;
};


}




