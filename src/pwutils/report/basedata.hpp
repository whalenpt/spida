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
        ~ReportDataBase1D() override = default;
        const std::vector<T1>& getX() const {return m_x;}
        const std::vector<T2>& getY() const {return m_y;}
	private:
        const std::vector<T1>& m_x;
        const std::vector<T2>& m_y;
};

template<class T1,class T2>
class ReportComplexDataBase1D : public ReportDataBase1D<T1,std::complex<T2>>
{
	public:
        ReportComplexDataBase1D(const std::string& name,
            const std::vector<T1>& x,const std::vector<std::complex<T2>>& y) : 
                ReportDataBase1D<T1,std::complex<T2>>(name,x,y) {}
        ~ReportComplexDataBase1D() override = default;
		void setPower(bool val) {m_power= val;}  
		void setPhase(bool val) {m_phase = val;}  
		bool getPower() const {return m_power;}
		bool getPhase() const {return m_phase;}
	private:
  		bool m_power{false};
  		bool m_phase{false};
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
        ~ReportDataBase2D() override = default;
        const std::vector<T1>& getX() const {return m_x;}
        const std::vector<T2>& getY() const {return m_y;}
        const std::vector<T3>& getZ() const {return m_z;}
    private:
        const std::vector<T1>& m_x;
        const std::vector<T2>& m_y;
        const std::vector<T3>& m_z;
}; 

template<class T1,class T2,class T3>
class ReportComplexDataBase2D : public ReportDataBase2D<T1,T2,std::complex<T3>>
{
    public:
        ReportComplexDataBase2D(const std::string& name,
            const std::vector<T1>& x,
            const std::vector<T2>& y,
	        const std::vector<std::complex<T3>>& z) :
                ReportDataBase2D<T1,T2,std::complex<T3>>(name,x,y,z) {}
        ~ReportComplexDataBase2D() override = default;
        void setPower(bool val) {m_power= val;}  
        void setPhase(bool val) {m_phase = val;}  
        bool getPower() const {return m_power;}
        bool getPhase() const {return m_phase;}
    private:
        bool m_power{false};
        bool m_phase{false};
};

}
