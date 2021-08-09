#ifndef BASEDATA_HPP
#define BASEDATA_HPP

#include <string>
#include <fstream>
#include "pwutils/report/basereport.h"

namespace pw{

template<class T1,class T2>
class ReportDataBase1D : public ReportData1D
{
	public:
        ReportDataBase1D(const std::string& name,
            const std::vector<T1>& x, 
            const std::vector<T2>& y, 
            const std::string& x_label = "x",
            const std::string& y_label = "y") : 
                ReportData1D(name), m_x(x),m_y(y),
                m_xlabel(x_label),m_ylabel(y_label) {
                    ReportBase::setItem("xlabel",x_label);
                    ReportBase::setItem("ylabel",y_label);
                }
        virtual ~ReportDataBase1D() {};
        const std::vector<T1>& getX() const {return m_x;}
        const std::vector<T2>& getY() const {return m_y;}
		std::string getLabelX() const {return m_xlabel;}
		std::string getLabelY() const {return m_ylabel;}
		void setLabelX(const std::string& xlabel) {m_xlabel=xlabel;
            ReportBase::setItem("xlabel",xlabel);
		}
		void setLabelY(const std::string& ylabel) {m_ylabel=ylabel;
            ReportBase::setItem("ylabel",ylabel);
		}
	private:
        const std::vector<T1>& m_x;
        const std::vector<T2>& m_y;
		std::string m_xlabel;
		std::string m_ylabel;
		virtual void reportData(std::ofstream& os) const = 0;
};

template<class T1>
class ReportComplexDataBase1D : public ReportDataBase1D<T1,dcmplx>
{
	public:
        ReportComplexDataBase1D(const std::string& name,
            const std::vector<T1>& x,const std::vector<dcmplx>& y, 
            const std::string& x_label = "x",
            const std::string& y_label = "y") : 
                ReportDataBase1D<T1,dcmplx>(name,x,y,x_label,y_label),
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



}




#endif
