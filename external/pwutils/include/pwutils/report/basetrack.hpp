// basetrack.hpp
#pragma once

#include <string>
#include <fstream>
#include <complex>
#include "pwutils/report/basereport.h"
#include "pwutils/pwmath.hpp"

namespace pw{

// Templated Base class for TrackData, specialize to data type, i.e., json or whatever else
// by deriving from this class

template<class T>
class TrackDataBase : public TrackData
{
	public:
        TrackDataBase(const std::string& name,
            TrackType ttype,
            const std::vector<T>& data) :
                TrackData(name,ttype),m_data(data) {}
        virtual ~TrackDataBase() {};
        const std::vector<double>& getX() const {return m_x;}
        const std::vector<T>& getY() const {return m_y;}
        const std::vector<T>& getData() const {return m_data;}
		void updateTracker(double x); 
	private:
        const std::vector<T>& m_data;
        std::vector<double> m_x;
        std::vector<T> m_y;
		virtual void reportData(std::ofstream& os) const = 0;
};

template<class T>
void TrackDataBase<T>::updateTracker(double x)
{
    m_x.push_back(x);
    if(TrackData::getTrackType() == TrackType::Max) {
        m_y.push_back(pw::max(m_data));
    } else if(TrackData::getTrackType() == TrackType::Min){
        m_y.push_back(pw::min(m_data));
    }
}

template<class T>
class TrackComplexDataBase : public TrackDataBase<std::complex<T>>
{
	public:
        TrackComplexDataBase(const std::string& name,
            TrackType ttype,
            const std::vector<std::complex<T>>& data, 
            ComplexOp cmplxop = ComplexOp::None) : 
                TrackDataBase<std::complex<T>>(name,ttype,data),
                m_cmplxop(cmplxop) {}
        virtual ~TrackComplexDataBase() {};
        const std::vector<T>& getOpY() const {return m_opy;}
		void updateTracker(double x); 
		void setComplexOp(ComplexOp cmplxop) {m_cmplxop=cmplxop;}
		ComplexOp getComplexOp() const {return m_cmplxop;}
	private:
	    ComplexOp m_cmplxop;
        std::vector<T> m_opy;
		virtual void reportData(std::ofstream& os) const = 0;
};

template<class T>
void TrackComplexDataBase<T>::updateTracker(double x)
{
    TrackType ttype = TrackData::getTrackType();
    if(m_cmplxop == ComplexOp::None){
        TrackDataBase<std::complex<T>>::updateTracker(x);
        return;
    }
    else if(m_cmplxop == ComplexOp::Power){
        if(ttype == TrackType::Max){
            T val = pw::max(TrackDataBase<std::complex<T>>::getData());
            m_opy.push_back(std::norm(val));
        }
        else if(ttype == TrackType::Min){
            T val = pw::min(TrackDataBase<std::complex<T>>::getData());
            m_opy.push_back(std::norm(val));
        }
    }
}



}


