
#include "pwutils/report/basetrack.hpp"

namespace pw{

void TrackComplexDataBase::updateTracker(double x)
{
    TrackType ttype = getTrackType();
    if(m_cmplxop == ComplexOp::None){
        TrackDataBase<dcmplx>::updateTracker(x);
        return;
    }
    else if(m_cmplxop == ComplexOp::Power){
        if(ttype == TrackType::Max){
            dcmplx val = pw::max(getData());
            m_opy.push_back(std::norm(val));
        }
        else if(ttype == TrackType::Min){
            dcmplx val = pw::min(getData());
            m_opy.push_back(std::norm(val));
        }
    }
}

}








