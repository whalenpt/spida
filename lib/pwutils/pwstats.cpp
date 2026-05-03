
#include "pwutils/pwstats.h"
#include <map>
#include <string>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <ctime>
#include <cstdio>
#include <stdexcept>

namespace pw{

StatCenter::StatCenter(const std::string& name,unsigned steps_per_log) :
    m_name(name)
{
    setLogFrequency(steps_per_log);
}

void StatCenter::setLogFrequency(unsigned val)
{
    if(val < 1)
        throw std::invalid_argument("Error in StatCenter::setReportFrequency(val): "
                "val must be an integer greater than 0");
    m_steps_per_log = val;
}

void StatCenter::statUpdate(std::ostream& os)
{
  m_stat_requests++;
  if(!(m_stat_requests % m_steps_per_log))
      report(os);
}

void StatCenter::report(std::ostream& os) const
{
    os <<  m_name << std::endl;
    m_timer.report(os);
    m_counter.report(os);
    m_tracker.report(os);
    m_clocker.report(os);
    os << std::endl;
}
 
void Counter::addCounter(const std::string& str,unsigned initial_count)
{
    m_map.try_emplace(str,initial_count);
}

void Counter::increment(const std::string& str,unsigned incr_amount)
{
    auto it = m_map.find(str);
    if(it != m_map.end())
        (*it).second += incr_amount;
    else
        addCounter(str,incr_amount);
}

void Counter::report(std::ostream& os) const
{
    for(auto it = m_map.cbegin(); it != m_map.cend(); it++){
        os << "  " << std::setiosflags(std::ios::left) <<  std::setw(50) << (*it).first + ":";
        os << std::setw(16) << (*it).second << std::endl;
    }
}

void Timer::addTimer(const std::string& str)
{
    using std::chrono::system_clock;
    std::chrono::time_point<system_clock> start;
    std::chrono::duration<double> netTime(0.0);
    m_map.try_emplace(str,netTime);
    m_st.try_emplace(str,start);
}

void Timer::startTimer(const std::string& str)
{
    // Find starting timer element for str
    auto it = m_st.find(str);
    if(it != m_st.end())
        (*it).second = std::chrono::system_clock::now();
    else{
        // No starting timer element for str, so create one
        addTimer(str);
        m_st[str] = std::chrono::system_clock::now();
    }
}

bool Timer::endTimer(const std::string& str)
{
    auto st_it = m_st.find(str);
    if(st_it == m_st.end())
        return false;
    auto map_it = m_map.find(str);
    if(map_it != m_map.end())
        return false;
    std::chrono::duration<double> elapsed_time = std::chrono::system_clock::now() - st_it->second;
    map_it->second += elapsed_time;
    return true;
}

void Timer::report(std::ostream& os) const
{
    for(auto it = m_map.cbegin(); it != m_map.cend(); it++){
        os << "  " << std::setiosflags(std::ios::left) <<  std::setw(50) << (*it).first + ":";
        os << std::setw(16) << std::fixed << (*it).second.count() << std::endl;
    }
}

void Tracker::addTracker(const std::string& str,double val)
{
    m_map.try_emplace(str,val);
}

void Tracker::updateTracker(const std::string& str,double val)
{
    m_map[str] = val;
}

void Tracker::report(std::ostream& os) const
{
    for(auto it = m_map.cbegin(); it != m_map.cend(); it++){
        os << "  " << std::setiosflags(std::ios::left) <<  std::setw(50) << (*it).first + ":";
        os << std::setw(16) << std::scientific << std::setprecision(3) << (*it).second << std::endl;
    }
}

void Clocker::addClock(const std::string& str)
{
    m_map.try_emplace(str,0.0);
    m_st.try_emplace(str,clock());
}

void Clocker::startClock(const std::string& str)
{
    auto it = m_st.find(str);
    if(it != m_st.end())
        (*it).second = std::clock();
    else{
        addClock(str);
        m_st[str] = std::clock();
    }
}

bool Clocker::endClock(const std::string& str)
{
    auto st_it = m_st.find(str);
    if(st_it == m_st.end())
        return false;
    auto map_it = m_map.find(str);
    if(map_it == m_map.end())
        return false;
  
    std::clock_t clockTicks = std::clock() - (*st_it).second;
    auto d_time = static_cast<double>(clockTicks) / static_cast<double>(CLOCKS_PER_SEC);  
    map_it->second += d_time;
    return true;
}

void Clocker::report(std::ostream& os) const
{
    for(auto it = m_map.cbegin(); it != m_map.cend(); it++){
        os << "  " << std::setiosflags(std::ios::left) <<  std::setw(50) << (*it).first + ":";
        os << std::setw(16) << std::fixed << (*it).second << std::endl;
    }
}

}


