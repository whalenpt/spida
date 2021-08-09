
#include "pwutils/pwstats.h"
#include <map>
#include <string>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <ctime>
#include <cstdio>

namespace pw{

StatCenter::StatCenter(std::string nm,int stPerRp)
{
  statRequests = 0;
  statReports = 0;
  name = nm;
  stepsPerReport = stPerRp;
}

void StatCenter::statUpdate(std::ostream& os)
{
  statRequests++;
  if(stepsPerReport == 0)
    return;
  else{
    if(!(statRequests % stepsPerReport))
      report(os);
  }
}

void StatCenter::report(std::ostream& os) const
{
  os <<  name << std::endl;
  timer.report(os);
  counter.report(os);
  tracker.report(os);
  clocker.report(os);
  os << std::endl;
}
 
void Counter::addCounter(std::string str,int val)
{
  map.insert(intPair(str,0));
  sz.insert(intPair(str,val));
}
bool Counter::increment(std::string str)
{
  countMap::iterator it;
  it = map.find(str);
  if(it != map.end()){
    int val = (*it).second;
    countMap::iterator itb;
    itb = sz.find(str);
    if(itb != sz.end()){
      int incr = (*itb).second;
      val+=incr; 
      map[str] = val;
      return true;
    }
    else
      return false;
  }
  else
    return false;
}

void Counter::incrementAll()
{
  countMap::iterator it = map.begin();
  while(it != map.end()){
    std::string str = (*it).first;
    int val = (*it).second;
    int incr = sz[str];
    val+=incr; 
    map[str] = val;
    it++;
  }
}

void Counter::report(std::ostream& os) const
{
  countMap::const_iterator it = map.begin();
  while(it != map.end()){
    os << "  " << std::setiosflags(std::ios::left) <<  std::setw(50) << (*it).first + ":";
    os << std::setw(16) << (*it).second << std::endl;
    it++;
  }
}

/*
void Timer::addTimer(std::string str)
{
  map.insert(dblPair(str,0.0));
  st.insert(timePair(str,time(0)));
  end.insert(timePair(str,time(0)));
}

bool Timer::startTimer(std::string str)
{
  timeMap::iterator it;
  it = st.find(str);
  if(it != st.end()){
    time(&(*it).second);
    return true;
  }
  else
    return false;
}

bool Timer::endTimer(std::string str)
{
  sizet sttime;
  sizet end_time;
  timeMap::iterator it;
  it = end.find(str);
  if(it != end.end()){
    time(&(*it).second);
    end_time = (*it).second;
    timeMap::iterator its;
    its = st.find(str);
    if(its != st.end())
      st_time = (*its).second;
    else
      return false;
    double dTime = difftime(end_time,st_time);
    double netTime = map[str];
    netTime += dTime;
    map[str] = netTime;
    return true;
  }
  else
    return false;
}

void Timer::report(std::ostream& os) const
{
  dblMap::const_iterator it = map.begin();
  while(it != map.end()){
    os << "  " << std::setiosflags(std::ios::left) <<  std::setw(50) << (*it).first + ":";
    os << std::setw(16) << std::fixed << (*it).second << std::endl;
    it++;
  }
}


*/

void Timer::addTimer(std::string str)
{
  using std::chrono::system_clock;
  std::chrono::time_point<system_clock> start,endt;
  std::chrono::duration<double> netTime(0.0);
  map.insert(chronoDurPair(str,netTime));
  st.insert(chronoPair(str,start));
  end.insert(chronoPair(str,endt));
}

bool Timer::startTimer(std::string str)
{
  chronoMap::iterator it;
  it = st.find(str);
  if(it != st.end()){
    (*it).second = std::chrono::system_clock::now();
    return true;
  }
  else
    return false;
}

bool Timer::endTimer(std::string str)
{
  std::chrono::duration<double> elapsed_time;
  chronoMap::iterator end_it = end.find(str);
  if(end_it != end.end()){
    (*end_it).second = std::chrono::system_clock::now();
    chronoMap::iterator st_it = st.find(str);
    if(st_it != st.end())
      elapsed_time = (*end_it).second - (*st_it).second;
    else
      return false;
    map[str] = map[str] + elapsed_time;
    return true;
  }
  else
    return false;
}

void Timer::report(std::ostream& os) const
{
  chronoDurMap::const_iterator it = map.begin();
  while(it != map.end()){
    os << "  " << std::setiosflags(std::ios::left) <<  std::setw(50) << (*it).first + ":";
    os << std::setw(16) << std::fixed << (*it).second.count() << std::endl;
    it++;
  }
}

void Tracker::addTracker(std::string str,double val)
{
  map.insert(dblPair(str,val));
}

void Tracker::updateTracker(std::string str,double val)
{
  map[str] = val;
}

void Tracker::report(std::ostream& os) const
{
  dblMap::const_iterator it = map.begin();
  while(it != map.end()){
    os << "  " << std::setiosflags(std::ios::left) <<  std::setw(50) << (*it).first + ":";
    os << std::setw(16) << std::scientific << std::setprecision(3) << (*it).second << std::endl;
    it++;
  }
}

void Clocker::addClock(std::string str)
{
  map.insert(dblPair(str,0.0));
  st.insert(clockPair(str,clock()));
}

bool Clocker::startClock(std::string str)
{
  clockMap::iterator it;
  it = st.find(str);
  if(it != st.end()){
    (*it).second = std::clock();
    return true;
  }
  else
    return false;
}

bool Clocker::endClock(std::string str)
{
  clockMap::iterator it;
  it = st.find(str);
  if(it != st.end()){
    std::clock_t clockTicks = std::clock() - (*it).second;
    double dTime = clockTicks / (double) CLOCKS_PER_SEC;  
    double netTime = map[str];
    netTime += dTime;
    map[str] = netTime;
    return true;
  }
  else
    return false;
}

void Clocker::report(std::ostream& os) const
{
  dblMap::const_iterator it = map.begin();
  while(it != map.end()){
    os << "  " << std::setiosflags(std::ios::left) <<  std::setw(50) << (*it).first + ":";
    os << std::setw(16) << std::fixed << (*it).second << std::endl;
    it++;
  }
}

}


