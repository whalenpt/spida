
#ifndef STATCENTER_H_
#define STATCENTER_H_ 

#include <ctime>
#include <map>
#include <string>
#include <iostream>
#include <chrono>

namespace pw{

typedef std::pair<std::string,int> intPair;
typedef std::map<std::string,int> countMap;

typedef std::pair<std::string,time_t> timePair;
typedef std::map<std::string,time_t> timeMap;

typedef std::pair<std::string, std::chrono::time_point<std::chrono::system_clock> > chronoPair;
typedef std::map<std::string,std::chrono::time_point<std::chrono::system_clock> > chronoMap;
typedef std::pair<std::string, std::chrono::duration<double> > chronoDurPair;
typedef std::map<std::string,std::chrono::duration<double> > chronoDurMap;

typedef std::pair<std::string,double> dblPair;
typedef std::map<std::string,double> dblMap;
typedef std::pair<std::string,clock_t> clockPair;
typedef std::map<std::string,clock_t> clockMap;

class Counter{
  public:
    Counter() {}
    void addCounter(std::string str,int val);
    bool increment(std::string str);
    void incrementAll();
    void report(std::ostream& os = std::cout) const;
  private:
    countMap map; 
    countMap sz; 
};

class Timer{
  public:
    Timer() {} 
    void addTimer(std::string str);
    bool startTimer(std::string str);
    bool endTimer(std::string str);
    void report(std::ostream& os = std::cout) const;
  private:
    chronoDurMap map; 
    chronoMap st; 
    chronoMap end; 
};

class Clocker{
  public:
    Clocker() {} 
    void addClock(std::string str);
    bool startClock(std::string str);
    bool endClock(std::string str);
    void report(std::ostream& os = std::cout) const;
  private:
    dblMap map; 
    clockMap st; 
};

class Tracker{
  public:
    Tracker() {}
    void addTracker(std::string str,double val = 0.0);
    void updateTracker(std::string str,double val);
    void report(std::ostream& os = std::cout) const;
  private:
    dblMap map; 
};

class StatCenter{
  public:
    StatCenter(std::string nm = "STATS",int stPerRp = 0);
    ~StatCenter() {}
    void statUpdate(std::ostream& os = std::cout);
    void report(std::ostream& os = std::cout) const;
    void addCounter(std::string str,int val) {counter.addCounter(str,val);}
    void incrementCounter(std::string str) {counter.increment(str);}
    void incrementAllCounters() {counter.incrementAll();}
    void addTimer(std::string str) {timer.addTimer(str);}
    void startTimer(std::string str) {timer.startTimer(str);}
    void endTimer(std::string str) {timer.endTimer(str);}
    void addClock(std::string str) {clocker.addClock(str);}
    void startClock(std::string str) {clocker.startClock(str);}
    void endClock(std::string str) {clocker.endClock(str);}
    void addTracker(std::string str,double val = 0.0) {tracker.addTracker(str,val);}
    void updateTracker(std::string str,double val) {tracker.updateTracker(str,val);}
    void setReportFrequency(int val) {stepsPerReport = val;}
    void setHeader(std::string str) {name = str;}
  private:
    Timer timer;
    Counter counter;
    Tracker tracker;
    Clocker clocker;
    int statRequests;
    int statReports;
    int stepsPerReport;
    std::string name;
};



}

#endif


