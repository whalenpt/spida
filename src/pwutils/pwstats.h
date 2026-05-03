// pwstats.h
#pragma once

#include <chrono>
#include <ctime>
#include <iostream>
#include <map>
#include <string>
#include <string_view>

namespace pw{

using unsignedPair = std::pair<std::string,unsigned>;
using unsignedMap = std::map<std::string,unsigned,std::less<>>;
using dblPair = std::pair<std::string,double>;
using dblMap = std::map<std::string,double, std::less<>>;

using chronoPair = std::pair<std::string, std::chrono::time_point<std::chrono::system_clock> >;
using chronoMap = std::map<std::string,std::chrono::time_point<std::chrono::system_clock>, std::less<>>;
using chronoDurPair = std::pair<std::string, std::chrono::duration<double> >;
using chronoDurMap = std::map<std::string,std::chrono::duration<double>, std::less<>>;
using clockPair = std::pair<std::string,clock_t>;
using clockMap = std::map<std::string,clock_t, std::less<>>;

class Counter{
    public:
        Counter() = default;
        void addCounter(const std::string& str,unsigned initial_count=0);
        void increment(const std::string& str,unsigned incr_amount=1);
        void report(std::ostream& os = std::cout) const;
    private:
        unsignedMap m_map; 
};

class Timer{
    public:
        Timer() = default;
        void addTimer(const std::string& str);
        void startTimer(const std::string& str);
        bool endTimer(const std::string& str);
        void report(std::ostream& os = std::cout) const;
    private:
        chronoDurMap m_map; 
        chronoMap m_st; 
};

class Clocker{
    public:
        Clocker() = default;
        void addClock(const std::string& str);
        void startClock(const std::string& str);
        bool endClock(const std::string& str);
        void report(std::ostream& os = std::cout) const;
    private:
        dblMap m_map; 
        clockMap m_st; 
};

class Tracker{
    public:
        Tracker() = default;
        void addTracker(const std::string& str,double val = 0.0);
        void updateTracker(const std::string& str,double val);
        void report(std::ostream& os = std::cout) const;
    private:
        dblMap m_map; 
};

class StatCenter{
    public:
        StatCenter(const std::string& name = "STATS",unsigned steps_per_log = 1);
        ~StatCenter() = default;
        void statUpdate(std::ostream& os = std::cout);
        void report(std::ostream& os = std::cout) const;
        void addCounter(const std::string& str,unsigned start_count = 0) {
            m_counter.addCounter(str,start_count);}
        void incrementCounter(const std::string& str,unsigned incr_amount = 1) {
            m_counter.increment(str,incr_amount);
        }
        void addTimer(std::string const& str) {m_timer.addTimer(str);}
        void startTimer(const std::string& str) {m_timer.startTimer(str);}
        bool endTimer(const std::string& str) {return m_timer.endTimer(str);}

        void addClock(const std::string& str) {m_clocker.addClock(str);}
        void startClock(const std::string& str) {m_clocker.startClock(str);}
        bool endClock(const std::string& str) {return m_clocker.endClock(str);}

        void addTracker(const std::string& str,double val = 0.0) {m_tracker.addTracker(str,val);}
        void updateTracker(const std::string& str,double val) {m_tracker.updateTracker(str,val);}

        void setLogFrequency(unsigned val); 
        void setHeader(std::string_view str) {m_name = str;}
    private:
        Timer m_timer;
        Counter m_counter;
        Tracker m_tracker;
        Clocker m_clocker;
        unsigned m_stat_requests{0};
        unsigned m_steps_per_log;
        std::string m_name;
};

}