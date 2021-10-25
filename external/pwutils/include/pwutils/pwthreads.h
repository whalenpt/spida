//pwthreads.h
#pragma once
#include <vector>

namespace pw{

class ThreadManager{

    public:
        explicit ThreadManager(unsigned num_threads);  
        std::vector<unsigned> getBounds(unsigned size) const;
        void setNumThreads(unsigned num_threads);
        unsigned getNumThreads() const {return m_num_threads;}
    private:
        unsigned m_num_threads;
};


}




