#ifndef PWTHREADS_H
#define PWTHREADS_H

#include <vector>

namespace pw{

class ThreadManager{

    public:
        explicit ThreadManager(unsigned int cnum_threads);  
        std::vector<unsigned int> getBounds(unsigned int size);
        void setNumThreads(unsigned int cnum_threads);
        unsigned int getNumThreads() {return num_threads;}
    private:
        unsigned int num_threads;
};


}

#endif



