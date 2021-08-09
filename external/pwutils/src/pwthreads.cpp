
#include "pwutils/pwthreads.h"
#include "pwutils/pwexcept.h"
#include <vector>

namespace pw{

ThreadManager::ThreadManager(unsigned int cnum_threads)
{
    if(cnum_threads < 1){
        throw pw::Exception("ThreadManager constructor failure: "\
           + std::to_string(cnum_threads) + " specified, at least one thread"\
           " is needed to construct ThreadManager.");
    }
    num_threads = cnum_threads;
}


void ThreadManager::setNumThreads(unsigned int cnum_threads)
{
    if(cnum_threads < 1){
        throw pw::Exception("ThreadManager::setNumThreads failure: "\
           + std::to_string(cnum_threads) + " specified, at least one thread"\
           " is needed to for ThreadManager to function properly.");
    }
    num_threads = cnum_threads;
}

std::vector<unsigned int> ThreadManager::getBounds(unsigned int size)
{
    std::vector<unsigned int> bounds(num_threads+1);
    bounds[0] = 0;
    unsigned int chunk_size = size/num_threads;
    for(unsigned int i = 1; i < num_threads; i++)
        bounds[i] = bounds[i-1] + chunk_size;
    bounds[num_threads] = size;
    return bounds;
}

}



