
#include "pwutils/pwthreads.h"
#include "pwutils/pwexcept.h"
#include <vector>

namespace pw{

ThreadManager::ThreadManager(unsigned num_threads)
{
    if(num_threads < 1){
        throw pw::Exception("ThreadManager constructor failure: "\
           + std::to_string(num_threads) + " specified, at least one thread"\
           " is needed to construct ThreadManager.");
    }
    m_num_threads = num_threads;
}


void ThreadManager::setNumThreads(unsigned num_threads)
{
    if(num_threads < 1){
        throw pw::Exception("ThreadManager::setNumThreads failure: "\
           + std::to_string(num_threads) + " specified, at least one thread"\
           " is needed to for ThreadManager to function properly.");
    }
    m_num_threads = num_threads;
}

std::vector<unsigned> ThreadManager::getBounds(unsigned size) const
{
    std::vector<unsigned> bounds(m_num_threads+1);
    bounds[0] = 0;
    unsigned chunk_size = size/m_num_threads;
    for(unsigned i = 1; i < m_num_threads; i++)
        bounds[i] = bounds[i-1] + chunk_size;
    bounds[m_num_threads] = size;
    return bounds;
}

}



