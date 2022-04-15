#pragma once
#include <thread>
#ifdef GLMNETPP_HAS_PTHREAD
#include <pthread.h>
#include <iostream>
#endif

namespace glmnetpp {

template <class PIDType>
void set_affinity(int i, PIDType id)
{
    static_cast<void>(i);
    static_cast<void>(id);
    auto n_cpus = std::thread::hardware_concurrency();
    if (n_cpus <= 1) return;

#if defined(GLMNETPP_HAS_PTHREAD) && defined(__linux__)
    cpu_set_t mask;
    int status;
    CPU_ZERO(&mask);
    CPU_SET(i % n_cpus, &mask);
    status = pthread_setaffinity_np(id, sizeof(mask), &mask);
    if (status != 0)
    {
        std::cerr << "Error calling pthread_setaffinity_np: " << status << "\n";
    }
#endif
}

} // namespace glmnetpp
