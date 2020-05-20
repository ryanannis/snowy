#pragma once

#include "Common.hpp"
#include <vector>
#include <thread>
#include <cmath> 
#include <functional>

template<typename T, typename Func>
void IterateThread(std::vector<T>& data, Uint low, Uint high, Func f) {
    for(Uint i = low; i < high; i++) 
    {
        f(data[i]);
    }
}

template<typename T, typename Func>
void IterateThreadConst(const std::vector<T>& data, Uint low, Uint high, Func f) {
    for(Uint i = low; i < high; i++) 
    {
        f(data[i]);
    }
}


class MTIterator
{
public:
    MTIterator(Uint numthreads);
    
    // Todo:  Yeah, this should be able to use things other than vectors... but we're not for now
    template<typename T, typename Func>
    void IterateOverVector(std::vector<T>& data, Func f) {
        // Split into NumThreads buckets without creating empty buckets 
        // (unless thre are less elements than there are buckets)
        std::vector<std::thread> threads;

        Uint idx = 0;
        Uint partsLeft = data.size();
        for (int threadsLeft = mNumThreads; threadsLeft > 0 ; threadsLeft--)
        {
            Uint parts = std::ceil(Float(partsLeft) / Float(threadsLeft));
            partsLeft -= parts;
            threads.push_back(std::thread(IterateThread<T, Func>, std::ref(data), idx, idx + parts, f));
            idx += parts;
        }

        // todo:  ya I'm creating empty threads and joining them for the edge cases
        for(std::thread& t : threads)
        {
            t.join();
        }
    }

    // Todo:  Remove this
    template<typename T, typename Func>
    void IterateOverVector(const std::vector<T>& data, Func f) {
        // Split into NumThreads buckets without creating empty buckets 
        // (unless thre are less elements than there are buckets)
        std::vector<std::thread> threads;

        Uint idx = 0;
        Uint partsLeft = data.size();
        for (int threadsLeft = mNumThreads; threadsLeft > 0 ; threadsLeft--)
        {
            Uint parts = std::ceil(Float(partsLeft) / Float(threadsLeft));
            partsLeft -= parts;
            threads.push_back(std::thread(IterateThreadConst<T, Func>, std::ref(data), idx, idx + parts, f));
            idx += parts;
        }

        // todo:  ya I'm creating empty threads and joining them for the edge cases
        for(std::thread& t : threads)
        {
            t.join();
        }
    }

private:
    const Uint mNumThreads;
};
