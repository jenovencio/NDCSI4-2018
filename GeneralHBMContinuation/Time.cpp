
#include <iostream>
#include "Time.h"


void Time::Start()
{
    mCurrentStart = std::chrono::high_resolution_clock::now();
}
double Time::Stop(const bool& aPrint)
{
    std::chrono::system_clock::time_point lNow = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> lDurationMS = lNow - mCurrentStart;
    double lMs = lDurationMS.count();
    
    if (aPrint) std::cout << "Elapsed time: " << lMs << " ms" << std::endl;
    
    return lMs;
}
