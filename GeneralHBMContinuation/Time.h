
class Time;

#pragma once
#include <chrono>

class Time
{
private:
    std::chrono::system_clock::time_point mCurrentStart = std::chrono::high_resolution_clock::now();
    
public:
    void Start();
    double Stop(const bool& aPrint = false);
};
