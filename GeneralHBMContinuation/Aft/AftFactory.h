
#pragma once
#include <functional>
#include <map>
#include <string>

#include "AftBase.h"
#include "AftSimple.h"
#include "AftFftw.h"

// just simple creation functions (should do new because the returned pointer will be deleted later)
const std::map<std::string, std::function<AftBase*(const int&, const int&, const int&)>> C_AftFactory = 
{
    { std::string("Simple"),    [](const int& aIntegrationPointCount, const int& aHarmonicWaveCount, const int& aDofCountTime) 
        { 
            return new AftSimple(aIntegrationPointCount, aHarmonicWaveCount, aDofCountTime);
        } },
    { std::string("Fftw"),      [](const int& aIntegrationPointCount, const int& aHarmonicWaveCount, const int& aDofCountTime) 
        { 
            return new AftFftw(aIntegrationPointCount, aHarmonicWaveCount, aDofCountTime);
        } },
};
