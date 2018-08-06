
#pragma once
#include <functional>
#include <map>
#include <string>

#include "AftBase.h"
#include "AftSimple.h"
#include "AftFftw.h"
#include "../ProblemParams.h"

// just simple creation functions (should do new because the returned pointer will be deleted later)
const std::map<std::string, std::function<AftBase*(const int&, const ProblemParams&)>> C_AftFactory = 
{
    { std::string("Simple"),    [](const int& aIntegrationPointCount, const ProblemParams& aParams) 
        { 
            return new AftSimple(aIntegrationPointCount, aParams);
        } },
    { std::string("Fftw"),      [](const int& aIntegrationPointCount, const ProblemParams& aParams) 
        { 
            return new AftFftw(aIntegrationPointCount, aParams);
        } },
};
