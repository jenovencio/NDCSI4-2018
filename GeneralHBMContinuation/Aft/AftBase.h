
class AftBase;

#pragma once
#include <vector>

#include "NOX_LAPACK_Matrix.H"
#include "NOX_LAPACK_Vector.H"

#include "../ProblemParams.h"

class AftBase
{
private:
protected:
    const int cIntegrationPointCount;
    const ProblemParams cProblemParams;
    
public:
    AftBase(const int& aIntegrationPointCount, const ProblemParams& aProblemParams)
        : cIntegrationPointCount(aIntegrationPointCount), cProblemParams(aProblemParams)
        {
            if (cIntegrationPointCount <= 0) throw "Number of integration points must be positive!";
            if (aProblemParams.HarmonicCount <= 0) throw "Number of harmonic coefficients must be positive!";
        }
    
    virtual const std::vector<NOX::LAPACK::Vector>& FrequencyToTime(const NOX::LAPACK::Vector& aXFreq, const double& aFrequency) = 0;
    
    virtual const NOX::LAPACK::Vector& TimeToFrequency(const std::vector<NOX::LAPACK::Vector>& aXTime, const double& aFrequency) = 0;
    virtual const NOX::LAPACK::Matrix<double>& TimeToFrequency(const std::vector<NOX::LAPACK::Matrix<double>>& aXTime, const double& aFrequency) = 0;
};
