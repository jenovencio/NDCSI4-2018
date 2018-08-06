
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
    
    virtual std::vector<NOX::LAPACK::Vector> FrequencyToTime(const NOX::LAPACK::Vector& aXFreq, const double& aFrequency) = 0;
    
    // dft of results in time domain
    virtual NOX::LAPACK::Vector TimeToFrequency(const std::vector<NOX::LAPACK::Vector>& aXTime, const double& aFrequency) = 0;
    // dft of results in time domain
    // the matrices are expected to be derivatives of forces in time domain by harmonic coefficients
    // so their size should be: (dof) x (dof * harm)
    virtual NOX::LAPACK::Matrix<double> TimeToFrequency(const std::vector<NOX::LAPACK::Matrix<double>>& aXTime, const double& aFrequency) = 0;
    // get value of a harmonic function for a particular time point
    virtual double GetBValue(const int& aBIndex, const int& aTimePointIndex) = 0;
};
