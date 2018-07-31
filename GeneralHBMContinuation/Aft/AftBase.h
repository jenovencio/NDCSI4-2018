
class AftBase;

#pragma once
#include <vector>

#include "NOX_LAPACK_Matrix.H"
#include "NOX_LAPACK_Vector.H"

class AftBase
{
private:
protected:
    const int cIntegrationPointCount;
    const int cHarmonicCount;
    const int cHarmonicWaveCount;
    const int cDofCountTime;
    
public:
    AftBase(const int& aIntegrationPointCount, const int& aHarmonicWaveCount, const int& aDofCountTime)
        : cIntegrationPointCount(aIntegrationPointCount), cHarmonicCount(2 * aHarmonicWaveCount - 1), cHarmonicWaveCount(aHarmonicWaveCount),
        cDofCountTime(aDofCountTime) { }
    
    virtual const std::vector<NOX::LAPACK::Vector>& FrequencyToTime(const NOX::LAPACK::Vector& aXFreq, const double& aFrequency) = 0;
    
    virtual const NOX::LAPACK::Vector& TimeToFrequency(const std::vector<NOX::LAPACK::Vector>& aXTime, const double& aFrequency) = 0;
    virtual const NOX::LAPACK::Matrix<double>& TimeToFrequency(const std::vector<NOX::LAPACK::Matrix<double>>& aXTime, const double& aFrequency) = 0;
};
