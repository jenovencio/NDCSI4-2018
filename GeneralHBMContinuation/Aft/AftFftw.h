
class AftFftw;

#pragma once
#include "AftBase.h"

class AftFftw : public AftBase
{
public:
    AftFftw(const int& aIntegrationPointCount, const ProblemParams& aParams);
    
    virtual std::vector<NOX::LAPACK::Vector> FrequencyToTime(const NOX::LAPACK::Vector& aXFreq, const double& aFrequency) override;
        
    virtual NOX::LAPACK::Vector TimeToFrequency(const std::vector<NOX::LAPACK::Vector>& aXTime, const double& aFrequency) override;
    virtual NOX::LAPACK::Matrix<double> TimeToFrequency(const std::vector<NOX::LAPACK::Matrix<double>>& aXTime, const double& aFrequency) override;
    // get value of a harmonic function for a particular time point
    virtual double GetBValue(const int& aBIndex, const int& aTimePointIndex) override;
};
