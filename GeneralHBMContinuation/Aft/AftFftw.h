
class AftFftw;

#pragma once
#include "AftBase.h"

class AftFftw : public AftBase
{
public:
    AftFftw(const int& aIntegrationPointCount, const ProblemParams& aParams);
    
    virtual const std::vector<NOX::LAPACK::Vector>& FrequencyToTime(const NOX::LAPACK::Vector& aXFreq, const double& aFrequency) override;
        
    virtual const NOX::LAPACK::Vector& TimeToFrequency(const std::vector<NOX::LAPACK::Vector>& aXTime, const double& aFrequency) override;
    virtual const NOX::LAPACK::Matrix<double>& TimeToFrequency(const std::vector<NOX::LAPACK::Matrix<double>>& aXTime, const double& aFrequency) override;
};
