
#include "AftFftw.h"


AftFftw::AftFftw(const int& aIntegrationPointCount, const int& aHarmonicWaveCount, const int& aDofCountTime)
    : AftBase(aIntegrationPointCount, aHarmonicWaveCount, aDofCountTime)
{
    
}

const std::vector<NOX::LAPACK::Vector>& AftFftw::FrequencyToTime(const NOX::LAPACK::Vector& aXFreq, const double& aFrequency)
{
}

const NOX::LAPACK::Vector& AftFftw::TimeToFrequency(const std::vector<NOX::LAPACK::Vector>& aXTime, const double& aFrequency)
{
}
const NOX::LAPACK::Matrix<double>& AftFftw::TimeToFrequency(const std::vector<NOX::LAPACK::Matrix<double>>& aXTime, const double& aFrequency)
{
}
