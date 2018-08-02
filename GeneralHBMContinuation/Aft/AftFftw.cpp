
#include "AftFftw.h"


AftFftw::AftFftw(const int& aIntegrationPointCount, const ProblemParams& aParams)
    : AftBase(aIntegrationPointCount, aParams)
{
    
}

const std::vector<NOX::LAPACK::Vector>& AftFftw::FrequencyToTime(const NOX::LAPACK::Vector& aXFreq, const double& aFrequency)
{
    throw "AftFftw is not implemented yet!";
}

const NOX::LAPACK::Vector& AftFftw::TimeToFrequency(const std::vector<NOX::LAPACK::Vector>& aXTime, const double& aFrequency)
{
    throw "AftFftw is not implemented yet!";
}
const NOX::LAPACK::Matrix<double>& AftFftw::TimeToFrequency(const std::vector<NOX::LAPACK::Matrix<double>>& aXTime, const double& aFrequency)
{
    throw "AftFftw is not implemented yet!";
}
