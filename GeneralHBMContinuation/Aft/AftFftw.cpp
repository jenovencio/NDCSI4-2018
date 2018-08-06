
#include "AftFftw.h"


AftFftw::AftFftw(const int& aIntegrationPointCount, const ProblemParams& aParams)
    : AftBase(aIntegrationPointCount, aParams)
{
    
}

std::vector<NOX::LAPACK::Vector> AftFftw::FrequencyToTime(const NOX::LAPACK::Vector& aXFreq, const double& aFrequency)
{
    throw "AftFftw is not implemented yet!";
}

NOX::LAPACK::Vector AftFftw::TimeToFrequency(const std::vector<NOX::LAPACK::Vector>& aXTime, const double& aFrequency)
{
    throw "AftFftw is not implemented yet!";
}
NOX::LAPACK::Matrix<double> AftFftw::TimeToFrequency(const std::vector<NOX::LAPACK::Matrix<double>>& aXTime, const double& aFrequency)
{
    throw "AftFftw is not implemented yet!";
}
double AftFftw::GetBValue(const int& aBIndex, const int& aTimePointIndex)
{
    if (aBIndex < 0 || aBIndex >= cProblemParams.HarmonicCount) throw "Invalid B index!";
    if (aTimePointIndex < 0 || aTimePointIndex >= cIntegrationPointCount) throw "Invalid time point index!";
    
    throw "AftFftw is not implemented yet!";
}
