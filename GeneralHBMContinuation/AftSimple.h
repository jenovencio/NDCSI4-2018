
class AftSimple;

#pragma once


#include "AftBase.h"

class AftSimple : public AftBase
{
private:
    const int cIntegrationPointCount;
    const int cHarmonicCount;
    const int cHarmonicWaveCount;
    // coordinates of the integration points, relative to a time period T
    // these values would effectively be for period T = 1
    std::vector<double> mIntPointsRelative;
    
    std::vector<double> mBValues;
    std::vector<double> mBProducts;
    
public:
    AftSimple(const int& aIntegrationPointCount, const int& aHarmonicWaveCount);
    virtual const std::vector<NOX::LAPACK::Vector>& FrequencyToTime(const NOX::LAPACK::Vector& aXFreq, const double& aFrequency) override;
    virtual const NOX::LAPACK::Vector& TimeToFrequency(const std::vector<NOX::LAPACK::Vector>& aXTime, const double& aFrequency) override;
};
