
class AftSimple;

#pragma once
#include "AftBase.h"

class AftSimple : public AftBase
{
private:
    // coordinates of the integration points, relative to a time period T
    // these values would effectively be for period T = 1
    std::vector<double> mIntPointsRelative;
    
    std::vector<double> mBValues;
    std::vector<double> mBProducts;
    
    // storage for frequency to time transform
    std::vector<NOX::LAPACK::Vector> mTimeVectors;
    
    // storage for time to frequency transform (of F values)
    NOX::LAPACK::Vector mFreqVector;
    
    // storage for time to frequency transform (of J values)
    NOX::LAPACK::Matrix<double> mFreqMatrix;
    
public:
    AftSimple(const int& aIntegrationPointCount, const int& aHarmonicWaveCount, const int& aDofCountTime);
    
    virtual const std::vector<NOX::LAPACK::Vector>& FrequencyToTime(const NOX::LAPACK::Vector& aXFreq, const double& aFrequency) override;
    
    virtual const NOX::LAPACK::Vector& TimeToFrequency(const std::vector<NOX::LAPACK::Vector>& aXTime, const double& aFrequency) override;    virtual const NOX::LAPACK::Matrix<double>& TimeToFrequency(const std::vector<NOX::LAPACK::Matrix<double>>& aXTime, const double& aFrequency) override;
};
