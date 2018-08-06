
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
    
//     // storage for frequency to time transform
//     std::vector<NOX::LAPACK::Vector> mTimeVectors;
//     
//     // storage for time to frequency transform (of F values)
//     NOX::LAPACK::Vector mFreqVector;
//     
//     // storage for time to frequency transform (of J values)
//     NOX::LAPACK::Matrix<double> mFreqMatrix;
    
public:
    AftSimple(const int& aIntegrationPointCount, const ProblemParams& aParams);
    
    virtual std::vector<NOX::LAPACK::Vector> FrequencyToTime(const NOX::LAPACK::Vector& aXFreq, const double& aFrequency) override;
    
    virtual NOX::LAPACK::Vector TimeToFrequency(const std::vector<NOX::LAPACK::Vector>& aXTime, const double& aFrequency) override;    
    virtual NOX::LAPACK::Matrix<double> TimeToFrequency(const std::vector<NOX::LAPACK::Matrix<double>>& aXTime, const double& aFrequency) override;
    // get value of a harmonic function for a particular time point
    virtual double GetBValue(const int& aBIndex, const int& aTimePointIndex) override;
    
private:
    // index of a harmonic function and a time point into a serial index
    // aHarmonicCount here means total number of waves, so cosines and sines together
    inline int GetBValuesIndex(const int& aBIndex, const int& aIntPointIndex, const int& aHarmonicCount, const int& aIntPointCount)
    {
        return (aBIndex * aIntPointCount) + aIntPointIndex;
    }
    // index of two harmonic functions and a time point into a serial index
    // aHarmonicCount here means total number of waves, so cosines and sines together
    inline int GetBProductIndex(const int& aBIndex1, const int& aBIndex2, const int& aIntPointIndex, const int& aHarmonicCount, const int& aIntPointCount)
    {
        return (aBIndex1 * aHarmonicCount + aBIndex2) * aIntPointCount + aIntPointIndex;
    }
};
