
#include <functional>

#include "AftSimple.h"
#include "Functions.h"
#include "Misc.h"

AftSimple::AftSimple(const int& aIntegrationPointCount, const int& aHarmonicWaveCount)
    : cIntegrationPointCount(aIntegrationPointCount), cHarmonicCount(2 * aHarmonicWaveCount - 1), cHarmonicWaveCount(aHarmonicWaveCount)
{
    if (cIntegrationPointCount <= 0) throw "Number of integration points must be positive!";
    if (aHarmonicWaveCount <= 0) throw "Number of harmonic waves must be positive!";
    
    mIntPointsRelative = GetRelativeTimePoints(cIntegrationPointCount);
    
    double* lTempArray1 = new double[cIntegrationPointCount * cHarmonicCount];
    
    std::function<double(double)> lB;
    
    for (int i = 0; i < cHarmonicCount; i++)
    {
        int lWaveNumber = (i + 1) / 2;
        
        if (i == 0)             lB = [](double aX) { return 1.0; };
        else if (i % 2 == 1)    lB = [lWaveNumber](double aX) { return std::cos(2 * PI * lWaveNumber * aX); };
        else                    lB = [lWaveNumber](double aX) { return std::sin(2 * PI * lWaveNumber * aX); };
        
        for (int iIntPoint = 0; iIntPoint < cIntegrationPointCount; iIntPoint++)
        {
            int lIndex = GetBValuesIndex(i, iIntPoint, cHarmonicCount, cIntegrationPointCount);
            double lIntPointPos = mIntPointsRelative[iIntPoint];
            
            double lValue = lB(lIntPointPos);
            
            lTempArray1[lIndex] = lValue;
        }
    }
    
    mBValues.assign(lTempArray1, lTempArray1 + cIntegrationPointCount * cHarmonicCount);
    
    double* lTempArray2 = new double[cIntegrationPointCount * cHarmonicCount * cHarmonicCount];
    
    std::function<double(double)> lB1;
    std::function<double(double)> lB2;
    
    for (int i = 0; i < cHarmonicCount; i++)
    {
        int lWaveNumber = (i + 1) / 2;
        
        if (i == 0)             lB1 = [](double aX) { return 1.0; };
        else if (i % 2 == 1)    lB1 = [lWaveNumber](double aX) { return std::cos(2 * PI * lWaveNumber * aX); };
        else                    lB1 = [lWaveNumber](double aX) { return std::sin(2 * PI * lWaveNumber * aX); };
        
        for (int j = 0; j < cHarmonicCount; j++)
        {
            int lWaveNumber2 = (j + 1) / 2;
            
            if (j == 0)             lB2 = [](double aX) { return 1.0; };
            else if (j % 2 == 1)    lB2 = [lWaveNumber2](double aX) { return std::cos(2 * PI * lWaveNumber2 * aX); };
            else                    lB2 = [lWaveNumber2](double aX) { return std::sin(2 * PI * lWaveNumber2 * aX); };
            
            for (int iIntPoint = 0; iIntPoint < cIntegrationPointCount; iIntPoint++)
            {
                int lIndex = GetBProductIndex(i, j, iIntPoint, cHarmonicCount, cIntegrationPointCount);
                double lIntPointPos = mIntPointsRelative[iIntPoint];
                
                double lProduct = lB1(lIntPointPos) * lB2(lIntPointPos);
                lTempArray2[lIndex] = lProduct;
            }
        }
    }
    mBProducts.assign(lTempArray2, lTempArray2 + cIntegrationPointCount * cHarmonicCount * cHarmonicCount);
}

const std::vector<NOX::LAPACK::Vector>& FrequencyToTime(const NOX::LAPACK::Vector& aXFreq, const double& aFrequency)
{
    
}

const NOX::LAPACK::Vector& TimeToFrequency(const std::vector<NOX::LAPACK::Vector>& aXTime, const double& aFrequency)
{
    
}
