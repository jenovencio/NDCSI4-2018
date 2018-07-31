
#include <functional>

#include "AftSimple.h"
#include "../Functions.h"
#include "../Misc.h"

AftSimple::AftSimple(const int& aIntegrationPointCount, const int& aHarmonicWaveCount, const int& aDofCountTime)
    : AftBase(aIntegrationPointCount, aHarmonicWaveCount, aDofCountTime)
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
    
    mTimeVectors.reserve(cIntegrationPointCount);
    
    for (int i = 0; i < cIntegrationPointCount; i++)
    {
        mTimeVectors.push_back(NOX::LAPACK::Vector(cDofCountTime));
    }
    
    mFreqVector = NOX::LAPACK::Vector(cDofCountTime * cHarmonicCount);
    mFreqMatrix = NOX::LAPACK::Matrix<double>(cDofCountTime * cHarmonicCount, cDofCountTime * cHarmonicCount);
}

const std::vector<NOX::LAPACK::Vector>& AftSimple::FrequencyToTime(const NOX::LAPACK::Vector& aXFreq, const double& aFrequency)
{
    for (int iTime = 0; iTime < cIntegrationPointCount; iTime++)
    {
        for (int iDof = 0; iDof < cDofCountTime; iDof++)
        {
            mTimeVectors[iTime](iDof) = 0.0;
            
            for (int iHarm = 0; iHarm < cHarmonicCount; iHarm++)
            {
                int lHarmIndex = GetHBMDofIndex(iDof, iHarm, cHarmonicCount);
                int lBIndex = GetBValuesIndex(iHarm, iTime, cHarmonicCount, cIntegrationPointCount);
                
                mTimeVectors[iTime](iDof) += aXFreq(lHarmIndex) * mBValues[lBIndex];
            }
        }
    }
    
    return mTimeVectors;
}

const NOX::LAPACK::Vector& AftSimple::TimeToFrequency(const std::vector<NOX::LAPACK::Vector>& aXTime, const double& aFrequency)
{
    // check 
    if (aXTime.size() != cIntegrationPointCount) throw "Number of vectors in the aXTime and number of time integration points don't match!";
    
    double lPeriod = 2 * PI / aFrequency;
    
    // reset the "return" vector
    for (int i = 0; i < mFreqVector.length(); i++) mFreqVector(i) = 0.0;
    
    for (int iTimePoint = 0; iTimePoint < cIntegrationPointCount; iTimePoint++)
    {
        const NOX::LAPACK::Vector& lTimeVector = aXTime[iTimePoint];
        
        for (int iDof = 0; iDof < cDofCountTime; iDof++)
        {
            double lTimeValue = lTimeVector(iDof);
            
            for (int iHarm = 0; iHarm < cHarmonicCount; iHarm++)
            {
                double lHBMIndex = GetHBMDofIndex(iDof, iHarm, cHarmonicCount);
                double lBIndex = GetBValuesIndex(iHarm, iTimePoint, cHarmonicCount, cIntegrationPointCount);
                double lProjValue = lTimeValue * mBValues[lBIndex];
                
                mFreqVector(lHBMIndex) += lProjValue;
            }
        }
    }
    
    mFreqVector.scale(lPeriod / cIntegrationPointCount);
    
    return mFreqVector;
}

const NOX::LAPACK::Matrix<double>& AftSimple::TimeToFrequency(const std::vector<NOX::LAPACK::Matrix<double>>& aXTime, const double& aFrequency)
{
    // check 
    if (aXTime.size() != cIntegrationPointCount) throw "Number of vectors in the aXTime and number of time integration points don't match!";
    
    double lPeriod = 2 * PI / aFrequency;
    
    // reset the "return" matrix
    for (int j = 0; j < mFreqMatrix.numCols(); j++)
        for (int i = 0; i < mFreqMatrix.numRows(); i++) 
            mFreqMatrix(i, j) = 0.0;
        
    for (int iTimePoint = 0; iTimePoint < cIntegrationPointCount; iTimePoint++)
    {
        const NOX::LAPACK::Matrix<double>& lTimeMatrix = aXTime[iTimePoint];
        
        for (int jDof = 0; jDof < cDofCountTime; jDof++)
        {
            for (int iDof = 0; iDof < cDofCountTime; iDof++)
            {
                double lTimeValue = lTimeMatrix(iDof, jDof);
                
                for (int jHarm = 0; jHarm < cHarmonicCount; jHarm++)
                {
                    for (int iHarm = 0; iHarm < cHarmonicCount; iHarm++)
                    {
                        double lHBMIndex1 = GetHBMDofIndex(iDof, iHarm, cHarmonicCount);
                        double lHBMIndex2 = GetHBMDofIndex(jDof, jHarm, cHarmonicCount);
                        double lBProdIndex = GetBProductIndex(iHarm, jHarm, iTimePoint, cHarmonicCount, cIntegrationPointCount);
                        double lProjValue = lTimeValue * mBProducts[lBProdIndex];
                        
                        mFreqMatrix(lHBMIndex1, lHBMIndex2) += lProjValue;
                    }
                }
            }
        }
    }
    
    mFreqMatrix.scale(lPeriod / cIntegrationPointCount);
    
    return mFreqMatrix;
}
