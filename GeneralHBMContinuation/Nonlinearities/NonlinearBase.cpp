
#include "NonlinearBase.h"
#include "../Functions.h"

void NonlinearBase::Init(const std::vector<double>& aIntegrationPoints, const std::vector<double>& aBValues, const std::vector<double>& aBProducts, const int& aHarmonicCoeffCount)
{
    cIntegrationPoints = &aIntegrationPoints;
    cBValues = &aBValues;
    cBProducts = &aBProducts;
    mHarmonicCoeffCount = aHarmonicCoeffCount;
    
    // check
    int lCheck2 = aHarmonicCoeffCount * aIntegrationPoints.size();
    if (lCheck2 != aBValues.size()) throw "Number of harmonic coeffs does not pass the check! (first)";
    
    int lCheck = aHarmonicCoeffCount * aHarmonicCoeffCount * aIntegrationPoints.size();
    if (lCheck != aBProducts.size()) throw "Number of harmonic coeffs does not pass the check! (second)";
}

NOX::LAPACK::Matrix<double> NonlinearBase::ComputeJacobianFiniteDifference(const NOX::LAPACK::Vector& aX, const std::function<NOX::LAPACK::Vector(const NOX::LAPACK::Vector&)>& aRHSEval, double aStep) const
{
    if (aStep <= 0) throw "Step must be a positive value!";
    
    int lSize = aX.length();
    NOX::LAPACK::Matrix<double> lSteps(lSize, lSize);
    
    for (int j = 0; j < lSize; j++)
        for (int i = 0; i < lSize; i++) // loop order switched because the matrix is column major
            lSteps(i, j) = aStep;
        
    return ComputeJacobianFiniteDifference(aX, aRHSEval, lSteps);
}
NOX::LAPACK::Matrix<double> NonlinearBase::ComputeJacobianFiniteDifference(const NOX::LAPACK::Vector& aX, const std::function<NOX::LAPACK::Vector(const NOX::LAPACK::Vector&)>& aRHSEval, const NOX::LAPACK::Matrix<double>& aSteps) const
{
    int lSize = aX.length();
    NOX::LAPACK::Matrix<double> lReturnMatrix(lSize, lSize);
    
    for (int j = 0; j < lSize; j++)
    {
        for (int i = 0; i < lSize; i++) // loop order switched because the matrix is column major
        {
            double lStep = aSteps(i, j);
            
            if (lStep <= 0) throw "Step must be a positive value! (" + std::to_string(i) + ", " + std::to_string(j) + ")";
            
            NOX::LAPACK::Vector lVec1(aX);
            NOX::LAPACK::Vector lVec2(aX);
            
            lVec1(j) += lStep;
            lVec2(j) -= lStep;
            
            NOX::LAPACK::Vector lRHS1 = aRHSEval(lVec1);
            NOX::LAPACK::Vector lRHS2 = aRHSEval(lVec2);
            
            double lDerivative = (lRHS1(i) - lRHS2(i)) / 2.0 / lStep;
            
            lReturnMatrix(i, j) = lDerivative;
        }
    }
    
    return lReturnMatrix;
}

NOX::LAPACK::Vector NonlinearBase::FreqToTime(const NOX::LAPACK::Vector& aX, const int& aIntegrationPointIndex) const
{
    if (aIntegrationPointIndex < 0 || aIntegrationPointIndex >= cIntegrationPoints->size())
        throw "Invalid integration point index value!";
    
    const std::vector<double>& lBValues = *cBValues;
    // in time domain
    int lDofCount = aX.length() / mHarmonicCoeffCount;
    
    NOX::LAPACK::Vector lReturnVector(lDofCount);
        
    for (int iDof = 0; iDof < lDofCount; iDof++)
    {
        double lValue = 0.0;
        for (int iHarm = 0; iHarm < mHarmonicCoeffCount; iHarm++)
        {
            int lHarmIndex = GetHBMDofIndex(iDof, iHarm, mHarmonicCoeffCount);
            int lBValIndex = GetBValuesIndex(iHarm, aIntegrationPointIndex, mHarmonicCoeffCount, cIntegrationPoints->size());
            
            lValue += aX(lHarmIndex) * lBValues[lBValIndex];
        }
        lReturnVector(iDof) = lValue;
    }
    
    return lReturnVector;
}
