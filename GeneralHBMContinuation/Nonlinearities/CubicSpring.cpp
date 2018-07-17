
#include "CubicSpring.h"
#include "../Functions.h"
#include "../Misc.h"

NOX::LAPACK::Vector CubicSpring::ComputeRHS(const NOX::LAPACK::Vector& aX, const double& aFrequency) const
{
    NOX::LAPACK::Vector lReturnVector(aX.length());
    const std::vector<double>& lIntPoints = *cIntegrationPoints;
    const std::vector<double>& lBValues = *cBValues;
    const std::vector<double>& lBProducts = *cBProducts;
    
    // check 
    if (aFrequency <= 0) throw "Frequency must be a positive value!";
    if (mHarmonicCoeffCount <= 0) throw "Number of harmonic coefficients must be a positive integer!";
    if (aX.length() % mHarmonicCoeffCount != 0) throw "Size of the problem in frequency domain (" + std::to_string(aX.length()) + ") is not divisible by number of harmonic coefficients (" + std::to_string(mHarmonicCoeffCount) + ")!";
    
    // in time domain
    int lDofCount = aX.length() / mHarmonicCoeffCount;
    // period
    double lT = 2 * PI / aFrequency;
    double lTimeStep = lT / cIntegrationPoints->size();
    
    for (int iIntPoint = 0; iIntPoint < cIntegrationPoints->size(); iIntPoint++)
    {
        NOX::LAPACK::Vector lXTime = FreqToTime(aX, iIntPoint);
        
        // calculate the nonlinearity ini the time domain
        NOX::LAPACK::Vector lNonlin(lXTime.length());
        
        for (int iSpring = 0; iSpring < mSprings.size(); iSpring++)
        {
            int lDof = mSprings[iSpring].DofIndex;
            double lStiffCoeff = mSprings[iSpring].StiffnessCoeff;
            
            lNonlin(lDof) += lStiffCoeff * lXTime(lDof) * lXTime(lDof) * lXTime(lDof);
        }
        
        for (int iHarm = 0; iHarm < mHarmonicCoeffCount; iHarm++)
        {
            double lBValIndex = GetBValuesIndex(iHarm, iIntPoint, mHarmonicCoeffCount, cIntegrationPoints->size());
            double lBValue = lBValues[lBValIndex];
            
            NOX::LAPACK::Vector lNonlinProjected = lNonlin;
            lNonlinProjected.scale(lBValue);
            
            for (int iDof = 0; iDof < lDofCount; iDof++)
            {
                double lHarmInd = GetHBMDofIndex(iDof, iHarm, mHarmonicCoeffCount);
                lReturnVector(lHarmInd) += lNonlinProjected(iDof);
            }
        }
    }
    
    lReturnVector.scale(lTimeStep);
    
    return lReturnVector;
}

NOX::LAPACK::Matrix<double> CubicSpring::ComputeJacobian(const NOX::LAPACK::Vector& aX, const double& aFrequency) const
{
//     NOX::LAPACK::Matrix<double> lReturnMatrix;
//     const std::vector<double>& lIntPoints = *cIntegrationPoints;
//     const std::vector<double>& lBValues = *cBValues;
//     const std::vector<double>& lBProducts = *cBProducts;
//     
//     return lReturnMatrix;
    
    std::function<NOX::LAPACK::Vector(const NOX::LAPACK::Vector&)> lRHSEval = [aFrequency, this](const NOX::LAPACK::Vector& aIn) { return ComputeRHS(aIn, aFrequency); };
    
    return ComputeJacobianFiniteDifference(aX, lRHSEval, 1e-8);
}

void CubicSpring::AddCubicSpring(const CubicSpringDef& aDef)
{
    if (aDef.DofIndex < 0) throw "Dof index can not be negative!";
    if (aDef.StiffnessCoeff < 0) throw "Cubic stiffness can not be negative!";
    if (aDef.StiffnessCoeff == 0) return;
    
    mSprings.push_back(aDef);
}
void CubicSpring::AddCubicSpring(const int& aDofIndex, const double& aStiffnessCoeff)
{
    CubicSpringDef lNewSpring;
    lNewSpring.DofIndex = aDofIndex;
    lNewSpring.StiffnessCoeff = aStiffnessCoeff;
    
    AddCubicSpring(lNewSpring);
}
void CubicSpring::ClearSprings()
{
    mSprings.clear();
}
