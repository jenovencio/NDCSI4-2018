
#include "fftw3.h"

#include "NonlinearBase.h"
#include "../Functions.h"
#include "../Misc.h"
#include "../Time.h"

void NonlinearBaseTD::Finalise()
{
    NonlinearBaseH::Finalise();
    
    // this is the final time this gets executed
//     mNonzeroFPositions = NonzeroFPositions();
    
//     // check dof validity
//     for (int iDof = 0; iDof < mNonzeroFPositions.size(); iDof++)
//     {
//         int lDof = mNonzeroFPositions[iDof];
//         if (lDof < 0) throw "Dof index can not be negative!";
//         // we can do this check here because at this point the object is already initialised
//         if (lDof >= mProblemParams->DofCountPhysical) throw "Dof index (" + std::to_string(lDof) + ") exceeds the size of the problem! (" + std::to_string(mProblemParams->DofCountPhysical) + ")";
//     }
}

NOX::LAPACK::Matrix<double> NonlinearBaseTD::ComputeDFDH(const NOX::LAPACK::Vector& aX, const int& aTimePointIndex)
{
    // check
//     if (aFrequency <= 0) throw "Frequency must be a positive value!";
    if (mProblemParams.HarmonicCount <= 0) throw "Number of harmonic coefficients must be a positive integer!";
    
    int lDofCount = mProblemParams.DofCountPhysical;
    
    NOX::LAPACK::Matrix<double> lReturnMatrix(lDofCount, mProblemParams.DofCountHBM);
    
    const NOX::LAPACK::Matrix<double>& lJTimeDomain = ComputeJacobianTimeDomain(aX);
    
    for (int jDof = 0; jDof < lDofCount; jDof++)
    {
        for (int iDof = 0; iDof < lDofCount; iDof++)
        {
            double lTimeValue = lJTimeDomain(iDof, jDof);
            
            for (int jHarm = 0; jHarm < mProblemParams.HarmonicCount; jHarm++)
            {
                double lBValue = cAft->GetBValue(jHarm, aTimePointIndex);
                double lHarmValue = lTimeValue * lBValue;
                
                int lColIndexHBM = GetHBMDofIndex(jDof, jHarm, mProblemParams.HarmonicCount);
                lReturnMatrix(iDof, lColIndexHBM) = lHarmValue;
            }
        }
    }
    
    return lReturnMatrix;
}





















// NOX::LAPACK::Matrix<double> NonlinearBaseTD::ComputeJacobianFiniteDifference(const NOX::LAPACK::Vector& aX, const std::function<NOX::LAPACK::Vector(const NOX::LAPACK::Vector&)>& aFEval, double aStep) const
// {
//     if (aStep <= 0) throw "Step must be a positive value!";
//     
//     int lSize = aX.length();
//     NOX::LAPACK::Matrix<double> lSteps(lSize, lSize);
//     
//     for (int j = 0; j < lSize; j++)
//         for (int i = 0; i < lSize; i++) // loop order switched because the matrix is column major
//             lSteps(i, j) = aStep;
//         
//     return ComputeJacobianFiniteDifference(aX, aFEval, lSteps);
// }
// 
// NOX::LAPACK::Matrix<double> NonlinearBaseTD::ComputeJacobianFiniteDifference(const NOX::LAPACK::Vector& aX, const std::function<NOX::LAPACK::Vector(const NOX::LAPACK::Vector&)>& aFEval, const NOX::LAPACK::Matrix<double>& aSteps) const
// {
//     int lSize = aX.length();
//     NOX::LAPACK::Matrix<double> lReturnMatrix(lSize, lSize);
//     
//     for (int j = 0; j < lSize; j++)
//     {
//         for (int i = 0; i < lSize; i++) // loop order switched because the matrix is column major
//         {
//             double lStep = aSteps(i, j);
//             
//             if (lStep <= 0) throw "Step must be a positive value! (" + std::to_string(i) + ", " + std::to_string(j) + ")";
//             
//             NOX::LAPACK::Vector lVec1(aX);
//             NOX::LAPACK::Vector lVec2(aX);
//             
//             lVec1(j) += lStep;
//             lVec2(j) -= lStep;
//             
//             NOX::LAPACK::Vector lF1 = aFEval(lVec1);
//             NOX::LAPACK::Vector lF2 = aFEval(lVec2);
//             
//             double lDerivative = (lF1(i) - lF2(i)) / 2.0 / lStep;
//             
//             lReturnMatrix(i, j) = lDerivative;
//         }
//     }
//     
//     return lReturnMatrix;
// }
// 
// NOX::LAPACK::Matrix<double> NonlinearBaseTD::ComputeJacobianFiniteDifference(const NOX::LAPACK::Vector& aX, const double& aFrequency, double aStep) const
// {
//     std::function<NOX::LAPACK::Vector(const NOX::LAPACK::Vector&)> lFEval = [aFrequency, this](const NOX::LAPACK::Vector& aIn) { return ComputeF(aIn, aFrequency); };
//     
//     return ComputeJacobianFiniteDifference(aX, lFEval, aStep);
// }
// 
// NOX::LAPACK::Matrix<double> NonlinearBaseTD::ComputeJacobianFiniteDifference(const NOX::LAPACK::Vector& aX, const double& aFrequency, const NOX::LAPACK::Matrix<double>& aSteps) const
// {
//     std::function<NOX::LAPACK::Vector(const NOX::LAPACK::Vector&)> lFEval = [aFrequency, this](const NOX::LAPACK::Vector& aIn) { return ComputeF(aIn, aFrequency); };
//     
//     return ComputeJacobianFiniteDifference(aX, lFEval, aSteps);
// }
// // Same functions as above, but using the ComputeF as the aFEval function (in time domain)
// NOX::LAPACK::Matrix<double> NonlinearBaseTD::ComputeJacobianFiniteDifferenceTD(const NOX::LAPACK::Vector& aX, const NOX::LAPACK::Vector& aXPrev, double aStep) const
// {
//     std::function<NOX::LAPACK::Vector(const NOX::LAPACK::Vector&)> lFEval = [aXPrev, this](const NOX::LAPACK::Vector& aIn) 
//     {
//         FResult lRes = ComputeFTimeDomain(aIn, aXPrev);
//         return lRes.FValues;
//     };
//     
//     return ComputeJacobianFiniteDifference(aX, lFEval, aStep);
// }
// NOX::LAPACK::Matrix<double> NonlinearBaseTD::ComputeJacobianFiniteDifferenceTD(const NOX::LAPACK::Vector& aX, const NOX::LAPACK::Vector& aXPrev, const NOX::LAPACK::Matrix<double>& aSteps) const
// {
//     std::function<NOX::LAPACK::Vector(const NOX::LAPACK::Vector&)> lFEval = [aXPrev, this](const NOX::LAPACK::Vector& aIn) 
//     {
//         FResult lRes = ComputeFTimeDomain(aIn, aXPrev);
//         return lRes.FValues;
//     };
//     
//     return ComputeJacobianFiniteDifference(aX, lFEval, aSteps);
// }
