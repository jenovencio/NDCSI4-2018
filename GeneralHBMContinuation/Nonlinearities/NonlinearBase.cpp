
#include "fftw3.h"

#include "NonlinearBase.h"
#include "../Functions.h"
#include "../Misc.h"
#include "../Time.h"

NonlinearBase::~NonlinearBase()
{
    
}
void NonlinearBase::LoadFromFile(const std::string& aFilePath)
{
    std::string lTypeName = ClassName();
    
    throw "Can not load a nonlinearity from file \"" + aFilePath + "\"! This class does not support loading from files.";
}

void NonlinearBase::Init(AftBase* const aAft, const int& aHarmonicCoeffCount, const int& aDofCountTimeDomain)
{
    if (mIsInitialised) throw "This nonlinearity has already been initialised, can not initialise again!";
    if (mIsFinalised) throw "This nonlinearity has already been finalised! Fix the way you treat the object. It first needs to be initialised, then finalised!";
    
    cAft = aAft;
    mHarmonicCoeffCount = aHarmonicCoeffCount;
        
    mIsInitialised = true;
    
    mDofCountTimeDomain = aDofCountTimeDomain;
}

// frequency domain to frequency domain
const NOX::LAPACK::Vector& NonlinearBase::ComputeF(const NOX::LAPACK::Vector& aX, const double& aFrequency) const
{
    CheckStatus();
    
    AftBase& lAft = *cAft;
    
// check 
//     if (aFrequency <= 0) throw "Frequency must be a positive value!";
    if (mHarmonicCoeffCount <= 0) throw "Number of harmonic coefficients must be a positive integer!";
    if (aX.length() % mHarmonicCoeffCount != 0) throw "Size of the problem in frequency domain (" + std::to_string(aX.length()) + ") is not divisible by number of harmonic coefficients (" + std::to_string(mHarmonicCoeffCount) + ")!";
    
    // in time domain
    int lDofCount = aX.length() / mHarmonicCoeffCount;
    
    const std::vector<NOX::LAPACK::Vector>& lXTimeAll = lAft.FrequencyToTime(aX, aFrequency);
    
    std::vector<NOX::LAPACK::Vector> lFTimeAll;
    lFTimeAll.reserve(lXTimeAll.size());
    
    NOX::LAPACK::Vector lXTimePrev(lDofCount);
    
    int lLoopCount = NumberOfPrepLoops() + 1; // number of loop over the period
    
    NOX::LAPACK::Vector lXTimeAvg(lDofCount);
        
    for (int iDof = 0; iDof < lDofCount; iDof++)
    {
        int lHarmIndex = GetHBMDofIndex(iDof, 0, mHarmonicCoeffCount);
        lXTimeAvg(iDof) = aX(lHarmIndex);
    }
//     std::cout << "XTimeAvg: " << lXTimeAvg << std::endl;
    
    for (int iLoop = 0; iLoop < lLoopCount; iLoop++)
    {
        for (int iTimePoint = 0; iTimePoint < lXTimeAll.size(); iTimePoint++)
        {
            const NOX::LAPACK::Vector& lXTime = lXTimeAll[iTimePoint];
            if (iTimePoint == 0 && iLoop == 0) lXTimePrev = lXTimeAvg;
            
            // calculate the nonlinearity in time domain
            FResult lNonlinResult = ComputeFTimeDomain(lXTime, lXTimePrev);
            
            if (!IsCorrectingX() && lNonlinResult.XCorrSet)
                throw "The code logic of this class is wrong! The class says it's not correcting the X values, but at the same time sets corrections in it's F evaluations.";
            
            if (IsCorrectingX() && lNonlinResult.XCorrSet)
                lXTimePrev = lNonlinResult.XCorr;
            else
                lXTimePrev = lXTime;
            
            if (iLoop == lLoopCount - 1)
            {
                lFTimeAll.push_back(lNonlinResult.FValues);
            }
        }
    }
    
    const NOX::LAPACK::Vector& lReturnVector = lAft.TimeToFrequency(lFTimeAll, aFrequency);
    
//     double lMs = lTime.Stop();
//     std::cout << "Full F eval time: " << lMs << " ms" << std::endl;
    
    return lReturnVector;
}
// frequency domain to frequency domain
const NOX::LAPACK::Matrix<double>& NonlinearBase::ComputeJacobian(const NOX::LAPACK::Vector& aX, const double& aFrequency) const
{
    CheckStatus();
    
    AftBase& lAft = *cAft;
    
    // check
//     if (aFrequency <= 0) throw "Frequency must be a positive value!";
    if (mHarmonicCoeffCount <= 0) throw "Number of harmonic coefficients must be a positive integer!";
    if (aX.length() % mHarmonicCoeffCount != 0) throw "Size of the problem in frequency domain (" + std::to_string(aX.length()) + ") is not divisible by number of harmonic coefficients (" + std::to_string(mHarmonicCoeffCount) + ")!";
    
    // in time domain
    int lDofCount = aX.length() / mHarmonicCoeffCount;
    
    const std::vector<NOX::LAPACK::Vector>& lXTimeAll = lAft.FrequencyToTime(aX, aFrequency);
    
    std::vector<NOX::LAPACK::Matrix<double>> lJTimeAll;
    lJTimeAll.reserve(lXTimeAll.size());
    
    NOX::LAPACK::Vector lXTimePrev;
    
    int lLoopCount = NumberOfPrepLoops() + 1; // number of loop over the period
    
    NOX::LAPACK::Vector lXTimeAvg(lDofCount);
    
    for (int iDof = 0; iDof < lDofCount; iDof++)
    {
        int lHarmIndex = GetHBMDofIndex(iDof, 0, mHarmonicCoeffCount);
        lXTimeAvg(iDof) = aX(lHarmIndex);
    }
//     std::cout << "XTimeAvg: " << lXTimeAvg << std::endl;
    
    for (int iLoop = 0; iLoop < lLoopCount; iLoop++)
    {
        for (int iIntPoint = 0; iIntPoint < lXTimeAll.size(); iIntPoint++)
        {
//             lFftInvTime.Start();
            NOX::LAPACK::Vector lXTime = lXTimeAll[iIntPoint];
//             lFftInvTimeTotal += lFftInvTime.Stop();
            if (iIntPoint == 0 && iLoop == 0) lXTimePrev = lXTimeAvg;
            
            // calculate the nonlinearity jacobian in time domain
            NOX::LAPACK::Matrix<double> lNonlin = ComputeJacobianTimeDomain(lXTime, lXTimePrev);
            
            if (IsCorrectingX())
            {
                // we use the correction from the F evaluation because in that we use the actual lXTime instead of some "fictional variations"
                // like in finite difference.
                FResult lFResult = ComputeFTimeDomain(lXTime, lXTimePrev);
                if (lFResult.XCorrSet) lXTimePrev = lFResult.XCorr;
                else lXTimePrev = lXTime;
            }
            else lXTimePrev = lXTime;
            
            if (iLoop == lLoopCount - 1)
            {
                lJTimeAll.push_back(lNonlin);
//                 for (int iHarm = 0; iHarm < mHarmonicCoeffCount; iHarm++)
//                 {
//                     for (int jHarm = 0; jHarm < mHarmonicCoeffCount; jHarm++)
//                     {
//                         double lBProdIndex = GetBProductIndex(iHarm, jHarm, iIntPoint, mHarmonicCoeffCount, lIntPoints.size());
//                         double lBProdValue = lBProducts[lBProdIndex];
//                         
//                         // we iterate over only the nonlinear dofs, because the rest of values in the lNonlinResult.FValues will be zeros
//                         for (int iDof = 0; iDof < mNonzeroFPositions.size(); iDof++)
//                         {
//                             int lDof = mNonzeroFPositions[iDof];
//                             int lHarmInd1 = GetHBMDofIndex(iDof, iHarm, mHarmonicCoeffCount);
//                             
//                             for (int jDof = 0; jDof < lDofCount; jDof++)
//                             {
//                                 int lHarmInd2 = GetHBMDofIndex(jDof, jHarm, mHarmonicCoeffCount);
//                                 
//                                 lReturnMatrix(lHarmInd1, lHarmInd2) += lBProdValue * lNonlin(lDof, jDof);
//                             }
//                         }
//                     }
//                 }
            }
        }
    }
    
    const NOX::LAPACK::Matrix<double>& lReturnMatrix = lAft.TimeToFrequency(lJTimeAll, aFrequency);
    
    return lReturnMatrix;
}
void NonlinearBase::Finalise()
{
    if (!mIsInitialised) throw "This nonlinearity has not yet been initialised. Initialise first before calling finalise!";
    if (mIsFinalised) return;
    mIsFinalised = true;
    
    // this is the final time this gets executed
    mNonzeroFPositions = NonzeroFPositions();
    
    // check dof validity
    for (int iDof = 0; iDof < mNonzeroFPositions.size(); iDof++)
    {
        int lDof = mNonzeroFPositions[iDof];
        if (lDof < 0) throw "Dof index can not be negative!";
        // we can do this check here because at this point the object is already initialised
        if (lDof >= mDofCountTimeDomain) throw "Dof index (" + std::to_string(lDof) + ") exceeds the size of the problem! (" + std::to_string(mDofCountTimeDomain) + ")";
    }
}
bool NonlinearBase::IsFinalised() const
{
    return mIsFinalised;
}
int NonlinearBase::DofCountTimeDomain() const
{
    if (!mIsInitialised)
        throw "Dof count is not set yet because the nonlinearity is not initialised!";
    
    return mDofCountTimeDomain;
}

NOX::LAPACK::Matrix<double> NonlinearBase::ComputeJacobianFiniteDifference(const NOX::LAPACK::Vector& aX, const std::function<NOX::LAPACK::Vector(const NOX::LAPACK::Vector&)>& aFEval, double aStep) const
{
    if (aStep <= 0) throw "Step must be a positive value!";
    
    int lSize = aX.length();
    NOX::LAPACK::Matrix<double> lSteps(lSize, lSize);
    
    for (int j = 0; j < lSize; j++)
        for (int i = 0; i < lSize; i++) // loop order switched because the matrix is column major
            lSteps(i, j) = aStep;
        
    return ComputeJacobianFiniteDifference(aX, aFEval, lSteps);
}

NOX::LAPACK::Matrix<double> NonlinearBase::ComputeJacobianFiniteDifference(const NOX::LAPACK::Vector& aX, const std::function<NOX::LAPACK::Vector(const NOX::LAPACK::Vector&)>& aFEval, const NOX::LAPACK::Matrix<double>& aSteps) const
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
            
            NOX::LAPACK::Vector lF1 = aFEval(lVec1);
            NOX::LAPACK::Vector lF2 = aFEval(lVec2);
            
            double lDerivative = (lF1(i) - lF2(i)) / 2.0 / lStep;
            
            lReturnMatrix(i, j) = lDerivative;
        }
    }
    
    return lReturnMatrix;
}

NOX::LAPACK::Matrix<double> NonlinearBase::ComputeJacobianFiniteDifference(const NOX::LAPACK::Vector& aX, const double& aFrequency, double aStep) const
{
    std::function<NOX::LAPACK::Vector(const NOX::LAPACK::Vector&)> lFEval = [aFrequency, this](const NOX::LAPACK::Vector& aIn) { return ComputeF(aIn, aFrequency); };
    
    return ComputeJacobianFiniteDifference(aX, lFEval, aStep);
}

NOX::LAPACK::Matrix<double> NonlinearBase::ComputeJacobianFiniteDifference(const NOX::LAPACK::Vector& aX, const double& aFrequency, const NOX::LAPACK::Matrix<double>& aSteps) const
{
    std::function<NOX::LAPACK::Vector(const NOX::LAPACK::Vector&)> lFEval = [aFrequency, this](const NOX::LAPACK::Vector& aIn) { return ComputeF(aIn, aFrequency); };
    
    return ComputeJacobianFiniteDifference(aX, lFEval, aSteps);
}
// Same functions as above, but using the ComputeF as the aFEval function (in time domain)
NOX::LAPACK::Matrix<double> NonlinearBase::ComputeJacobianFiniteDifferenceTD(const NOX::LAPACK::Vector& aX, const NOX::LAPACK::Vector& aXPrev, double aStep) const
{
    std::function<NOX::LAPACK::Vector(const NOX::LAPACK::Vector&)> lFEval = [aXPrev, this](const NOX::LAPACK::Vector& aIn) 
    {
        FResult lRes = ComputeFTimeDomain(aIn, aXPrev);
        return lRes.FValues;
    };
    
    return ComputeJacobianFiniteDifference(aX, lFEval, aStep);
}
NOX::LAPACK::Matrix<double> NonlinearBase::ComputeJacobianFiniteDifferenceTD(const NOX::LAPACK::Vector& aX, const NOX::LAPACK::Vector& aXPrev, const NOX::LAPACK::Matrix<double>& aSteps) const
{
    std::function<NOX::LAPACK::Vector(const NOX::LAPACK::Vector&)> lFEval = [aXPrev, this](const NOX::LAPACK::Vector& aIn) 
    {
        FResult lRes = ComputeFTimeDomain(aIn, aXPrev);
        return lRes.FValues;
    };
    
    return ComputeJacobianFiniteDifference(aX, lFEval, aSteps);
}

void NonlinearBase::CheckStatus() const
{
    if (!mIsInitialised)    throw "Nonlinearity is not initialised!";
    if (!mIsFinalised)      throw "Nonlinearity is not finalised!";
}
