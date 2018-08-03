
#include "NonlinearBaseH.h"
#include "../Functions.h"

// frequency domain to frequency domain
NOX::LAPACK::Vector NonlinearBaseH::ComputeFInner(const NOX::LAPACK::Vector& aX, const double& aFrequency) const 
{
    AftBase& lAft = *cAft;
    
// check 
//     if (aFrequency <= 0) throw "Frequency must be a positive value!";
    if (mProblemParams->HarmonicCount <= 0) throw "Number of harmonic coefficients must be a positive integer!";
    if (aX.length() % mProblemParams->HarmonicCount != 0) throw "Size of the problem in frequency domain (" + std::to_string(aX.length()) + ") is not divisible by number of harmonic coefficients (" + std::to_string(mProblemParams->HarmonicCount) + ")!";
    
    // in time domain
    int lDofCount = aX.length() / mProblemParams->HarmonicCount;
    
    std::vector<NOX::LAPACK::Vector> lXTimeAll = lAft.FrequencyToTime(aX, aFrequency);
    
    std::vector<NOX::LAPACK::Vector> lFTimeAll;
    lFTimeAll.reserve(lXTimeAll.size());
    
    NOX::LAPACK::Vector lXTimePrev(lDofCount);
    
    
    int lLoopCount = NumberOfPrepLoops() + 1; // number of loop over the period
    
    NOX::LAPACK::Vector lXTimeAvg(lDofCount);
        
    for (int iDof = 0; iDof < lDofCount; iDof++)
    {
        int lHarmIndex = GetHBMDofIndex(iDof, 0, mProblemParams->HarmonicCount);
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
            FResult lNonlinResult = ComputeFTD(lXTime, lXTimePrev, iTimePoint);
            
            if (!IsHistoryDependent() && lNonlinResult.XCorrSet)
                throw "The code logic of this class is wrong! The class says it's not history dependent, but at the same time sets corrections to history in it's F evaluations.";
            
            if (IsHistoryDependent() && lNonlinResult.XCorrSet)
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
NOX::LAPACK::Matrix<double> NonlinearBaseH::ComputeJacobianInner(const NOX::LAPACK::Vector& aX, const double& aFrequency) const
{
    AftBase& lAft = *cAft;
    
    // check
//     if (aFrequency <= 0) throw "Frequency must be a positive value!";
    if (mProblemParams->HarmonicCount <= 0) throw "Number of harmonic coefficients must be a positive integer!";
    if (aX.length() % mProblemParams->HarmonicCount != 0) throw "Size of the problem in frequency domain (" + std::to_string(aX.length()) + ") is not divisible by number of harmonic coefficients (" + std::to_string(mProblemParams->HarmonicCount) + ")!";
    
    // in time domain
    int lDofCount = aX.length() / mProblemParams->HarmonicCount;
    
    const std::vector<NOX::LAPACK::Vector>& lXTimeAll = lAft.FrequencyToTime(aX, aFrequency);
    
    std::vector<NOX::LAPACK::Matrix<double>> lJTimeAll;
    lJTimeAll.reserve(lXTimeAll.size());
    
    NOX::LAPACK::Vector lXTimePrev;
    NOX::LAPACK::Matrix<double> lJTimePrev = GetFirstJ();
    
    int lLoopCount = NumberOfPrepLoops() + 1; // number of loop over the period
    
    NOX::LAPACK::Vector lXTimeAvg(lDofCount);
    
    for (int iDof = 0; iDof < lDofCount; iDof++)
    {
        int lHarmIndex = GetHBMDofIndex(iDof, 0, mProblemParams->HarmonicCount);
        lXTimeAvg(iDof) = aX(lHarmIndex);
    }
//     std::cout << "XTimeAvg: " << lXTimeAvg << std::endl;
    
    for (int iLoop = 0; iLoop < lLoopCount; iLoop++)
    {
        for (int iIntPoint = 0; iIntPoint < lXTimeAll.size(); iIntPoint++)
        {
//             lFftInvTime.Start();
            const NOX::LAPACK::Vector& lXTime = lXTimeAll[iIntPoint];
//             lFftInvTimeTotal += lFftInvTime.Stop();
            if (iIntPoint == 0 && iLoop == 0) lXTimePrev = lXTimeAvg;
            
            // calculate the nonlinearity jacobian in time domain (derivatives by fourier coefficients)
            NOX::LAPACK::Matrix<double> lNonlin = ComputeDFDH(lXTime, lXTimePrev, lJTimePrev, iIntPoint);
            
            if (iLoop == lLoopCount - 1)
            {
                lJTimeAll.push_back(lNonlin);
            }
            
            if (IsHistoryDependent())
            {
                FResult lFResult = ComputeFTD(lXTime, lXTimePrev, iIntPoint);
                if (lFResult.XCorrSet) lXTimePrev = lFResult.XCorr;
                else lXTimePrev = lXTime;
                
                lJTimePrev = lNonlin;
            }
            else lXTimePrev = lXTime;
        }
    }
    
    const NOX::LAPACK::Matrix<double>& lReturnMatrix = lAft.TimeToFrequency(lJTimeAll, aFrequency);
    
    return lReturnMatrix;
}

