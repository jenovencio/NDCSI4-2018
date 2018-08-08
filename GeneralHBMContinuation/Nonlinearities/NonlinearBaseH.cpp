
#include "NonlinearBaseH.h"
#include "../Functions.h"

// frequency domain to frequency domain
NOX::LAPACK::Vector NonlinearBaseH::ComputeFInner(const NOX::LAPACK::Vector& aX, const double& aFrequency) 
{
    AftBase& lAft = *cAft;
    
// check 
//     if (aFrequency <= 0) throw "Frequency must be a positive value!";
    if (mProblemParams.HarmonicCount <= 0) throw "Number of harmonic coefficients must be a positive integer!";
    if (aX.length() % mProblemParams.HarmonicCount != 0) throw "Size of the problem in frequency domain (" + std::to_string(aX.length()) + ") is not divisible by number of harmonic coefficients (" + std::to_string(mProblemParams.HarmonicCount) + ")!";
    
    // in time domain
    int lDofCount = mProblemParams.DofCountPhysical;
    
    std::vector<NOX::LAPACK::Vector> lXTimeAll = lAft.FrequencyToTime(aX, aFrequency);
    
    std::vector<NOX::LAPACK::Vector> lFTimeAll;
    lFTimeAll.reserve(lXTimeAll.size());
    
    NOX::LAPACK::Vector lXTimePrev(lDofCount);
    
    
    int lLoopCount = NumberOfPrepLoops() + 1; // number of loop over the period
    
    NOX::LAPACK::Vector lXTimeAvg(lDofCount);
        
    for (int iDof = 0; iDof < lDofCount; iDof++)
    {
        int lHarmIndex = GetHBMDofIndex(iDof, 0, mProblemParams.HarmonicCount);
        lXTimeAvg(iDof) = aX(lHarmIndex);
    }
//     std::cout << "XTimeAvg: " << lXTimeAvg << std::endl;
    
    InitFComputation(aX);
    
    for (int iLoop = 0; iLoop < lLoopCount; iLoop++)
    {
        for (int iTimePoint = 0; iTimePoint < lXTimeAll.size(); iTimePoint++)
        {
            const NOX::LAPACK::Vector& lXTime = lXTimeAll[iTimePoint];
            if (iTimePoint == 0 && iLoop == 0) lXTimePrev = lXTimeAvg;
            
            // calculate the nonlinearity in time domain
            NOX::LAPACK::Vector lNonlinResult = ComputeFTD(lXTime, iTimePoint);
                        
            if (iLoop == lLoopCount - 1)
            {
                lFTimeAll.push_back(lNonlinResult);
            }
        }
    }
    
    const NOX::LAPACK::Vector& lReturnVector = lAft.TimeToFrequency(lFTimeAll, aFrequency);
    
//     double lMs = lTime.Stop();
//     std::cout << "Full F eval time: " << lMs << " ms" << std::endl;
    
    return lReturnVector;
}
// frequency domain to frequency domain
NOX::LAPACK::Matrix<double> NonlinearBaseH::ComputeJacobianInner(const NOX::LAPACK::Vector& aX, const double& aFrequency)
{
    AftBase& lAft = *cAft;
    
    // check
//     if (aFrequency <= 0) throw "Frequency must be a positive value!";
    if (mProblemParams.HarmonicCount <= 0) throw "Number of harmonic coefficients must be a positive integer!";
    if (aX.length() % mProblemParams.HarmonicCount != 0) throw "Size of the problem in frequency domain (" + std::to_string(aX.length()) + ") is not divisible by number of harmonic coefficients (" + std::to_string(mProblemParams.HarmonicCount) + ")!";
        
    const std::vector<NOX::LAPACK::Vector>& lXTimeAll = lAft.FrequencyToTime(aX, aFrequency);
    
    std::vector<NOX::LAPACK::Matrix<double>> lJTimeAll;
    
    lJTimeAll.reserve(lXTimeAll.size());
        
    int lLoopCount = NumberOfPrepLoops() + 1; // number of loop over the period
    
    InitJComputation(aX);
    
    for (int iLoop = 0; iLoop < lLoopCount; iLoop++)
    {
        for (int iIntPoint = 0; iIntPoint < lXTimeAll.size(); iIntPoint++)
        {
//             lFftInvTime.Start();
            const NOX::LAPACK::Vector& lXTime = lXTimeAll[iIntPoint];
            
            // calculate the nonlinearity jacobian in time domain (derivatives by fourier coefficients)
            NOX::LAPACK::Matrix<double> lNonlin = ComputeDFDH(lXTime, iIntPoint);
            
            if (iLoop == lLoopCount - 1)
            {
                lJTimeAll.push_back(lNonlin);
            }
        }
    }
    
    const NOX::LAPACK::Matrix<double>& lReturnMatrix = lAft.TimeToFrequency(lJTimeAll, aFrequency);
    
    return lReturnMatrix;
}

