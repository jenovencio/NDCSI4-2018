
#include "NonlinearBaseFD.h"

void NonlinearBaseFD::Init(AftBase* const aAft, const ProblemParams& aProblemParams)
{
    NonlinearBase::Init(aAft, aProblemParams);
    
    mFiniteDifferenceSteps = NOX::LAPACK::Matrix<double>(aProblemParams.DofCountPhysical, aProblemParams.DofCountPhysical);
    
    SetFiniteDifferenceSteps(DEFAULT_FD_STEP);
}

NOX::LAPACK::Matrix<double> NonlinearBaseFD::ComputeJacobianTimeDomain(const NOX::LAPACK::Vector& aX, const NOX::LAPACK::Vector& aXPrev) const
{
    return ComputeJacobianFiniteDifferenceTD(aX, aXPrev, mFiniteDifferenceSteps);
}
void NonlinearBaseFD::SetFiniteDifferenceSteps(const double& aStep)
{
    for (int j = 0; j < mFiniteDifferenceSteps.numCols(); j++)
        for (int i = 0; i < mFiniteDifferenceSteps.numRows(); i++)
            mFiniteDifferenceSteps(i, j) = aStep;
}
void NonlinearBaseFD::SetFiniteDifferenceSteps(const NOX::LAPACK::Matrix<double>& aSteps)
{
    if (mFiniteDifferenceSteps.numRows() != aSteps.numRows()) 
        throw "Input matrix has different number of rows (" + std::to_string(aSteps.numRows()) + ") that the step matrix (" + std::to_string(aSteps.numRows()) + ")!";
    
    if (mFiniteDifferenceSteps.numCols() != aSteps.numCols()) 
        throw "Input matrix has different number of columns (" + std::to_string(aSteps.numCols()) + ") that the step matrix (" + std::to_string(aSteps.numCols()) + ")!";
    
    for (int j = 0; j < mFiniteDifferenceSteps.numCols(); j++)
        for (int i = 0; i < mFiniteDifferenceSteps.numRows(); i++)
            mFiniteDifferenceSteps(i, j) = aSteps(i, j);
}
