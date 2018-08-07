
#include "NonlinearBase.h"

void NonlinearBase::LoadFromFile(const std::string& aFilePath)
{
    std::string lTypeName = ClassName();
    
    throw "Can not load a nonlinearity from file \"" + aFilePath + "\"! This class does not support loading from files.";
}

void NonlinearBase::Init(AftBase* const aAft, const ProblemParams& aProblemParams)
{
    if (mIsInitialised) throw "This nonlinearity has already been initialised, can not initialise again!";
    if (mIsFinalised) throw "This nonlinearity has already been finalised! Fix the way you treat the object. It first needs to be initialised, then finalised!";
    
    cAft = aAft;
    mProblemParams = aProblemParams;
        
    mIsInitialised = true;
}
// frequency domain to frequency domain
NOX::LAPACK::Vector NonlinearBase::ComputeF(const NOX::LAPACK::Vector& aX, const double& aFrequency)
{
    CheckStatus();
    
    return ComputeFInner(aX, aFrequency);
}
// frequency domain to frequency domain
NOX::LAPACK::Matrix<double> NonlinearBase::ComputeJacobian(const NOX::LAPACK::Vector& aX, const double& aFrequency)
{
    CheckStatus();
    
    return ComputeJacobianInner(aX, aFrequency);
}

void NonlinearBase::Finalise()
{
    if (!mIsInitialised) throw "This nonlinearity has not yet been initialised. Initialise first before calling finalise!";
    if (mIsFinalised) return;
    mIsFinalised = true;
}
bool NonlinearBase::IsFinalised() const
{
    return mIsFinalised;
}
int NonlinearBase::DofCountTimeDomain() const
{
    if (!mIsInitialised)
        throw "Dof count is not set yet because the nonlinearity is not initialised!";
    
    return mProblemParams.DofCountPhysical;
}

void NonlinearBase::CheckStatus() const
{
    if (!mIsInitialised)    throw "Nonlinearity is not initialised!";
    if (!mIsFinalised)      throw "Nonlinearity is not finalised!";
}
