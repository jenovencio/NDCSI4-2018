
class NonlinearBase;

#pragma once
#include <string>
#include <functional>

#include "NOX_LAPACK_Matrix.H"
#include "NOX_LAPACK_Vector.H"

#include "../Aft/AftBase.h"
#include "../ProblemParams.h"

class NonlinearBase
{
    // LIFTIME of this class:
    // 1) create
    // 2) initialise
    // 3) make modifications
    // 4) finalise
    // 5) use, don't make any more modifications
    
protected:
    
    // do not delete this pointer
    // can not be const because it will be modifying itself inside (it's internal structures) during the transformations
    AftBase* cAft;
    
    // this is set in the Init call, do not modify it elsewhere!
    ProblemParams mProblemParams;
    
private:
    bool mIsFinalised = false;
    bool mIsInitialised = false;
    
protected:
    // frequency domain to frequency domain
    virtual NOX::LAPACK::Vector ComputeFInner(const NOX::LAPACK::Vector& aX, const double& aFrequency) = 0;
    // frequency domain to frequency domain
    virtual NOX::LAPACK::Matrix<double> ComputeJacobianInner(const NOX::LAPACK::Vector& aX, const double& aFrequency) = 0;
public:
    virtual ~NonlinearBase() { }
    virtual void LoadFromFile(const std::string& aFilePath);
    // returns the name of the class
    virtual std::string ClassName() const = 0;
    virtual void Init(AftBase* const aAft, const ProblemParams& aProblemParams);
    // frequency domain to frequency domain
    NOX::LAPACK::Vector ComputeF(const NOX::LAPACK::Vector& aX, const double& aFrequency);
    // frequency domain to frequency domain
    NOX::LAPACK::Matrix<double> ComputeJacobian(const NOX::LAPACK::Vector& aX, const double& aFrequency);
    // "locks" the nonlinearity to any changes (logically, not actually of course)
    // the point is, after this method is called, the nonlinearity should not change anymore (in terms of it's mathematical definition at least),
    // so for instance for a cubic spring nonlinearity, no springs should be added or removed, etc.
    // changes to the nonlinearity after this function call might not be taken into consideration when the nonlinearity is not eveluated, that's 
    // the point of this mechanism
    virtual void Finalise();
    bool IsFinalised() const;
    int DofCountTimeDomain() const;
    
private:
    
    void CheckStatus() const;
};
