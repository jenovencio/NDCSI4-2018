
struct FResult;
struct JResult;
class NonlinearBase;

#pragma once
#include <functional>

#include "NOX_LAPACK_Matrix.H"
#include "NOX_LAPACK_Vector.H"

#define DEFAULT_FD_STEP 1e-8

class NonlinearBase
{
    // LIFTIME of this class:
    // 1) create
    // 2) initialise
    // 3) make modifications
    // 4) finalise
    // 5) use, don't make any more modifications
    
private:
    // indices of nonlinear dofs in time domain, not to be touched outside the Init and Finalise methods,
    // 
    std::vector<int> mNonzeroFPositions;
    
    // do not delete these pointers
    const std::vector<double>* cIntegrationPoints;
    const std::vector<double>* cBValues;
    const std::vector<double>* cBProducts;
    int mHarmonicCoeffCount;
    bool mIsFinalised = false;
    bool mIsInitialised = false;
    
    // this is set in the Init call, do not touch it elswhere!
    int mDofCountTimeDomain;
    
public:
    virtual void LoadFromFile(const std::string& aFilePath);
    // returns the name of the class
    virtual std::string ClassName() const = 0;
    // the aDofCountTimeDomain argument was added for the derived classes, in case they need it. We don't really need it here
    virtual void Init(const std::vector<double>& aIntegrationPoints, const std::vector<double>& aBValues, const std::vector<double>& aBProducts, const int& aHarmonicCoeffCount, const int& aDofCountTimeDomain);
    // frequency domain to frequency domain
    NOX::LAPACK::Vector ComputeF(const NOX::LAPACK::Vector& aX, const double& aFrequency) const;
    // frequency domain to frequency domain
    NOX::LAPACK::Matrix<double> ComputeJacobian(const NOX::LAPACK::Vector& aX, const double& aFrequency) const;
    // "locks" the nonlinearity to any changes (logically, not actually of course)
    // the point is, after this method is called, the nonlinearity should not change anymore (in terms of it's mathematical definition at least),
    // so for instance for a cubic spring nonlinearity, no springs should be added or removed, etc.
    // changes to the nonlinearity after this function call might not be taken into consideration when the nonlinearity is not eveluated, that's 
    // the point of this mechanism
    void Finalise();
    bool IsFinalised() const;
    int DofCountTimeDomain() const;
    
protected:
    // time domain to time domain
    virtual FResult ComputeFTimeDomain(const NOX::LAPACK::Vector& aX, const NOX::LAPACK::Vector& aXPrev) const = 0;
    // time domain to time domain
    virtual NOX::LAPACK::Matrix<double> ComputeJacobianTimeDomain(const NOX::LAPACK::Vector& aX, const NOX::LAPACK::Vector& aXPrev) const = 0;
    
    virtual int NumberOfPrepLoops() const = 0;
    // returns indices of elements in the F vector (in time domain) that are nonzero, i.e.
    // where some nonlinearity occurs
    // this is just to speed up the F and jacobian computations in the frequency domain (so we don't uselessly iterate over all elements)
    virtual std::vector<int> NonzeroFPositions() const = 0;
    
    // general jacobian evaluation functions, using arbitrary F evaluation functions
    NOX::LAPACK::Matrix<double> ComputeJacobianFiniteDifference(const NOX::LAPACK::Vector& aX, const std::function<NOX::LAPACK::Vector(const NOX::LAPACK::Vector&)>& aFEval, double aStep = DEFAULT_FD_STEP) const;
    NOX::LAPACK::Matrix<double> ComputeJacobianFiniteDifference(const NOX::LAPACK::Vector& aX, const std::function<NOX::LAPACK::Vector(const NOX::LAPACK::Vector&)>& aFEval, const NOX::LAPACK::Matrix<double>& aSteps) const;
    
    // Same functions as above, but using the ComputeF as the aFEval function (in frequency domain)
    NOX::LAPACK::Matrix<double> ComputeJacobianFiniteDifference(const NOX::LAPACK::Vector& aX, const double& aFrequency, double aStep = DEFAULT_FD_STEP) const;
    NOX::LAPACK::Matrix<double> ComputeJacobianFiniteDifference(const NOX::LAPACK::Vector& aX, const double& aFrequency, const NOX::LAPACK::Matrix<double>& aSteps) const;
    
    // Same functions as above, but using the ComputeF as the aFEval function (in time domain)
    NOX::LAPACK::Matrix<double> ComputeJacobianFiniteDifferenceTD(const NOX::LAPACK::Vector& aX, const NOX::LAPACK::Vector& aXPrev, double aStep = DEFAULT_FD_STEP) const;
    NOX::LAPACK::Matrix<double> ComputeJacobianFiniteDifferenceTD(const NOX::LAPACK::Vector& aX, const NOX::LAPACK::Vector& aXPrev, const NOX::LAPACK::Matrix<double>& aSteps) const;
    
    virtual bool IsCorrectingX() const = 0;
    
private:
    // relative time point value (considering period = 1)
    NOX::LAPACK::Vector FreqToTime(const NOX::LAPACK::Vector& aX, const int& aIntegrationPointIndex) const;
    
    void CheckStatus() const;
};

struct FResult
{
public:
    NOX::LAPACK::Vector FValues;
    // correction of the provided X
    NOX::LAPACK::Vector XCorr;
    bool XCorrSet = false;
};

