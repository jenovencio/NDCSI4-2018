
class NonlinearBase;

#pragma once
#include <functional>

#include "NOX_LAPACK_Matrix.H"
#include "NOX_LAPACK_Vector.H"

class NonlinearBase
{
protected:
    // do not delete these pointers
    const std::vector<double>* cIntegrationPoints;
    const std::vector<double>* cBValues;
    const std::vector<double>* cBProducts;
    int mHarmonicCoeffCount;
    
public:
    void Init(const std::vector<double>& aIntegrationPoints, const std::vector<double>& aBValues, const std::vector<double>& aBProducts, const int& aHarmonicCoeffCount);
    // frequency domain to frequency domain
    virtual NOX::LAPACK::Vector ComputeRHS(const NOX::LAPACK::Vector& aX, const double& aFrequency) = 0;
    // frequency domain to frequency domain
    virtual NOX::LAPACK::Matrix<double> ComputeJacobian(const NOX::LAPACK::Vector& aX, const double& aFrequency) = 0;
    
protected:
    NOX::LAPACK::Matrix<double> ComputeJacobianFiniteDifference(const NOX::LAPACK::Vector& aX, const std::function<NOX::LAPACK::Vector(const NOX::LAPACK::Vector&)>& aRHSEval, double aStep = 1e-5);
    NOX::LAPACK::Matrix<double> ComputeJacobianFiniteDifference(const NOX::LAPACK::Vector& aX, const std::function<NOX::LAPACK::Vector(const NOX::LAPACK::Vector&)>& aRHSEval, const NOX::LAPACK::Matrix<double>& aSteps);
    
    // relative time point value (considering period = 1)
    NOX::LAPACK::Vector FreqToTime(const NOX::LAPACK::Vector& aX, const int& aIntegrationPointIndex);
};
