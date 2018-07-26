
class FricElem3D;

#pragma once
#include "NonlinearBase.h"

class FricElem3D : public NonlinearBase
{
private:
protected:
    // time domain to time domain
    virtual FResult ComputeFTimeDomain(const NOX::LAPACK::Vector& aX, const NOX::LAPACK::Vector& aXPrev) const override;
    // time domain to time domain
    virtual NOX::LAPACK::Matrix<double> ComputeJacobianTimeDomain(const NOX::LAPACK::Vector& aX, const NOX::LAPACK::Vector& aXPrev) const override;
    virtual int NumberOfPrepLoops() const override;
    virtual std::vector<int> NonzeroFPositions() const override;
    virtual bool IsCorrectingX() const override;
public:
    virtual void LoadFromFile(const std::string& aFilePath) override;
    virtual std::string ClassName() const override;
};
