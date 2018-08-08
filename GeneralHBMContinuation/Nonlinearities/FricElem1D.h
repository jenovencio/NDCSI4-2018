
class FricElem1D;

#pragma once
#include "NonlinearBaseTDHD.h"

class FricElem1D : public NonlinearBaseTDHD
{
private:
    int mDofID = 0;
    double mN0 = 0.001;
    double mMu = 0.2;
    double mKtx = 1;
    
    NOX::LAPACK::Vector mXPrev;
    
protected:
    virtual FResult ComputeFAndC(const NOX::LAPACK::Vector& aX, const NOX::LAPACK::Vector& aC) const override;
    // 4 derivatives necessary to compute the DFDH for the history dependent nonlinearity
    // the aC parameter comes from the previous F evaluation (for previous time step)
    virtual JResult ComputeDerivatives(const NOX::LAPACK::Vector& aX, const NOX::LAPACK::Vector& aC) const override;
    virtual NOX::LAPACK::Vector InitC(const NOX::LAPACK::Vector& aX) const  override;
    virtual NOX::LAPACK::Matrix<double> InitDGbyDH(const NOX::LAPACK::Vector& aX) const override;
    
    virtual int NumberOfPrepLoops() const override { return 1; }
public:
    
    virtual void LoadFromFile(const std::string& aFilePath) override { }
    virtual std::string ClassName() const override { return "Friction element 1D"; }
};
