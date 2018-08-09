
class FricElem1D;
struct FricElem1DDef;

#pragma once
#include <set>
#include "NonlinearBaseTDHD.h"

class FricElem1D : public NonlinearBaseTDHD
{
private:
    
    std::set<int> mFrictionDofs;
    std::vector<FricElem1DDef> mFrictionDefinitions;
    
    NOX::LAPACK::Vector mXPrev;
    
protected:
    virtual FResult ComputeFAndC(const NOX::LAPACK::Vector& aX, const NOX::LAPACK::Vector& aC) const override;
    // 4 derivatives necessary to compute the DFDH for the history dependent nonlinearity
    // the aC parameter comes from the previous F evaluation (for previous time step)
    virtual JResult ComputeDerivatives(const NOX::LAPACK::Vector& aX, const NOX::LAPACK::Vector& aC) const override;
    virtual NOX::LAPACK::Vector InitC(const NOX::LAPACK::Vector& aX) const  override;
    
    virtual int NumberOfPrepLoops() const override { return 2; }
public:
    
    virtual void LoadFromFile(const std::string& aFilePath) override;
    virtual std::string ClassName() const override { return "Friction element 1D"; }
    
    void AddFrictionDef(const FricElem1DDef& aDef);
};


struct FricElem1DDef
{
public:
    // physical dof index
    int DofID;
    // normal load
    double N;
    // contact stiffness
    double K;
    // friction coefficient
    double Mu;
};
