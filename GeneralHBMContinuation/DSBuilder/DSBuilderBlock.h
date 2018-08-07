
class DSBuilderBlock;

#pragma once
#include "DSBuilderBase.h"

class DSBuilderBlock : public DSBuilderBase
{
private:
    // for the 0th harmonic there will only be one stiffness matrix (size dofphys x dofphys)
    NOX::LAPACK::Matrix<double> mStiffnessMatrixDC;
    // matrices for the nonzero harmonics, these will be "2x2" block matrices (i.e. size (2*dofphys) x (2*dofphys))
    std::vector<NOX::LAPACK::Matrix<double>> mMassMatricesDyn;
    // damping matrices are optional - if this vector is empty, damping matrices will just not be used 
    // (it's not an error unline in case of stiffness or mass matrices)
    std::vector<NOX::LAPACK::Matrix<double>> mDampingMatricesDyn;
    std::vector<NOX::LAPACK::Matrix<double>> mStiffnessMatricesDyn;
    
    // just a zero matrix of a proper size ((2*dofphys) x (2*dofphys)) to which we redirect
    // the ourselves when we don't have actual damping matrices
    NOX::LAPACK::Matrix<double> mDummyDamping;
    
    int mHarmonicWaveCount = 0;
    
public:
    virtual void Init(const Config& aConfig, const ProblemParams& aParams) override;
    
protected:
    virtual NOX::LAPACK::Matrix<double>* CreateDynamicStiffnessMatrixInner(const double& aFrequency) override;
    
};
