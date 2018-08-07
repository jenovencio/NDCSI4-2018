
class DSBuilderSimple;

#pragma once
#include "DSBuilderBase.h"

class DSBuilderSimple : public DSBuilderBase
{
private:
    NOX::LAPACK::Matrix<double> mMassMatrix;
    NOX::LAPACK::Matrix<double> mDampingMatrix;
    NOX::LAPACK::Matrix<double> mStiffnessMatrix;
    
public:
    virtual void Init(const Config& aConfig, const ProblemParams& aParams) override;
    
protected:
    virtual NOX::LAPACK::Matrix<double>* CreateDynamicStiffnessMatrixInner(const double& aFrequency) override;
};
