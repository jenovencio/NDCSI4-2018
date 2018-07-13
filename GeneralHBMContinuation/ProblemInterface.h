
class ProblemInterface;

#pragma once
#include "LOCA_LAPACK_Interface.H"
#include "LOCA_Parameter_Vector.H"
#include "NOX_LAPACK_Matrix.H"
#include "NOX_LAPACK_Vector.H"
#include "Teuchos_RCP.hpp"

#include "Config.h"

class ProblemInterface : public LOCA::LAPACK::Interface
{
private:
    // number of standard degrees of freedom (in physical domain, not HBM)
    const int cDOFCount;
    // total number of harmonic waves
    const int cHarmonicCount;
    // in other words, the total number of degrees of freedom in the HBM system
    const int cDOFCountHBM;
    double mFrequency = 0;
    bool mFrequencyChanged = false;
    // coordinates of the integration points, relative to a time period T
    std::vector<double> mIntPointsRelative;
    
    double* mBProducts = nullptr;
    
    NOX::LAPACK::Matrix<double> mMassMatrix;
    NOX::LAPACK::Matrix<double> mDampingMatrix;
    NOX::LAPACK::Matrix<double> mStiffnessMatrix;
    
    NOX::LAPACK::Vector mInitGuess;
    
public:
    ProblemInterface(const Config& aConfig);
    ~ProblemInterface();
    
    virtual const NOX::LAPACK::Vector& getInitialGuess() override;
    virtual bool computeF(NOX::LAPACK::Vector& aRhs, const NOX::LAPACK::Vector& aX) override;
    virtual bool computeJacobian(NOX::LAPACK::Matrix<double>& aJ, const NOX::LAPACK::Vector& aX) override;
    virtual void setParams(const LOCA::ParameterVector& aParams) override;
    
private:
    NOX::LAPACK::Matrix<double> CreateDynamicStiffnessMatrix(double aFrequency);
};
