
class ProblemInterface;

#pragma once
#include "LOCA_LAPACK_Interface.H"
#include "LOCA_Parameter_Vector.H"
#include "NOX_LAPACK_Matrix.H"
#include "NOX_LAPACK_Vector.H"
#include "Teuchos_RCP.hpp"

#include "Config.h"
#include "Nonlinearities/NonlinearBase.h"

class ProblemInterface : public LOCA::LAPACK::Interface
{
private:
    // number of standard degrees of freedom (in physical domain, not HBM)
    // not constant because it's determined from the system files (mass, stiffness, damping matrix)
    int mDOFCount;
    // total number of harmonic waves
    const int cHarmonicCount;
    // in other words, the total number of degrees of freedom in the HBM system
    int mDOFCountHBM;
    double mFrequency = 0;
    bool mRecomputeDynamicStiffness = true;
    bool mRecomputeExcitationRHS = true;
    // current dynamic matrix (for the last used frequency)
    NOX::LAPACK::Matrix<double>* mDynamicStiffnessMatrix = nullptr;
    // current excitation RHS (for the last used frequency)
    NOX::LAPACK::Vector* mExcitationRHS = nullptr;
    // coordinates of the integration points, relative to a time period T
    // these values would effectively be for period T = 1
    std::vector<double> mIntPointsRelative;
    
    std::vector<double> mBValues;
    std::vector<double> mBProducts;
    
    NOX::LAPACK::Matrix<double> mMassMatrix;
    NOX::LAPACK::Matrix<double> mDampingMatrix;
    NOX::LAPACK::Matrix<double> mStiffnessMatrix;
    // amplitudes of the excitation force, for each dof and harmonic wave (DC, cos, sin)
    std::vector<double> mExcitationAmp;
    // number of excitation coeffs for each physical dof (can be less than the number of coeffs used for the HBM)
    int mExcitationCoeffCount;
    
    NOX::LAPACK::Vector mInitGuess;
    
    std::vector<double> mSolutionFrequencies;
    std::vector<double> mSolutionNorms;
    std::vector<NOX::LAPACK::Vector> mSolutions;
    
    // we won't delete this, it's up to the outside caller to delete the pointers inside the vector
    const std::vector<NonlinearBase*> cNonlinearities;
    
    const bool cSaveWholeSolutions;
    
public:
    static const std::string cFrequencyName;
    
public:
    ProblemInterface(const Config& aConfig, const std::vector<NonlinearBase*>& aNonlinearities, bool aSaveWholeSolutions = false);
    ~ProblemInterface();
    
    virtual const NOX::LAPACK::Vector& getInitialGuess() override;
    virtual bool computeF(NOX::LAPACK::Vector& aRhs, const NOX::LAPACK::Vector& aX) override;
    virtual bool computeJacobian(NOX::LAPACK::Matrix<double>& aJ, const NOX::LAPACK::Vector& aX) override;
    virtual void setParams(const LOCA::ParameterVector& aParams) override;
    // this one is used to actually store the solutions
    virtual void printSolution(const NOX::LAPACK::Vector& aX, const double aConParam) override;
    
    void ClearSolutions();
    void WriteSolutionNorms(std::ostream& aStream);
    void WriteWholeSolutions(std::ostream& aStream);
    
private:
    NOX::LAPACK::Matrix<double>* CreateDynamicStiffnessMatrix(double aFrequency);
    NOX::LAPACK::Vector* CreateExcitationRHS(double aFrequency);
};
