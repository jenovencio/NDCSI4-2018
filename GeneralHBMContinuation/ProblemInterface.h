
class ProblemInterface;

#pragma once
#include "LOCA_LAPACK_Interface.H"
#include "LOCA_Parameter_Vector.H"
#include "NOX_LAPACK_Matrix.H"
#include "NOX_LAPACK_Vector.H"
#include "Teuchos_RCP.hpp"

#include "Config.h"
#include "ProblemParams.h"
#include "Nonlinearities/NonlinearBase.h"
#include "Aft/AftBase.h"
#include "DSBuilder/DSBuilderBase.h"

class ProblemInterface : public LOCA::LAPACK::Interface
{
private:
    const Config cConfig;
    const ProblemParams cProblemParams;
    
    double mFrequency = 0;
    double mCurrentContParam = -1;
    bool mRecomputeDynamicStiffness = true;
    bool mRecomputeExcitationRHS = true;
    // current dynamic matrix (for the last used frequency)
    NOX::LAPACK::Matrix<double>* mDynamicStiffnessMatrix = nullptr;
    // current excitation RHS (for the last used frequency)
    NOX::LAPACK::Vector* mExcitationRHS = nullptr;
    
    // we can have either one matrix for all harmonic or one matrix for each harmonic (nothing in between)
//     std::vector<NOX::LAPACK::Matrix<double>> mMassMatrices;
//     std::vector<NOX::LAPACK::Matrix<double>> mDampingMatrices;
//     std::vector<NOX::LAPACK::Matrix<double>> mStiffnessMatrices;
    // amplitudes of the excitation force, for each dof and harmonic wave (DC, cos, sin)
    std::vector<double> mExcitationAmp;
    // number of excitation coeffs for each physical dof (can be less than the number of coeffs used for the HBM)
    int mExcitationCoeffCount;
    
    NOX::LAPACK::Vector mInitGuess;
    
    std::vector<double> mSolutionFrequencies;
    std::vector<double> mSolutionNorms;
    std::vector<NOX::LAPACK::Vector> mSolutions;
    
    // nonlinearities coming from the "outside" (from the constructor)
    // defined in the outer code
    // we won't delete this (the pointers inside the vector),
    // it's up to the outside caller to delete the pointers
    const std::vector<NonlinearBase*> cNonlinearitiesOuter;
    // nonlinearities from the "inside", loaded from files specified in the provided config
    // we will delete these pointers because they will be created inside this class
    std::vector<NonlinearBase*> mNonlinearitiesInner;
    // all nonlinearities, outer and inner combined
    std::vector<NonlinearBase*> mNonlinearities;
    
    bool mHasWholeSolutions = false;
    
    // aft implementation
    // delete this pointer
    AftBase* mAft = nullptr;
    // dyn stiff matrix creation implementation
    // delete this pointer
    DSBuilderBase* mDSBuilder = nullptr;
    
public:
    // name of the generic continuation parameter that will go between 0 and 1
    static const std::string cContParameterName;
    
public:
    ProblemInterface(const Config& aConfig, const std::vector<NonlinearBase*>& aNonlinearities);
    ~ProblemInterface();
    
    virtual const NOX::LAPACK::Vector& getInitialGuess() override;
    virtual bool computeF(NOX::LAPACK::Vector& aRhs, const NOX::LAPACK::Vector& aX) override;
    virtual bool computeJacobian(NOX::LAPACK::Matrix<double>& aJ, const NOX::LAPACK::Vector& aX) override;
    virtual void setParams(const LOCA::ParameterVector& aParams) override;
    // this one is used to actually store the solutions
    virtual void printSolution(const NOX::LAPACK::Vector& aX, const double aConParam) override;
    
    void ClearSolutions();
    void WriteSolutionNorms(std::ostream& aStream) const;
    void WriteWholeSolutions(std::ostream& aStream) const;
    void PlotSolutionNorms() const;
    bool HasWholeSolutions() const;
    
private:
    NOX::LAPACK::Matrix<double>* CreateDynamicStiffnessMatrix(double aFrequency);
    NOX::LAPACK::Vector* CreateExcitationRHS(double aFrequency);
    
    void CheckMatrixCount(const std::string& aMatrixName, const int& aMatrixCount, const int& aHarmonicWaveCount);
    void CheckMatrixSize(const NOX::LAPACK::Matrix<double>& aMatrix, const std::string& aMatrixName);
};
