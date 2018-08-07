
class DSBuilderBase;

#pragma once
#include "NOX_LAPACK_Matrix.H"

#include "../ProblemParams.h"
#include "../MatrixDefinition.h"
#include "../Functions.h"
#include "../Misc.h"

class DSBuilderBase
{
private:
    bool mIsInitialised = false;
protected:
    // assigned in Init, do not modify anywhere else
    ProblemParams mProblemParams;
public:
    virtual ~DSBuilderBase() { }
    virtual void Init(const Config& aConfig, const ProblemParams& aParams)
    {
        mProblemParams = aParams;
        mIsInitialised = true;
    }
    
    NOX::LAPACK::Matrix<double>* CreateDynamicStiffnessMatrix(const double& aFrequency)
    {
        if (!mIsInitialised) throw "DS builder is not initialised!";
        
        return CreateDynamicStiffnessMatrixInner(aFrequency);
    }
    
protected:
    virtual NOX::LAPACK::Matrix<double>* CreateDynamicStiffnessMatrixInner(const double& aFrequency) = 0;
    
    // check size of a matrix (if it has certain number of rows and columns) (we consider just one size for both rows and columns)
    void CheckMatrixSize(const NOX::LAPACK::Matrix<double>& aMatrix, const int& aSize, const std::string& aMatrixName)
    {
        if (aMatrix.numCols() != aSize)
            throw "Number of columns of \"" + aMatrixName + "\" matrix does not equal to the number of dofs specified in the configuration!";
        
        if (aMatrix.numRows() != aSize)
            throw "Number of rows of \"" + aMatrixName + "\" matrix does not equal to the number of dofs specified in the configuration!";
    }
};
