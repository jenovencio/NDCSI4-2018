
#include "FunctionsAlgebra.h"

NOX::LAPACK::Matrix<double> Multiply(const NOX::LAPACK::Matrix<double>& aA, const NOX::LAPACK::Matrix<double>& aB)
{
    if (aA.numCols() != aB.numRows())
        throw "Incompatible matrix sizes!";
    
    NOX::LAPACK::Matrix<double> lReturnMatrix(aA.numRows(), aB.numCols());
    
    for (int iCol = 0; iCol < aB.numCols(); iCol++)
    {
        for (int iRow = 0; iRow < aA.numRows(); iRow++)
        {
            for (int iComp = 0; iComp < aA.numCols(); iComp++)
            {
                lReturnMatrix(iRow, iCol) += aA(iRow, iComp) * aB(iComp, iCol);
            }
        }
    }
    
    return lReturnMatrix;
}

NOX::LAPACK::Matrix<double> Add(const NOX::LAPACK::Matrix<double>& aA, const NOX::LAPACK::Matrix<double>& aB)
{
    if (aA.numRows() != aB.numRows() || aA.numCols() != aB.numCols())
        throw "Incompatible matrix sizes!";
    
    NOX::LAPACK::Matrix<double> lReturnMatrix(aA.numRows(), aA.numCols());
    
    for (int iCol = 0; iCol < aA.numCols(); iCol++)
    {
        for (int iRow = 0; iRow < aA.numRows(); iRow++)
        {
            lReturnMatrix(iRow, iCol) = aA(iRow, iCol) + aB(iRow, iCol);
        }
    }
    
    return lReturnMatrix;
}
void Add(NOX::LAPACK::Matrix<double>& aTarget, const NOX::LAPACK::Matrix<double>& aAddedMatrix)
{
    if (aTarget.numRows() != aAddedMatrix.numRows() || aTarget.numCols() != aAddedMatrix.numCols())
        throw "Incompatible matrix sizes!";
    
    for (int iCol = 0; iCol < aTarget.numCols(); iCol++)
    {
        for (int iRow = 0; iRow < aTarget.numRows(); iRow++)
        {
            aTarget(iRow, iCol) += aAddedMatrix(iRow, iCol);
        }
    }
}
