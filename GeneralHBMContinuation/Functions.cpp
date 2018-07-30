
#include <algorithm>
#include <string>

#include "Functions.h"
#include "Misc.h"
#include "Nonlinearities/NonlinearitiesFactory.h"

std::map<std::string, std::vector<std::string>> ParseGeneralConfigFile(const std::string& aFilePath)
{
    std::ifstream lInputFile(aFilePath);
    if (!lInputFile.is_open()) throw std::string("Unable to open file \"" + aFilePath + "\"!");
    
    std::map<std::string, std::vector<std::string>> lReturnMap;
    
    bool lIsValidLine = true;
    
    std::string lLastKey = "";
    
    while (true)
    {
        std::string lLine = GetNextValidLine(lInputFile, lIsValidLine);
        if (!lIsValidLine) break;
        
        if (lLine[0] == KEY_PREFIX_CHAR)
        {
            // this line is a key
            // skip the key prefix character
            lLine = lLine.substr(1);
            
            // skip possible whitespace after the key prefix character
            lLine = SkipWhiteSpaces(lLine);
            
            if (lReturnMap.find(lLine) != lReturnMap.end()) throw "Key \"" + lLine + "\" occurs in the file \"" + aFilePath + "\" multiple times!";
                        
            lReturnMap[lLine] = std::vector<std::string>();
            
            lLastKey = lLine;
        }
        else
        {
            if (lLastKey.size() == 0) throw "First valid line in the configuration file must be a key definition!";
            
            lReturnMap[lLastKey].push_back(lLine);
        }
    }
    
    return lReturnMap;
}



NOX::LAPACK::Matrix<double> LoadSquareMatrixFull(const std::string& aFilePath, int& aDim)
{
    std::ifstream lInputFile(aFilePath);
    if (!lInputFile.is_open()) throw "Unable to open file \"" + aFilePath + "\"!";
    
    int lDim;
    
    lInputFile >> lDim;
    
    aDim = lDim;
    
    NOX::LAPACK::Matrix<double> lReturnMatrix(lDim, lDim);
    
    int lRow = 0;
    int lCol = 0;
    for (int i = 0; i < lDim * lDim; i++)
    {
        lInputFile >> lReturnMatrix(lRow, lCol);
        
        lRow++;
        if (lRow == lDim)
        {
            lCol++;
            lRow = 0;
        }
    }
    
    lInputFile.close();
    
    return lReturnMatrix;
}
NOX::LAPACK::Matrix<double> LoadSquareMatrixSparse(const std::string& aFilePath, int& aDim)
{
    std::ifstream lInputFile(aFilePath);
    if (!lInputFile.is_open()) throw "Unable to open file \"" + aFilePath + "\"!";
    
    int lDim;
    lInputFile >> lDim;
    
    aDim = lDim;
        
    NOX::LAPACK::Matrix<double> lReturnMatrix(lDim, lDim);
    
    int lRow;
    int lCol;
    double lValue;
    
    int lEntryNumber = 0;
    
    while (!lInputFile.eof())
    {
        lInputFile >> lRow;
        lInputFile >> lCol;
        
        if (lRow < 0 || lRow >= lDim) throw "Invalid row value! (entry number " + std::to_string(lEntryNumber) + ")";
        if (lCol < 0 || lCol >= lDim) throw "Invalid col value! (entry number " + std::to_string(lEntryNumber) + ")";
        
        lInputFile >> lValue;
        
        lReturnMatrix(lRow, lCol) = lValue;
        
        lEntryNumber++;
    }
    
    lInputFile.close();
    
    return lReturnMatrix;
}
NOX::LAPACK::Matrix<double> LoadSquareMatrix(const std::string& aFilePath, const std::string& aMatrixType, int& aDim)
{
    if (aMatrixType == FULL_STRING)
        return LoadSquareMatrixFull(aFilePath, aDim);
    else if (aMatrixType == SPARSE_STRING)
        return LoadSquareMatrixSparse(aFilePath, aDim);
    else throw "Unknown matrix type: " + aMatrixType + " for file \"" + aFilePath + "\"!";
}
std::vector<double> LoadExcitationForce(const std::string& aFilePath, int& aDim, int& aHarmCoeffCount)
{
    std::ifstream lInputFile(aFilePath);
    if (!lInputFile.is_open()) throw "Unable to open file \"" + aFilePath + "\"!";
    
    int lDim;
    lInputFile >> lDim;
    
    aDim = lDim;
        
    int lHarmCount;
    lInputFile >> lHarmCount;
    
    aHarmCoeffCount = lHarmCount;
    
    std::vector<double> lReturnVector(lDim * lHarmCount);
    
    for (int iDof = 0; iDof < lDim; iDof++)
    {
        for (int iHarm = 0; iHarm < lHarmCount; iHarm++)
        {
            int lIndex = GetHBMDofIndex(iDof, iHarm, lHarmCount);
            lInputFile >> lReturnVector[lIndex];
        }
    }
    
    lInputFile.close();
    
    return lReturnVector;
}

std::string GetNextValidLine(std::ifstream& aFile, bool& aIsValid)
{
    std::string lReturnValue = "";
    
    bool lFoundValidLine = false;
    
    while (!aFile.eof())
    {
        std::getline(aFile, lReturnValue);
        
        lReturnValue = SkipWhiteSpaces(lReturnValue);
        
        if (lReturnValue.size() == 0) continue;
        if (lReturnValue[0] == '#') continue;
        
        lFoundValidLine = true;
        break;
    }
    
    aIsValid = lFoundValidLine;
    
    return lReturnValue;
}

void CheckString(const std::string& aString, const std::vector<std::string>& aPossibilities, const std::string& aGroupName)
{
    bool lIsIn = std::find(aPossibilities.begin(), aPossibilities.end(), aString) != aPossibilities.end();
    
    if (!lIsIn)
    {
        std::stringstream lStringBuilder;
        lStringBuilder << "Value \"" + aString + "\" is not a valid option for \"" + aGroupName + "\"! Valid options are: " << std::endl;
        
        if (aPossibilities.size() == 0) lStringBuilder << " none (you are fucked)";
        
        for (int i = 0; i < aPossibilities.size(); i++)
        {
            lStringBuilder << aPossibilities[i];
            if (i < aPossibilities.size() - 1)
                lStringBuilder << ", ";
        }
        
        lStringBuilder << std::endl;
        
        throw lStringBuilder.str();
    }
}

std::vector<double> GetRelativeTimePoints(const int& aTimePointCount)
{
    double lRelativeStep = 1.0 / aTimePointCount;
    std::vector<double> lReturnVector;
    
    lReturnVector.reserve(aTimePointCount);
    
    double lCurrentPos = lRelativeStep / 2;
    
    for (int i = 0; i < aTimePointCount; i++)
    {
        lReturnVector.push_back(lCurrentPos);
        lCurrentPos += lRelativeStep;
    }
    
    return lReturnVector;
}
NOX::LAPACK::Vector GetAverage(const std::vector<NOX::LAPACK::Vector>& aVectors)
{
    if (aVectors.size() <= 0) throw "No vectors to average!";
    
    NOX::LAPACK::Vector lReturnVector(aVectors[0].length());
    
    for (int iVec = 0; iVec < aVectors.size(); iVec++)
    {
        for (int iElem = 0; iElem < aVectors[0].length(); iElem++)
        {
            lReturnVector(iElem) += aVectors[iVec](iElem);
        }
    }
    
    lReturnVector.scale(1.0 / aVectors.size());
    
    return lReturnVector;
}
std::string SkipWhiteSpaces(const std::string& aString)
{
    int lFirstNonWhitespace = 0;
    while (lFirstNonWhitespace < aString.size())
    {
        char lChar = aString[lFirstNonWhitespace];
        if (lChar != ' ' && lChar != '\t') break;
        lFirstNonWhitespace++;
    }
    
    if (lFirstNonWhitespace >= aString.size()) return "";
    
    std::string lReturnValue = aString.substr(lFirstNonWhitespace);
    
    return lReturnValue;
}
