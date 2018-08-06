

// for templated functions

#pragma once
#include <string>
#include <vector>
#include <map>

template <class T>
void CheckString(const std::string& aString, const std::map<std::string, T>& aPossibilities, const std::string& aGroupName)
{
    bool lIsIn = aPossibilities.find(aString) != aPossibilities.end();
    
    if (!lIsIn)
    {
        std::stringstream lStringBuilder;
        lStringBuilder << "Value \"" + aString + "\" is not a valid option for \"" + aGroupName + "\"! Valid options are: " << std::endl;
        
        if (aPossibilities.size() == 0) lStringBuilder << " none (you are fucked)";
        
        int lCount = 0;
        for (auto nIt = aPossibilities.begin(); nIt != aPossibilities.end(); nIt++)
        {
            lStringBuilder << nIt->first;
            if (lCount < aPossibilities.size() - 1)
                lStringBuilder << ", ";
            
            lCount++;
        }
        
        lStringBuilder << std::endl;
        
        throw lStringBuilder.str();
    }
}

