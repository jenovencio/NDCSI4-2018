
#include "ProblemInterface.h"

ProblemInterface::ProblemInterface(const Config& aConfig)
 : cDOFCount(aConfig.DOFCount)
{
}

void ProblemInterface::setParams(const LOCA::ParameterVector& aParams)
{
    mFrequency = aParams.getValue("frequency");
}

