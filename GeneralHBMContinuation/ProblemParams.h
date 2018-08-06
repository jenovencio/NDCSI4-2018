
struct ProblemParams;

#pragma once

struct ProblemParams
{
public:
    const int DofCountPhysical;
    const int HarmonicCount;
    const int DofCountHBM;
    
    ProblemParams(const int& aDofPhysicalCount, const int& aHarmonicCount)
        : DofCountPhysical(aDofPhysicalCount), HarmonicCount(aHarmonicCount), DofCountHBM(aDofPhysicalCount * aHarmonicCount)
    {
    }
    
    ProblemParams(const ProblemParams& aOther)
        : ProblemParams(aOther.DofCountPhysical, aOther.HarmonicCount)
    {
    }
};
