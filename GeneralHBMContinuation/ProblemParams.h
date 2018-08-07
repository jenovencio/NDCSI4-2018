
struct ProblemParams;

#pragma once

struct ProblemParams
{
public:
    int DofCountPhysical;
    int HarmonicCount;
    int DofCountHBM;
    
    ProblemParams(const int& aDofPhysicalCount, const int& aHarmonicCount)
        : DofCountPhysical(aDofPhysicalCount), HarmonicCount(aHarmonicCount), DofCountHBM(aDofPhysicalCount * aHarmonicCount)
    {
    }
    
    ProblemParams(const ProblemParams& aOther)
        : ProblemParams(aOther.DofCountPhysical, aOther.HarmonicCount)
    {
    }
    ProblemParams()
        : ProblemParams(0, 0)
    {
    }
};
