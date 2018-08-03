
class PenaltyWall;
struct PenaltyWallDef;

#pragma once
#include "NonlinearBaseTD.h"

class PenaltyWall : public NonlinearBaseTD
{
private:
    std::vector<PenaltyWallDef> mWalls;
protected:
    
    // time domain to time domain
    virtual NOX::LAPACK::Vector ComputeFTDInner(const NOX::LAPACK::Vector& aX) const override;
    // time domain to time domain
    virtual NOX::LAPACK::Matrix<double> ComputeJacobianTimeDomain(const NOX::LAPACK::Vector& aX) const override;
//     virtual std::vector<int> NonzeroFPositions() const override;
public:
    virtual void LoadFromFile(const std::string& aFilePath) override;
    virtual std::string ClassName() const override;
    void AddWall(const PenaltyWallDef& aWall);
    void AddWall(const int& aDofIndex, const double& aDirStiffness, const double& aOffset);
};

struct PenaltyWallDef
{
public:
    int DofIndex;
    // determines the stiffness by it's magnitude and the orientation of the wall by it's sign
    // the sign is the sign of the inner normal of the wall
    double DirectionalStiffness;
    // offset of the walls (with respect to the 0 displacement)
    double Offset;
};
