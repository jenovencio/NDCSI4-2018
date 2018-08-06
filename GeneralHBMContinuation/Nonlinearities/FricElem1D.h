
// class FricElem3D;
// 
// #pragma once
// #include "NonlinearBaseTD.h"
// 
// class FricElem1D : public NonlinearBaseTD
// {
// private:
//     int mDofID = 0;
//     double mN0 = 0.001;
//     double mMu = 0.2;
//     double mKtx = 1;
// protected:
//     // time domain to time domain
//     virtual FResult ComputeFTimeDomain(const NOX::LAPACK::Vector& aX, const NOX::LAPACK::Vector& aXPrev) const override;
//     // time domain to time domain
//     virtual NOX::LAPACK::Matrix<double> ComputeJacobianTimeDomain(const NOX::LAPACK::Vector& aX, const NOX::LAPACK::Vector& aXPrev) const override;
//     virtual int NumberOfPrepLoops() const override;
//     virtual std::vector<int> NonzeroFPositions() const override;
//     virtual bool IsHistoryDependent() const override;
// public:
//     virtual void LoadFromFile(const std::string& aFilePath) override;
//     virtual std::string ClassName() const override;
// };
