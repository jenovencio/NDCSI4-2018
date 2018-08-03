
#pragma once
#include "Nonlinearities/NonlinearBase.h"
#include "Nonlinearities/FricElem3D.h"
#include "Nonlinearities/FResult.h"

FResult ComputeFTimeDomain(const NOX::LAPACK::Vector& aX, const NOX::LAPACK::Vector& aXPrev, const double& aNormalStaticLoad);
std::vector<NOX::LAPACK::Vector> GetForcesForMotion(const std::vector<NOX::LAPACK::Vector> aX, const double& aNormalStaticLoad);
