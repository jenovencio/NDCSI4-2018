
#pragma once
#include <string>
#include <map>
#include <functional>

#include "NonlinearBase.h"
#include "CubicSpring.h"
#include "FricElem1D.h"
#include "FricElem3D.h"
#include "PenaltyWall.h"

// just simple creation functions (should do new because the returned pointer will be deleted later)
const std::map<std::string, std::function<NonlinearBase*()>> C_NonlinearitiesFactory = 
{
    { std::string("cubic spring"), []() { return new CubicSpring(); } },
    { std::string("fric elem 1D"), []() { return new FricElem1D(); } },
//     { std::string("fric elem 3D"), []() { return new FricElem3D(); } },
    { std::string("penalty wall"), []() { return new PenaltyWall(); } },
};
