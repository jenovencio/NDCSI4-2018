
#pragma once
#include <string>
#include <map>
#include <functional>

#include "NonlinearBase.h"
#include "CubicSpring.h"

// just simple creation functions (should do new because the returned pointer will be deleted later
const std::map<std::string, std::function<NonlinearBase*()>> C_NonlinearitiesFactory = 
{
    { std::string("cubic spring"), []() { return new CubicSpring(); } } 
};
