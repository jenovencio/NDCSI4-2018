
#pragma once
#include <string>
#include <map>
#include <functional>

#include "Nonlinearities/NonlinearBase.h"
#include "Nonlinearities/CubicSpring.h"

#define PI 3.14159265358979323846264338327950288419716939937510582097494459230781
#define BORDER "-------------------------------------------------------------------"
#define OUT_EXTENSION (std::string(".out"))
#define STOP { getchar(); }

#define FULL_STRING (std::string("full"))
#define SPARSE_STRING (std::string("sparse"))

#define KEY_PREFIX_CHAR ('>')

const std::vector<std::string> C_MatrixTypes = { FULL_STRING, SPARSE_STRING };
const std::vector<std::string> C_PredictorTypes = { "Secant", "Tangent", "Constant" };

template <class T> 
inline int Sgn(T aVal)
{
    return (T(0) < aVal) - (aVal < T(0));
}
