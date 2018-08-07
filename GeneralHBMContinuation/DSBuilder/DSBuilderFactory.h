
#pragma once
#include <functional>
#include <map>
#include <string>

#include "DSBuilderBase.h"
#include "DSBuilderSimple.h"

// just simple creation functions (should do new because the returned pointer will be deleted later)
const std::map<std::string, std::function<DSBuilderBase*()>> C_DSBuilderFactory = 
{
    { 
        std::string("Simple"),    
        []()
        { 
            return new DSBuilderSimple();
        }
    },
};

