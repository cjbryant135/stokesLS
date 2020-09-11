#include "miscio.h"
#include <sstream>
#include <iomanip>
#include <iostream>
#include <string.h>

namespace levelset {
        
// file name creator
    const char* FileName(const char* base, const int num, const char* suffix, char* result)
    {
        result[0] = '\0';
        std::ostringstream s(result);
        s << base << '.' << std::setw(6) << std::setfill('0') << num << suffix << std::ends;
        strcpy(result,s.str().c_str());
        return result;
    }

    std::string FileName(std::string base, const int num, std::string suffix)
    {
        std::ostringstream s;
        s << base << '.' << std::setw(6) << std::setfill('0') << num << suffix << std::ends;
        return s.str();
    }

    const char* dtoa(const double a, char* result)
    {
        std::ostringstream s(result);
        s << a << std::ends;
        return s.str().c_str();
    }

}
