/*
 *  unittest.cpp
 *  levelset
 *
 *  Created by David Chopp on 4/2/10.
 *  Copyright 2010 Northwestern University. All rights reserved.
 *
 */


#include <iostream>
#include "level/units.h"
#include <stdlib.h>

using namespace levelset;

int main(int argc, char* argv[])
{
    Units units;

    double val = atof(argv[1]);
    double ans = units.Convert(val, std::string(argv[2]), std::string(argv[3]));

    std::cout << argv[1] << " " << argv[2] << " = " << ans << " " << argv[3] << '\n';
    return 0;
}
