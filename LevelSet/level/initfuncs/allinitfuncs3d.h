#ifndef __ALL3DINITFUNCS_H__
#define __ALL3DINITFUNCS_H__

#include "initfuncs/mballs.h"
#include "initfuncs/mcubes.h"
#include "initfuncs/mound.h"
#include "initfuncs/ovoid.h"
#include "initfuncs/sheet.h"
#include "initfuncs/snake.h"
#include "initfuncs/wobble.h"
#include "initfuncs/mcylinders.h"
#include "initfuncs/negative3d.h"

#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include "inputparams.h"

namespace levelset {
        
    static InitialFunc3D* PickInitialFunc(const char* name, InputParams* p)
    {
        char lower[256];
        for (unsigned int i=0; i<strlen(name); ++i)
            lower[i] = tolower(name[i]);
        lower[strlen(name)] = '\0';
        if (!strcmp(lower,"mballs")) return new Mballs(p);
        if (!strcmp(lower,"mcubes")) return new Mcubes(p);
        if (!strcmp(lower,"mound")) return new Mound(p);
        if (!strcmp(lower,"ovoid")) return new Ovoid(p);
        if (!strcmp(lower,"sheet")) return new Sheet(p);
        if (!strcmp(lower,"snake")) return new Snake(p);
        if (!strcmp(lower,"wobble")) return new Wobble(p);
        if (!strcmp(lower,"cylinders")) return new Mcylinders(p);
        std::cerr << "Error in PickInitialFunc: no such function: " << name << "\n";
        exit(1);
    }
        
}

#endif
