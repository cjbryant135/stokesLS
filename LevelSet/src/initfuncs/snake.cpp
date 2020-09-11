/*************************************************
    snake.cpp

    $Header: snake.cpp,v 2.1 99/01/06 13:59:51 chopp Exp $

    $Log:       snake.cpp,v $
Revision 2.1  99/01/06  13:59:51  13:59:51  chopp (David Chopp)
*** none ***

Revision 1.1  98/03/04  14:59:49  14:59:49  chopp (David Chopp)
Initial revision

Revision 1.1  98/03/02  12:58:14  12:58:14  chopp (David Chopp)
Initial revision


*************************************************/

#include "snake.h"
#include "utility.h"
#include <float.h>
#include <sstream>
#include <iostream>

namespace levelset {
        
    void Snake::SetParams(InputParams *params)
    {
        if (parameter) delete[] parameter;
        parameter = new double[3];
        parameter[Length] = params->GetDoubleParam("Snake Length", 1.);
        parameter[Width] = params->GetDoubleParam("Snake Width", 0.1);
        char dir[10];
        params->GetCharParam("Snake Direction", dir, "x");
        switch (dir[0]) {
        case 'X':
        case 'x':
            parameter[Direction] = 0.;
            break;
        case 'Y':
        case 'y':
            parameter[Direction] = 1.;
            break;
        case 'Z':
        case 'z':
            parameter[Direction] = 2.;
            break;
        default:
            std::cout << "Using default direction X\n";
            parameter[Direction] = 0.;
        }
        double ds = parameter[Length]/(Num-1);
        for (int i=0; i<3; ++i)
            snake[i] = new double[Num];
        for (int i=0; i<Num; ++i) {
            snake[0][i] = XPath(i*ds);
            snake[1][i] = YPath(i*ds);
            snake[2][i] = ZPath(i*ds);
        }
    }

    Snake::~Snake(void)
    {
        delete[] snake[0];
        delete[] snake[1];
        delete[] snake[2];
    }

    double Snake::XYZ(const double s, const double t, const double u) const
    {
        double answer = DBL_MAX;
        if (parameter[Direction] == 0.) {
            for (int i=0; i<Num; ++i)
                answer = min(answer,sqrt(sqr(snake[0][i]-s)+sqr(snake[1][i]-t)
                                         +sqr(snake[2][i]-u)));
        } else if (parameter[Direction] == 1.) {
            for (int i=0; i<Num; ++i)
                answer = min(answer,sqrt(sqr(snake[0][i]-t)+sqr(snake[1][i]-s)
                                         +sqr(snake[2][i]-u)));
        } else {
            for (int i=0; i<Num; ++i)
                answer = min(answer,sqrt(sqr(snake[0][i]-u)+sqr(snake[1][i]-t)
                                         +sqr(snake[2][i]-s)));
        }
   
        return parameter[Width]-answer;
    }

    double Snake::XPath(const double s)
    {
        return -1.+2.*s;
    }

    double Snake::YPath(const double s)
    {
        return sin(M_PI*s)*cos(4*M_PI*s);
    }

    double Snake::ZPath(const double s)
    {
        return sin(M_PI*s)*sin(4*M_PI*s);
    }

}






