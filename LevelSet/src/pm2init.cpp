//
//  pm2init.cpp
//  LevelSet
//
//  Created by David Chopp on 7/24/19.
//  Copyright © 2019 David Chopp. All rights reserved.
//

#include <math.h>
#include "defs.h"
#include "polarmesh2d.h"
#include "pm2boundary.h"

namespace levelset {
    
    void PolarMesh2D::SetValuesXY(const InitialFunc& f, const int wi, const double mult)
    {
        for (int i=0; i<maxi; ++i) {
            for (int j=0; j<maxj; ++j) {
                data_(i,j,wi) = f.XY(X(i,j),Y(i,j))*mult;
            }
        }
        bc->Apply(wi);
    }
    
    void PolarMesh2D::SetValuesRT(const InitialFunc& f, const int wi, const double mult)
    {
        for (int i=0; i<maxi; ++i) {
            double rho = r[i];
            for (int j=0; j<maxj; ++j) {
                double theta = j*dtheta;
                data_(i,j,wi) = f.XY(rho,theta)*mult;
            }
        }
        bc->Apply(wi);
    }
    
    void PolarMesh2D::SetValuesXY(double (*f)(double,double), const int wi, const double mult)
    {
        for (int i=0; i<maxi; ++i) {
            for (int j=0; j<maxj; ++j) {
                data_(i,j,wi) = (*f)(X(i,j),Y(i,j))*mult;
            }
        }
        bc->Apply(wi);
    }
    
    void PolarMesh2D::SetValuesRT(double (*f)(double,double), const int wi, const double mult)
    {
        for (int i=0; i<maxi; ++i) {
            double rho = r[i];
            for (int j=0; j<maxj; ++j) {
                double theta = j*dtheta;
                data_(i,j,wi) = (*f)(rho,theta)*mult;
            }
        }
        bc->Apply(wi);
    }
    
    void PolarMesh2D::SetValues(const double* f, const int wi, const double mult)
    {
        for (int i=0; i<maxi; ++i) {
            for (int j=0; j<maxj; ++j) {
                data_(i,j,wi) = f[i*maxj+j]*mult;
            }
        }
        bc->Apply(wi);
    }

    void PolarMesh2D::SetValues(const double c, const int wi)
    {
        for (int i=0; i<maxi; ++i) {
            for (int j=0; j<maxj; ++j) {
                data_(i,j,wi) = c;
            }
        }
        bc->Apply(wi);
    }

}
