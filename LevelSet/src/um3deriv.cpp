#include "uniformmesh3d.h"
#include "um3boundary.h"
#include "utility.h"

namespace levelset {
        
    void UniformMesh3D::NormUpwindGrad(const int v, const int kdxm, const int kdxp,
                                       const int kdym, const int kdyp,
                                       const int kdzm, const int kdzp,
                                       const int knorm)
    {
        int i, j, k;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                for (k=0; k<maxk; ++k)
                    if (data_(i,j,k,v) > 0.) 
                        data_(i,j,k,knorm) = sqrt(sqr(max(data_(i,j,k,kdxp),
                                                          -data_(i,j,k,kdxm),0.))
                                                  +sqr(max(data_(i,j,k,kdyp),
                                                           -data_(i,j,k,kdym),0.))
                                                  +sqr(max(data_(i,j,k,kdzp),
                                                           -data_(i,j,k,kdzm),0.)));
                    else
                        data_(i,j,k,knorm) = sqrt(sqr(max(-data_(i,j,k,kdxp),
                                                          data_(i,j,k,kdxm),0.))
                                                  +sqr(max(-data_(i,j,k,kdyp),
                                                           data_(i,j,k,kdym),0.))
                                                  +sqr(max(-data_(i,j,k,kdzp),
                                                           data_(i,j,k,kdzm),0.)));
        bc->Apply(knorm);
    }

    void UniformMesh3D::NormUpwindGrad(const int v, const int l, const int knorm, 
                                       const int dir)
    {
        int i, j, k;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                for (k=0; k<maxk; ++k) 
                    if (data_(i,j,k,v)*dir > 0.) 
                        data_(i,j,k,knorm) = sqrt(sqr(max((data_(i+1,j,k,l)-data_(i,j,k,l))/dx,
                                                          -(data_(i,j,k,l)-data_(i-1,j,k,l))/dx,0.))
                                                  +sqr(max((data_(i,j+1,k,l)-data_(i,j,k,l))/dy,
                                                           -(data_(i,j,k,l)-data_(i,j-1,k,l))/dy,0.))
                                                  +sqr(max((data_(i,j,k+1,l)-data_(i,j,k,l))/dz,
                                                           -(data_(i,j,k,l)-data_(i,j,k-1,l))/dz,0.)));
                    else
                        data_(i,j,k,knorm) = sqrt(sqr(max(-(data_(i+1,j,k,l)-data_(i,j,k,l))/dx,
                                                          (data_(i,j,k,l)-data_(i-1,j,k,l))/dx,0.))
                                                  +sqr(max(-(data_(i,j+1,k,l)-data_(i,j,k,l))/dy,
                                                           (data_(i,j,k,l)-data_(i,j-1,k,l))/dy,0.))
                                                  +sqr(max(-(data_(i,j,k+1,l)-data_(i,j,k,l))/dz,
                                                           (data_(i,j,k,l)-data_(i,j,k-1,l))/dz,0.)));
        bc->Apply(knorm);
    }

    void UniformMesh3D::NormGrad(const int kdx, const int kdy, const int kdz,
                                 const int knorm)
    {
        int i, j, k;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                for (k=0; k<maxk; ++k)
                    data_(i,j,k,knorm) = sqrt(sqr(data_(i,j,k,kdx))
                                              +sqr(data_(i,j,k,kdy))
                                              +sqr(data_(i,j,k,kdz)));
        bc->Apply(knorm);
    }

    double UniformMesh3D::AvgGrad(const int kdx, const int kdy, const int kdz)
    {
        double answer = 0;
        for (int i=0; i<maxi; ++i)
            for (int j=0; j<maxj; ++j)
                for (int k=0; k<maxk; ++k)
                {
                    answer += sqr(data_(i,j,k,kdx))+sqr(data_(i,j,k,kdy))
                        +sqr(data_(i,j,k,kdz));
                }
        return answer/(maxi*maxj*maxk);
    }

    double UniformMesh3D::MaxGrad(const int kdx, const int kdy, const int kdz)
    {
        double answer = 0;
        for (int i=0; i<maxi; ++i)
            for (int j=0; j<maxj; ++j)
                for (int k=0; k<maxk; ++k)
                {
                    double grad = sqr(data_(i,j,k,kdx))+sqr(data_(i,j,k,kdy))
                        +sqr(data_(i,j,k,kdz));
                    answer = max(answer, grad);
                }
        return sqrt(answer);
    }

    void UniformMesh3D::Curvature(const int kdx, const int kdy, const int kdz,
                                  const int kdxx, const int kdxy, const int kdxz,
                                  const int kdyy, const int kdyz, const int kdzz,
                                  const int kcurv) 
    {
        int i, j, k;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                for (k=0; k<maxk; ++k)
                    data_(i,j,k,kcurv) = (data_(i,j,k,kdxx)*(sqr(data_(i,j,k,kdy))
                                                             +sqr(data_(i,j,k,kdz)))
                                          +data_(i,j,k,kdyy)*(sqr(data_(i,j,k,kdx))
                                                              +sqr(data_(i,j,k,kdz)))
                                          +data_(i,j,k,kdzz)*(sqr(data_(i,j,k,kdx))
                                                              +sqr(data_(i,j,k,kdy)))
                                          -2.*data_(i,j,k,kdxy)*data_(i,j,k,kdx)
                                          *data_(i,j,k,kdy)
                                          -2.*data_(i,j,k,kdxz)*data_(i,j,k,kdx)
                                          *data_(i,j,k,kdz)
                                          -2.*data_(i,j,k,kdyz)*data_(i,j,k,kdy)
                                          *data_(i,j,k,kdz))
                        /pow(sqr(data_(i,j,k,kdx))+sqr(data_(i,j,k,kdy))
                             +sqr(data_(i,j,k,kdz)),1.5);
        bc->Apply(kcurv);
    }

    void UniformMesh3D::Weight(double w(const double theta, const double phi),
                               const int kdx, const int kdy, const int kdz,
                               const int kcurv) 
    {
        int i, j, k;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j) 
                for (k=0; k<maxk; ++k)
                    data_(i,j,k,kcurv) *= w(atan2(data_(i,j,k,kdy),data_(i,j,k,kdx)),
                                            atan2(data_(i,j,k,kdz),
                                                  sqrt(sqr(data_(i,j,k,kdx))+sqr(data_(i,j,k,kdy)))));
        bc->Apply(kcurv);
    }


    void UniformMesh3D::Ds_zero(const int kf, const int kdx, const int kdy,
                                const int kdz, const int kds, const int kdt,
                                const int kdu) 
    {
        int i, j, k;
        double norm, f_x, f_y, f_z;
   
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                for (k=0; k<maxk; ++k) {
                    norm = sqr(data_(i,j,k,kdx))+sqr(data_(i,j,k,kdy))
                        +sqr(data_(i,j,k,kdz));
                    f_x = (data_(i+1,j,k,kf)-data_(i-1,j,k,kf))/2./dx;
                    f_y = (data_(i,j+1,k,kf)-data_(i,j-1,k,kf))/2./dy;
                    f_z = (data_(i,j,k+1,kf)-data_(i,j,k-1,kf))/2./dz;
                    data_(i,j,k,kds) = (f_x*(sqr(data_(i,j,k,kdy))
                                             +sqr(data_(i,j,k,kdz)))
                                        -f_y*data_(i,j,k,kdx)*data_(i,j,k,kdy)
                                        -f_z*data_(i,j,k,kdx)*data_(i,j,k,kdz))/norm;
                    data_(i,j,k,kdt) = (f_y*(sqr(data_(i,j,k,kdx))
                                             +sqr(data_(i,j,k,kdz)))
                                        -f_x*data_(i,j,k,kdx)*data_(i,j,k,kdy)
                                        -f_z*data_(i,j,k,kdy)*data_(i,j,k,kdz))/norm;
                    data_(i,j,k,kdu) = (f_z*(sqr(data_(i,j,k,kdy))
                                             +sqr(data_(i,j,k,kdx)))
                                        -f_y*data_(i,j,k,kdz)*data_(i,j,k,kdy)
                                        -f_x*data_(i,j,k,kdx)*data_(i,j,k,kdz))/norm;
                }

        bc->Apply(kds);
    }

    void UniformMesh3D::Dss_zero(const int kf, const int kdx, const int kdy,
                                 const int kdz, const int kdss)
    {
        double fx, fy, fz, fxx, fyy, fzz, fxy, fxz, fyz, grad;
        int i, j, k;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                for (k=0; k<maxk; ++k) {
                    fx = (-data_(i+2,j,k,kf)+8*data_(i+1,j,k,kf)
                          -8*data_(i-1,j,k,kf)+data_(i-2,j,k,kf))/12./dx;
                    fy = (-data_(i,j+2,k,kf)+8*data_(i,j+1,k,kf)
                          -8*data_(i,j-1,k,kf)+data_(i,j-2,k,kf))/12./dy;
                    fz = (-data_(i,j,k+2,kf)+8*data_(i,j,k+1,kf)
                          -8*data_(i,j,k-1,kf)+data_(i,j,k-2,kf))/12./dz;
                    fxx = (-data_(i+2,j,k,kf)+16*data_(i+1,j,k,kf)
                           -30*data_(i,j,k,kf)+16*data_(i-1,j,k,kf)
                           -data_(i-2,j,k,kf))/12./dx/dx;
                    fyy = (-data_(i,j+2,k,kf)+16*data_(i,j+1,k,kf)
                           -30*data_(i,j,k,kf) +16*data_(i,j-1,k,kf)
                           -data_(i,j-2,k,kf))/12./dy/dy;
                    fzz = (-data_(i,j,k+2,kf)+16*data_(i,j,k+1,kf)
                           -30*data_(i,j,k,kf) +16*data_(i,j,k-1,kf)
                           -data_(i,j,k-2,kf))/12./dz/dz;
                    fxy = (-data_(i-2,j-2,k,kf)+16*data_(i-1,j-1,k,kf)
                           +data_(i+2,j-2,k,kf)-16*data_(i+1,j-1,k,kf)
                           +data_(i-2,j+2,k,kf)-16*data_(i-1,j+1,k,kf)
                           -data_(i+2,j+2,k,kf)+16*data_(i+1,j+1,k,kf))/48/dx/dy;
                    fxz = (-data_(i-2,j,k-2,kf)+16*data_(i-1,j,k-1,kf)
                           +data_(i+2,j,k-2,kf)-16*data_(i+1,j,k-1,kf)
                           +data_(i-2,j,k+2,kf)-16*data_(i-1,j,k+1,kf)
                           -data_(i+2,j,k+2,kf)+16*data_(i+1,j,k+1,kf))/48/dx/dz;
                    fyz = (-data_(i,j-2,k-2,kf)+16*data_(i,j-1,k-1,kf)
                           +data_(i,j+2,k-2,kf)-16*data_(i,j+1,k-1,kf)
                           +data_(i,j-2,k+2,kf)-16*data_(i,j-1,k+1,kf)
                           -data_(i,j+2,k+2,kf)+16*data_(i,j+1,k+1,kf))/48/dy/dz;
                    grad = sqrt(sqr(data_(i,j,k,kdx))+sqr(data_(i,j,k,kdy))
                                +sqr(data_(i,j,k,kdz)));
                    data_(i,j,k,kdss) = fxx*(1.-sqr(data_(i,j,k,kdx)/grad))
                        +fyy*(1.-sqr(data_(i,j,k,kdy)/grad))
                        +fzz*(1.-sqr(data_(i,j,k,kdz)/grad))
                        -(data_(i,j,k,kdx)*fx+data_(i,j,k,kdy)*fy
                          +data_(i,j,k,kdz)*fz)*data_(i,j,k,kf)/grad;
                }
#ifdef LEVEL_DEBUG
        TestValidity(__FILE__,__LINE__,kdss);
#endif
        bc->Apply(kdss);
    }

    void UniformMesh3D::InterpGrad(const double x, const double y, const double z, 
                                   const int l, double& gx, double& gy, double& gz) const
    {
        int i, j, k;
        double ifrac, jfrac, kfrac;
        LocToIndex(x, y, z, i, j, k, ifrac, jfrac, kfrac);
        double igrad[2][2][2][3];
        for (int ii=0; ii<=1; ++ii)
            for (int jj=0; jj<=1; ++jj)
                for (int kk=0; kk<=1; ++kk) {
                    igrad[ii][jj][kk][0] = (data_(i+ii+1,j+jj,k+kk,l)
                                            -data_(i+ii-1,j+jj,k+kk,l))/2/dx;
                    igrad[ii][jj][kk][1] = (data_(i+ii,j+jj+1,k+kk,l)
                                            -data_(i+ii,j+jj-1,k+kk,l))/2/dy;
                    igrad[ii][jj][kk][2] = (data_(i+ii,j+jj,k+kk+1,l)
                                            -data_(i+ii,j+jj,k+kk-1,l))/2/dz;
                }
        gx = (1.-ifrac)*(1.-jfrac)*(1.-kfrac)*igrad[0][0][0][0]
            +ifrac*(1.-jfrac)*(1.-kfrac)*igrad[1][0][0][0]
            +(1.-ifrac)*jfrac*(1.-kfrac)*igrad[0][1][0][0]
            +ifrac*jfrac*(1.-kfrac)*igrad[1][1][0][0]
            +(1.-ifrac)*(1.-jfrac)*kfrac*igrad[0][0][1][0]
            +ifrac*(1.-jfrac)*kfrac*igrad[1][0][1][0]
            +(1.-ifrac)*jfrac*kfrac*igrad[0][1][1][0]
            +ifrac*jfrac*kfrac*igrad[1][1][1][0];
        gy = (1.-ifrac)*(1.-jfrac)*(1.-kfrac)*igrad[0][0][0][1]
            +ifrac*(1.-jfrac)*(1.-kfrac)*igrad[1][0][0][1]
            +(1.-ifrac)*jfrac*(1.-kfrac)*igrad[0][1][0][1]
            +ifrac*jfrac*(1.-kfrac)*igrad[1][1][0][1]
            +(1.-ifrac)*(1.-jfrac)*kfrac*igrad[0][0][1][1]
            +ifrac*(1.-jfrac)*kfrac*igrad[1][0][1][1]
            +(1.-ifrac)*jfrac*kfrac*igrad[0][1][1][1]
            +ifrac*jfrac*kfrac*igrad[1][1][1][1];
        gz = (1.-ifrac)*(1.-jfrac)*(1.-kfrac)*igrad[0][0][0][2]
            +ifrac*(1.-jfrac)*(1.-kfrac)*igrad[1][0][0][2]
            +(1.-ifrac)*jfrac*(1.-kfrac)*igrad[0][1][0][2]
            +ifrac*jfrac*(1.-kfrac)*igrad[1][1][0][2]
            +(1.-ifrac)*(1.-jfrac)*kfrac*igrad[0][0][1][2]
            +ifrac*(1.-jfrac)*kfrac*igrad[1][0][1][2]
            +(1.-ifrac)*jfrac*kfrac*igrad[0][1][1][2]
            +ifrac*jfrac*kfrac*igrad[1][1][1][2];
    }


    void UniformMesh3D::Dx_zero(const int l, const int kdx)
    {
        int i, j, k;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                for (k=0; k<maxk; ++k)
                    data_(i,j,k,kdx) = (data_(i+1,j,k,l)-data_(i-1,j,k,l))/2./dx;
        bc->Apply(kdx);
    }

    void UniformMesh3D::Dx_minus(const int l, const int kdx)
    {
        int i, j, k;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                for (k=0; k<maxk; ++k)
                    data_(i,j,k,kdx) = (data_(i,j,k,l)-data_(i-1,j,k,l))/dx;
        bc->Apply(kdx);
    }

    void UniformMesh3D::Dx_plus(const int l, const int kdx)
    {
        int i, j, k;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                for (k=0; k<maxk; ++k)
                    data_(i,j,k,kdx) = (data_(i+1,j,k,l)-data_(i,j,k,l))/dx;
        bc->Apply(kdx);
    }

    void UniformMesh3D::Dx_upwind(const int v, const int l, const int kdx)
    {
        double dxp, dxm;
        int i, j, k;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                for (k=0; k<maxk; ++k)
                {
                    dxp = (data_(i+1,j,k,l)-data_(i,j,k,l))/dx;
                    dxm = (data_(i,j,k,l)-data_(i-1,j,k,l))/dx;
                    if (data_(i,j,k,v) >= 0.) 
                        data_(i,j,k,kdx) = min(dxm,0.)+max(dxp,0.);
                    else
                        data_(i,j,k,kdx) = max(dxm,0.)+min(dxp,0.);
                }
        bc->Apply(kdx);
    }
   

    void UniformMesh3D::Dxx_zero(const int l, const int kdx)
    {
        int i, j, k;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                for (k=0; k<maxk; ++k)
                    data_(i,j,k,kdx) = (data_(i+1,j,k,l)-2.*data_(i,j,k,l)
                                        +data_(i-1,j,k,l))/dx/dx;
        bc->Apply(kdx);
    }

    void UniformMesh3D::Dy_zero(const int l, const int kdy)
    {
        int i, j, k;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                for (k=0; k<maxk; ++k)
                    data_(i,j,k,kdy) = (data_(i,j+1,k,l)-data_(i,j-1,k,l))/2./dy;
        bc->Apply(kdy);
    }

    void UniformMesh3D::Dy_minus(const int l, const int kdy)
    {
        int i, j, k;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                for (k=0; k<maxk; ++k)
                    data_(i,j,k,kdy) = (data_(i,j,k,l)-data_(i,j-1,k,l))/dy;
        bc->Apply(kdy);
    }

    void UniformMesh3D::Dy_plus(const int l, const int kdy)
    {
        int i, j, k;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                for (k=0; k<maxk; ++k)
                    data_(i,j,k,kdy) = (data_(i,j+1,k,l)-data_(i,j,k,l))/dy;
        bc->Apply(kdy);
    }

    void UniformMesh3D::Dy_upwind(const int v, const int l, const int kdy)
    {
        double dyp, dym;
        int i, j, k;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                for (k=0; k<maxk; ++k) {
                    dyp = (data_(i,j+1,k,l)-data_(i,j,k,l))/dy;
                    dym = (data_(i,j,k,l)-data_(i,j-1,k,l))/dy;
                    if (data_(i,j,k,v) >= 0.) 
                        data_(i,j,k,kdy) = min(dym,0.)+max(dyp,0.);
                    else
                        data_(i,j,k,kdy) = max(dym,0.)+min(dyp,0.);
                }
        bc->Apply(kdy);
    }

    void UniformMesh3D::Dyy_zero(const int l, const int kdy)
    {
        int i, j, k;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                for (k=0; k<maxk; ++k)
                    data_(i,j,k,kdy) = (data_(i,j+1,k,l)-2.*data_(i,j,k,l)
                                        +data_(i,j-1,k,l))/dy/dy;
        bc->Apply(kdy);
    }

    void UniformMesh3D::Dz_zero(const int l, const int kdz)
    {
        int i, j, k;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                for (k=0; k<maxk; ++k)
                    data_(i,j,k,kdz) = (data_(i,j,k+1,l)-data_(i,j,k-1,l))/2./dz;
        bc->Apply(kdz);
    }

    void UniformMesh3D::Dz_minus(const int l, const int kdz)
    {
        int i, j, k;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                for (k=0; k<maxk; ++k)
                    data_(i,j,k,kdz) = (data_(i,j,k,l)-data_(i,j,k-1,l))/dz;
        bc->Apply(kdz);
    }

    void UniformMesh3D::Dz_plus(const int l, const int kdz)
    {
        int i, j, k;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                for (k=0; k<maxk; ++k)
                    data_(i,j,k,kdz) = (data_(i,j,k+1,l)-data_(i,j,k,l))/dz;
        bc->Apply(kdz);
    }

    void UniformMesh3D::Dz_upwind(const int v, const int l, const int kdz)
    {
        double dzp, dzm;
        int i, j, k;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                for (k=0; k<maxk; ++k) {
                    dzp = (data_(i,j,k+1,l)-data_(i,j,k,l))/dz;
                    dzm = (data_(i,j,k,l)-data_(i,j,k-1,l))/dz;
                    if (data_(i,j,k,v) >= 0.) 
                        data_(i,j,k,kdz) = min(dzm,0.)+max(dzp,0.);
                    else
                        data_(i,j,k,kdz) = max(dzm,0.)+min(dzp,0.);
                }
        bc->Apply(kdz);
    }

    void UniformMesh3D::Dzz_zero(const int l, const int kdz)
    {
        int i, j, k;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                for (k=0; k<maxk; ++k)
                    data_(i,j,k,kdz) = (data_(i,j,k+1,l)-2.*data_(i,j,k,l)
                                        +data_(i,j,k-1,l))/dz/dz;
        bc->Apply(kdz);
    }

    void UniformMesh3D::Dxy_zero(const int l, const int kdxy)
    {
        int i, j, k;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                for (k=0; k<maxk; ++k)
                    data_(i,j,k,kdxy) = (data_(i-1,j-1,k,l)-data_(i+1,j-1,k,l)
                                         -data_(i-1,j+1,k,l)+data_(i+1,j+1,k,l)
                        )/4./dx/dy;
        bc->Apply(kdxy);
    }

    void UniformMesh3D::Dxz_zero(const int l, const int kdxz)
    {
        int i, j, k;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                for (k=0; k<maxk; ++k)
                    data_(i,j,k,kdxz) = (data_(i-1,j,k-1,l)-data_(i+1,j,k-1,l)
                                         -data_(i-1,j,k+1,l)+data_(i+1,j,k+1,l)
                        )/4./dx/dz;
        bc->Apply(kdxz);
    }

    void UniformMesh3D::Dyz_zero(const int l, const int kdyz)
    {
        int i, j, k;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                for (k=0; k<maxk; ++k)
                    data_(i,j,k,kdyz) = (data_(i,j-1,k-1,l)-data_(i,j+1,k-1,l)
                                         -data_(i,j-1,k+1,l)+data_(i,j+1,k+1,l)
                        )/4./dy/dz;
        bc->Apply(kdyz);
    }

    void UniformMesh3D::Dx_zero_4(const int l, const int kdx)
    {
        int i, j, k;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                for (k=0; k<maxk; ++k)
                    data_(i,j,k,kdx) = (-data_(i+2,j,k,l)+8*data_(i+1,j,k,l)
                                        -8*data_(i-1,j,k,l)+data_(i-2,j,k,l))/12./dx;
        bc->Apply(kdx);
    }

    void UniformMesh3D::Dx_minus_4(const int l, const int kdx)
    {
        int i, j, k;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                for (k=0; k<maxk; ++k)
                    data_(i,j,k,kdx) = (3*data_(i,j,k,l)-4*data_(i-1,j,k,l)
                                        +data_(i-2,j,k,l))/2./dx;
        bc->Apply(kdx);
    }

    void UniformMesh3D::Dx_plus_4(const int l, const int kdx)
    {
        int i, j, k;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                for (k=0; k<maxk; ++k)
                    data_(i,j,k,kdx) = (-data_(i+2,j,k,l)+4*data_(i+1,j,k,l)
                                        -3*data_(i,j,k,l))/2./dx;
        bc->Apply(kdx);
    }


    void UniformMesh3D::Dxx_zero_4(const int l, const int kdx)
    {
        int i, j, k;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                for (k=0; k<maxk; ++k)
                    data_(i,j,k,kdx) = (-data_(i+2,j,k,l)+16*data_(i+1,j,k,l)
                                        -30*data_(i,j,k,l)+16*data_(i-1,j,k,l)
                                        -data_(i-2,j,k,l))/12./dx/dx;
        bc->Apply(kdx);
    }

    void UniformMesh3D::Dy_zero_4(const int l, const int kdy)
    {
        int i, j, k;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                for (k=0; k<maxk; ++k)
                    data_(i,j,k,kdy) = (-data_(i,j+2,k,l)+8*data_(i,j+1,k,l)
                                        -8*data_(i,j-1,k,l)+data_(i,j-2,k,l))/12./dy;
        bc->Apply(kdy);
    }

    void UniformMesh3D::Dy_minus_4(const int l, const int kdy)
    {
        int i, j, k;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                for (k=0; k<maxk; ++k)
                    data_(i,j,k,kdy) = (3*data_(i,j,k,l)-4*data_(i,j-1,k,l)
                                        +data_(i,j-2,k,l))/2./dy;
        bc->Apply(kdy);
    }

    void UniformMesh3D::Dy_plus_4(const int l, const int kdy)
    {
        int i, j, k;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                for (k=0; k<maxk; ++k)
                    data_(i,j,k,kdy) = (-data_(i,j+2,k,l)+4*data_(i,j+1,k,l)
                                        -3*data_(i,j,k,l))/2./dy;
        bc->Apply(kdy);
    }

    void UniformMesh3D::Dyy_zero_4(const int l, const int kdy)
    {
        int i, j, k;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                for (k=0; k<maxk; ++k)
                    data_(i,j,k,kdy) = (-data_(i,j+2,k,l)+16*data_(i,j+1,k,l)
                                        -30*data_(i,j,k,l)+16*data_(i,j-1,k,l)
                                        -data_(i,j-2,k,l))/12./dy/dy;
        bc->Apply(kdy);
    }

    void UniformMesh3D::Dz_zero_4(const int l, const int kdz)
    {
        int i, j, k;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                for (k=0; k<maxk; ++k)
                    data_(i,j,k,kdz) = (-data_(i,j,k+2,l)+8*data_(i,j,k+1,l)
                                        -8*data_(i,j,k-1,l)+data_(i,j,k-2,l))/12./dz;
        bc->Apply(kdz);
    }

    void UniformMesh3D::Dz_minus_4(const int l, const int kdz)
    {
        int i, j, k;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                for (k=0; k<maxk; ++k)
                    data_(i,j,k,kdz) = (3*data_(i,j,k,l)-4*data_(i,j,k-1,l)
                                        +data_(i,j,k-2,l))/2./dz;
        bc->Apply(kdz);
    }

    void UniformMesh3D::Dz_plus_4(const int l, const int kdz)
    {
        int i, j, k;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                for (k=0; k<maxk; ++k)
                    data_(i,j,k,kdz) = (-data_(i,j,k+2,l)+4*data_(i,j,k+1,l)
                                        -3*data_(i,j,k,l))/2./dz;
        bc->Apply(kdz);
    }

    void UniformMesh3D::Dzz_zero_4(const int l, const int kdz)
    {
        int i, j, k;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                for (k=0; k<maxk; ++k)
                    data_(i,j,k,kdz) = (-data_(i,j,k+2,l)+16*data_(i,j,k+1,l)
                                        -30*data_(i,j,k,l)+16*data_(i,j,k-1,l)
                                        -data_(i,j,k-2,l))/12./dz/dz;
        bc->Apply(kdz);
    }

    void UniformMesh3D::Dxy_zero_4(const int l, const int kdxy)
    {
        int i, j, k;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                for (k=0; k<maxk; ++k)
                    data_(i,j,k,kdxy) = (-data_(i-2,j-2,k,l)+16*data_(i-1,j-1,k,l)
                                         +data_(i+2,j-2,k,l)-16*data_(i+1,j-1,k,l)
                                         +data_(i-2,j+2,k,l)-16*data_(i-1,j+1,k,l)
                                         -data_(i+2,j+2,k,l)+16*data_(i+1,j+1,k,l))/48/dx/dy;
        bc->Apply(kdxy);
    }

    void UniformMesh3D::Dxz_zero_4(const int l, const int kdxz)
    {
        int i, j, k;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                for (k=0; k<maxk; ++k)
                    data_(i,j,k,kdxz) = (-data_(i-2,j,k-2,l)+16*data_(i-1,j,k-1,l)
                                         +data_(i+2,j,k-2,l)-16*data_(i+1,j,k-1,l)
                                         +data_(i-2,j,k+2,l)-16*data_(i-1,j,k+1,l)
                                         -data_(i+2,j,k+2,l)+16*data_(i+1,j,k+1,l))/48/dx/dz;
        bc->Apply(kdxz);
    }

    void UniformMesh3D::Dyz_zero_4(const int l, const int kdyz)
    {
        int i, j, k;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                for (k=0; k<maxk; ++k)
                    data_(i,j,k,kdyz) = (-data_(i,j-2,k-2,l)+16*data_(i,j-1,k-1,l)
                                         +data_(i,j+2,k-2,l)-16*data_(i,j+1,k-1,l)
                                         +data_(i,j-2,k+2,l)-16*data_(i,j-1,k+1,l)
                                         -data_(i,j+2,k+2,l)+16*data_(i,j+1,k+1,l))/48/dy/dz;
        bc->Apply(kdyz);
    }

}


