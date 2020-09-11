#include "uniformmesh2d.h"
#include "um2boundary.h"
#include "utility.h"

namespace levelset {
        
    void UniformMesh2D::UnitNormal(Interface& segs, const int nx, const int ny, 
                                   const int kdx, const int kdy)
    {
        int len = segs.Length();
        Interp(segs, nx, kdx);
        Interp(segs, ny, kdy);
        double norm;
        int i;
        for (i=0; i<len; ++i) {
            norm = sqrt(sqr(segs[i][nx])+sqr(segs[i][ny]));
            segs[i][nx] /= norm;
            segs[i][ny] /= norm;
        }
    }

    void UniformMesh2D::NormUpwindGrad(Interface& segs, const int ng, 
                                       const int knorm)
    {
        Interp(segs, ng, knorm);
    }

    void UniformMesh2D::NormUpwindGrad(const int v, const int kdxm, const int kdxp,
                                       const int kdym, const int kdyp,
                                       const int knorm)
    {
        int i, j;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                if (data_(i,j,v) > 0.) 
                    data_(i,j,knorm) = sqrt(sqr(max(data_(i,j,kdxp),
                                                    -data_(i,j,kdxm),0.))
                                            +sqr(max(data_(i,j,kdyp),
                                                     -data_(i,j,kdym),0.)));
                else
                    data_(i,j,knorm) = sqrt(sqr(max(-data_(i,j,kdxp),
                                                    data_(i,j,kdxm),0.))
                                            +sqr(max(-data_(i,j,kdyp),
                                                     data_(i,j,kdym),0.)));
        bc->Apply(knorm);
    }

    void UniformMesh2D::NormUpwindGrad(const int v, const int k, const int knorm, const int dir)
    {
        int i, j;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                if (data_(i,j,v)*dir > 0.) 
                    data_(i,j,knorm) = sqrt(sqr(max((data_(i+1,j,k)-data_(i,j,k))/dx,
                                                    -(data_(i,j,k)-data_(i-1,j,k))/dx,0.))
                                            +sqr(max((data_(i,j+1,k)-data_(i,j,k))/dy,
                                                     -(data_(i,j,k)-data_(i,j-1,k))/dy,0.)));
                else
                    data_(i,j,knorm) = sqrt(sqr(max(-(data_(i+1,j,k)-data_(i,j,k))/dx,
                                                    (data_(i,j,k)-data_(i-1,j,k))/dx,0.))
                                            +sqr(max(-(data_(i,j+1,k)-data_(i,j,k))/dy,
                                                     (data_(i,j,k)-data_(i,j-1,k))/dy,0.)));
        bc->Apply(knorm);
    }

    void UniformMesh2D::NormGrad(const int kdx, const int kdy, const int knorm)
    {
        int i, j;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                data_(i,j,knorm) = sqrt(sqr(data_(i,j,kdx))
                                        +sqr(data_(i,j,kdy)));
        bc->Apply(knorm);
    }

    double UniformMesh2D::AvgGrad(const int kdx, const int kdy)
    {
        double answer = 0;
        for (int i=0; i<maxi; ++i)
            for (int j=0; j<maxj; ++j) {
                answer += sqr(data_(i,j,kdx))+sqr(data_(i,j,kdy));
            }
        return answer/(maxi*maxj);
    }

    double UniformMesh2D::MaxGrad(const int kdx, const int kdy)
    {
        double answer = 0;
        for (int i=0; i<maxi; ++i)
            for (int j=0; j<maxj; ++j) {
                double grad = sqr(data_(i,j,kdx))+sqr(data_(i,j,kdy));
                answer = max(answer, grad);
            }
        return sqrt(answer);
    }

    void UniformMesh2D::Curvature(const int kdx, const int kdy, const int kdxx,
                                  const int kdxy, const int kdyy, const int kcurv) 
    {
        int i, j;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                data_(i,j,kcurv) = (data_(i,j,kdxx)*sqr(data_(i,j,kdy))
                                    +data_(i,j,kdyy)*sqr(data_(i,j,kdx))
                                    -2.*data_(i,j,kdxy)*data_(i,j,kdx)
                                    *data_(i,j,kdy))
                    /(pow(sqr(data_(i,j,kdx))+sqr(data_(i,j,kdy)),1.5)+DBL_EPSILON);
        bc->Apply(kcurv);
    }

    void UniformMesh2D::Weight(double w(const double theta), const int kdx, 
                               const int kdy, const int kcurv) 
    {
        int i, j;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j) 
                data_(i,j,kcurv) *= w(atan2(data_(i,j,kdy),data_(i,j,kdx)));
        bc->Apply(kcurv);
    }


    void UniformMesh2D::Ds_zero(const int kf, const int kdx, const int kdy,
                                const int kds) 
    {
        int i, j;
        double norm, f_x, f_y;
   
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j) {
                norm = sqrt(sqr(data_(i,j,kdx))+sqr(data_(i,j,kdy)));
                f_x = (data_(i+1,j,kf)-data_(i-1,j,kf))/2./dx;
                f_y = (data_(i,j+1,kf)-data_(i,j-1,kf))/2./dy;
                data_(i,j,kds) = (-f_x*data_(i,j,kdy)+f_y*data_(i,j,kdx))/norm;
            }
        bc->Apply(kds);
    }

    void UniformMesh2D::Ds_zero_4(const int kf, const int kdx, const int kdy,
                                  const int kds) 
    {
        int i, j;
        double norm, f_x, f_y;
   
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j) {
                norm = sqrt(sqr(data_(i,j,kdx))+sqr(data_(i,j,kdy)));
                f_x = (-data_(i+2,j,kf)+8*data_(i+1,j,kf)
                       -8*data_(i-1,j,kf)+data_(i-2,j,kf))/12./dx;
                f_y = (-data_(i,j+2,kf)+8*data_(i,j+1,kf)
                       -8*data_(i,j-1,kf)+data_(i,j-2,kf))/12./dy;
                data_(i,j,kds) = (-f_x*data_(i,j,kdy)+f_y*data_(i,j,kdx))/norm;
            }
        bc->Apply(kds);
    }

    void UniformMesh2D::Dss_zero(const int kf, const int kdx, const int kdy, 
                                 const int kdss)
    {
        double fx, fy, fxx, fyy, fxy, grad;
        int i, j;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j) {
                fx = (data_(i+1,j,kf)-data_(i-1,j,kf))/2./dx;
                fy = (data_(i,j+1,kf)-data_(i,j-1,kf))/2./dy;
                fxx = (data_(i+1,j,kf)-2*data_(i,j,kf)+data_(i-1,j,kf))/dx/dx;
                fyy = (data_(i,j+1,kf)-2*data_(i,j,kf)+data_(i,j-1,kf))/dy/dy;
                fxy = (data_(i-1,j-1,kf)-data_(i+1,j-1,kf)
                       -data_(i-1,j+1,kf)+data_(i+1,j+1,kf))/4/dx/dy;
                grad = sqrt(sqr(data_(i,j,kdx))+sqr(data_(i,j,kdy)));
                data_(i,j,kdss) = ((- fxx*sqr(data_(i,j,kdy))
                                    + 2.*fxy*data_(i,j,kdx)*data_(i,j,kdy)
                                    - fyy*sqr(data_(i,j,kdx)))/grad
                                   + (fx*data_(i,j,kdx)+fy*data_(i,j,kdy))
                                   *data_(i,j,kf))/grad;
            }
        bc->Apply(kdss);
    }

    void UniformMesh2D::Dss_zero_4(const int kf, const int kdx, const int kdy, 
                                   const int kdss)
    {
        double fx, fy, fxx, fyy, fxy, grad;
        int i, j;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j) {
                fx = (-data_(i+2,j,kf)+8*data_(i+1,j,kf)
                      -8*data_(i-1,j,kf)+data_(i-2,j,kf))/12./dx;
                fy = (-data_(i,j+2,kf)+8*data_(i,j+1,kf)
                      -8*data_(i,j-1,kf)+data_(i,j-2,kf))/12./dy;
                fxx = (-data_(i+2,j,kf)+16*data_(i+1,j,kf)
                       -30*data_(i,j,kf)+16*data_(i-1,j,kf)
                       -data_(i-2,j,kf))/12./dx/dx;
                fyy = (-data_(i,j+2,kf)+16*data_(i,j+1,kf)
                       -30*data_(i,j,kf) +16*data_(i,j-1,kf)
                       -data_(i,j-2,kf))/12./dy/dy;
                fxy = (-data_(i-2,j-2,kf)+16*data_(i-1,j-1,kf)
                       +data_(i+2,j-2,kf)-16*data_(i+1,j-1,kf)
                       +data_(i-2,j+2,kf)-16*data_(i-1,j+1,kf)
                       -data_(i+2,j+2,kf)+16*data_(i+1,j+1,kf))/48/dx/dy;
                grad = sqrt(sqr(data_(i,j,kdx))+sqr(data_(i,j,kdy)));
                data_(i,j,kdss) = ((- fxx*sqr(data_(i,j,kdy))
                                    + 2.*fxy*data_(i,j,kdx)*data_(i,j,kdy)
                                    - fyy*sqr(data_(i,j,kdx)))/grad
                                   + (fx*data_(i,j,kdx)+fy*data_(i,j,kdy))
                                   *data_(i,j,kf))/grad;
            }
        bc->Apply(kdss);
    }

    void UniformMesh2D::Dx_zero(const int k, const int kdx)
    {
        int i, j;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                data_(i,j,kdx) = (data_(i+1,j,k)-data_(i-1,j,k))/2./dx;
        bc->Apply(kdx);
    }

    void UniformMesh2D::Dx_minus(const int k, const int kdx)
    {
        int i, j;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                data_(i,j,kdx) = (data_(i,j,k)-data_(i-1,j,k))/dx;
        bc->Apply(kdx);
    }

    void UniformMesh2D::Dx_plus(const int k, const int kdx)
    {
        int i, j;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                data_(i,j,kdx) = (data_(i+1,j,k)-data_(i,j,k))/dx;
        bc->Apply(kdx);
    }

    void UniformMesh2D::Dx_upwind(const int v, const int k, const int kdx)
    {
        double dxp, dxm;
        int i, j;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j) {
                dxp = (data_(i+1,j,k)-data_(i,j,k))/dx;
                dxm = (data_(i,j,k)-data_(i-1,j,k))/dx;
                if (data_(i,j,v) >= 0.) 
                    data_(i,j,kdx) = min(dxm,0.)+max(dxp,0.);
                else
                    data_(i,j,kdx) = max(dxm,0.)+min(dxp,0.);
            }
        bc->Apply(kdx);
    }
   

    void UniformMesh2D::Dxx_zero(const int k, const int kdx)
    {
        int i, j;
			 for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                data_(i,j,kdx) = (data_(i+1,j,k)-2.*data_(i,j,k)
                                  +data_(i-1,j,k))/dx/dx;
        bc->Apply(kdx);
    }

    void UniformMesh2D::Dy_zero(const int k, const int kdy)
    {
        int i, j;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                data_(i,j,kdy) = (data_(i,j+1,k)-data_(i,j-1,k))/2./dy;
        bc->Apply(kdy);
    }

    void UniformMesh2D::Dy_minus(const int k, const int kdy)
    {
        int ii, j;
        for (ii=0; ii<maxi; ++ii)
            for (j=0; j<maxj; ++j)
                data_(ii,j,kdy) = (data_(ii,j,k)-data_(ii,j-1,k))/dy;
        bc->Apply(kdy);
    }

    void UniformMesh2D::Dy_plus(const int k, const int kdy)
    {
        int ii, j;
        for (ii=0; ii<maxi; ++ii)
            for (j=0; j<maxj; ++j)
                data_(ii,j,kdy) = (data_(ii,j+1,k)-data_(ii,j,k))/dy;
        bc->Apply(kdy);
    }

    void UniformMesh2D::Dy_upwind(const int v, const int k, const int kdy)
    {
        double dyp, dym;
        int ii, j;
        for (ii=0; ii<maxi; ++ii)
            for (j=0; j<maxj; ++j) {
                dyp = (data_(ii,j+1,k)-data_(ii,j,k))/dy;
                dym = (data_(ii,j,k)-data_(ii,j-1,k))/dy;
                if (data_(ii,j,v) >= 0.) 
                    data_(ii,j,kdy) = min(dym,0.)+max(dyp,0.);
                else
                    data_(ii,j,kdy) = max(dym,0.)+min(dyp,0.);
            }
        bc->Apply(kdy);
    }

    void UniformMesh2D::Dyy_zero(const int k, const int kdy)
    {
        int ii, j;
        for (ii=0; ii<maxi; ++ii)
            for (j=0; j<maxj; ++j)
                data_(ii,j,kdy) = (data_(ii,j+1,k)-2.*data_(ii,j,k)
                                   +data_(ii,j-1,k))/dy/dy;
        bc->Apply(kdy);
    }

    void UniformMesh2D::Dxy_zero(const int k, const int kdxy)
    {
        int ii, j;
        for (ii=0; ii<maxi; ++ii)
            for (j=0; j<maxj; ++j)
                data_(ii,j,kdxy) = (data_(ii-1,j-1,k)-data_(ii+1,j-1,k)
                                    -data_(ii-1,j+1,k)+data_(ii+1,j+1,k))/4./dx/dy;
        bc->Apply(kdxy);
    }

    void UniformMesh2D::Dx_zero_4(const int k, const int kdx)
    {
        int i, j;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                data_(i,j,kdx) = (-data_(i+2,j,k)+8*data_(i+1,j,k)
                                  -8*data_(i-1,j,k)+data_(i-2,j,k))/12./dx;
        bc->Apply(kdx);
    }

    void UniformMesh2D::Dx_minus_4(const int k, const int kdx)
    {
        int i, j;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                data_(i,j,kdx) = (3*data_(i,j,k)-4*data_(i-1,j,k)
                                  +data_(i-2,j,k))/2./dx;
        bc->Apply(kdx);
    }

    void UniformMesh2D::Dx_plus_4(const int k, const int kdx)
    {
        int i, j;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                data_(i,j,kdx) = (-data_(i+2,j,k)+4*data_(i+1,j,k)
                                  -3*data_(i,j,k))/2./dx;
        bc->Apply(kdx);
    }

    void UniformMesh2D::Dx_upwind_4(const int v, const int kdxm, const int kdxp,
                                    const int kdx)
    {
        int i, j;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j) {
                if (data_(i,j,v) >= 0.) 
                    data_(i,j,kdx) = min(data_(i,j,kdxm),0.)
                        +max(data_(i,j,kdxp),0.);
                else
                    data_(i,j,kdx) = max(data_(i,j,kdxm),0.)
                        +min(data_(i,j,kdxp),0.);
            }
        bc->Apply(kdx);
    }
   

    void UniformMesh2D::Dxx_zero_4(const int k, const int kdx)
    {
        int i, j;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                data_(i,j,kdx) = (-data_(i+2,j,k)+16*data_(i+1,j,k)
                                  -30*data_(i,j,k)+16*data_(i-1,j,k)
                                  -data_(i-2,j,k))/12./dx/dx;
        bc->Apply(kdx);
    }

    void UniformMesh2D::Dy_zero_4(const int k, const int kdy)
    {
        int i, j;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                data_(i,j,kdy) = (-data_(i,j+2,k)+8*data_(i,j+1,k)
                                  -8*data_(i,j-1,k)+data_(i,j-2,k))/12./dy;
        bc->Apply(kdy);
    }

    void UniformMesh2D::Dy_minus_4(const int k, const int kdy)
    {
        int i, j;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                data_(i,j,kdy) = (3*data_(i,j,k)-4*data_(i,j-1,k)
                                  +data_(i,j-2,k))/2./dy;
        bc->Apply(kdy);
    }

    void UniformMesh2D::Dy_plus_4(const int k, const int kdy)
    {
        int i, j;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                data_(i,j,kdy) = (-data_(i,j+2,k)+4*data_(i,j+1,k)
                                  -3*data_(i,j,k))/2./dy;
        bc->Apply(kdy);
    }

    void UniformMesh2D::Dy_upwind_4(const int v, const int kdym, const int kdyp,
                                    const int kdy)
    {
        int i, j;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j) {
                if (data_(i,j,v) >= 0.) 
                    data_(i,j,kdy) = min(data_(i,j,kdym),0.)
                        +max(data_(i,j,kdyp),0.);
                else
                    data_(i,j,kdy) = max(data_(i,j,kdym),0.)
                        +min(data_(i,j,kdyp),0.);
            }
        bc->Apply(kdy);
    }

    void UniformMesh2D::Dyy_zero_4(const int k, const int kdy)
    {
        int i, j;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                data_(i,j,kdy) = (-data_(i,j+2,k)+16*data_(i,j+1,k)
                                  -30*data_(i,j,k)+16*data_(i,j-1,k)
                                  -data_(i,j-2,k))/12./dy/dy;
        bc->Apply(kdy);
    }

    void UniformMesh2D::Dxy_zero_4(const int k, const int kdxy)
    {
        int i, j;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                data_(i,j,kdxy) = (-data_(i-2,j-2,k)+16*data_(i-1,j-1,k)
                                   +data_(i+2,j-2,k)-16*data_(i+1,j-1,k)
                                   +data_(i-2,j+2,k)-16*data_(i-1,j+1,k)
                                   -data_(i+2,j+2,k)+16*data_(i+1,j+1,k))/48/dx/dy;
        bc->Apply(kdxy);
    }

// Single node computations

    double UniformMesh2D::Dx_zero(const int i, const int j, const int k) const
    {
        return (data_(i+1,j,k)-data_(i-1,j,k))/2./dx;
    }

    double UniformMesh2D::Dx_minus(const int i, const int j, const int k) const
    {
        return (data_(i,j,k)-data_(i-1,j,k))/dx;
    }

    double UniformMesh2D::Dx_plus(const int i, const int j, const int k) const
    {
        return (data_(i+1,j,k)-data_(i,j,k))/dx;
    }

    double UniformMesh2D::Dx_upwind(const int i, const int j, const int k,
                                    const int kv) const
    {
        double dxp = (data_(i+1,j,k)-data_(i,j,k))/dx;
        double dxm = (data_(i,j,k)-data_(i-1,j,k))/dx;
        return (data_(i,j,kv) >= 0.) ? 
            min(dxm,0.)+max(dxp,0.) : max(dxm,0.)+min(dxp,0.);
    }
   

    double UniformMesh2D::Dxx_zero(const int i, const int j, const int k) const
    {
        return (data_(i+1,j,k)-2.*data_(i,j,k)+data_(i-1,j,k))/dx/dx;
    }

    double UniformMesh2D::Dy_zero(const int i, const int j, const int k) const
    {
        return (data_(i,j+1,k)-data_(i,j-1,k))/2./dy;
    }

    double UniformMesh2D::Dy_minus(const int i, const int j, const int k) const
    {
        return (data_(i,j,k)-data_(i,j-1,k))/dy;
    }

    double UniformMesh2D::Dy_plus(const int i, const int j, const int k) const
    {
        return (data_(i,j+1,k)-data_(i,j,k))/dy;
    }

    double UniformMesh2D::Dy_upwind(const int i, const int j, const int k,
                                    const int kv) const
    {
        double dyp = (data_(i,j+1,k)-data_(i,j,k))/dy;
        double dym = (data_(i,j,k)-data_(i,j-1,k))/dy;
        return (data_(i,j,kv) >= 0.) ?
            min(dym,0.)+max(dyp,0.) : max(dym,0.)+min(dyp,0.);
    }

    double UniformMesh2D::Dyy_zero(const int i, const int j, const int k) const
    {
        return (data_(i,j+1,k)-2.*data_(i,j,k)+data_(i,j-1,k))/dy/dy;
    }

    double UniformMesh2D::Dxy_zero(const int i, const int j, const int k) const
    {
        return (data_(i-1,j-1,k)-data_(i+1,j-1,k)
                -data_(i-1,j+1,k)+data_(i+1,j+1,k))/4./dx/dy;
    }

    double UniformMesh2D::Dx_zero_4(const int i, const int j, const int k) const
    {
        return (-data_(i+2,j,k)+8*data_(i+1,j,k)
                -8*data_(i-1,j,k)+data_(i-2,j,k))/12./dx;
    }

    double UniformMesh2D::Dx_minus_4(const int i, const int j, const int k) const
    {
        return (3*data_(i,j,k)-4*data_(i-1,j,k)+data_(i-2,j,k))/2./dx;
    }

    double UniformMesh2D::Dx_plus_4(const int i, const int j, const int k) const
    {
        return (-data_(i+2,j,k)+4*data_(i+1,j,k)-3*data_(i,j,k))/2./dx;
    }

#if 0
    double UniformMesh2D::Dx_upwind_4(const int i, const int j, const int k,
                                      const int kv) const
    {
        double dxm = Dx_minus_4(i,j,k);
        double dxp = Dx_plus_4(i,j,k);
        return (data_(i,j,kv) >= 0.) ?
            min(dxm,0.)+max(dxp,0.) : max(dxm,0.)+min(dxp,0.);
    }
#endif

    double UniformMesh2D::Dxx_zero_4(const int i, const int j, const int k) const
    {
        return (-data_(i+2,j,k)+16*data_(i+1,j,k)-30*data_(i,j,k)
                +16*data_(i-1,j,k)-data_(i-2,j,k))/12./dx/dx;
    }

#if 0
    double UniformMesh2D::Dxx_plus_4(const int i, const int j, const int k) const
    {
        return (11.*data_(i+4,j,k)-56.*data_(i+3,j,k)+144.*data_(i+2,j,k)
                -104.*data_(i+1,j,k)+35.data_(i,j,k))/12./dy;
    }

    double UniformMesh2D::Dxx_minus_4(const int i, const int j, const int k) const
    {
        return (11.*data_(i-4,j,k)-56.*data_(i-3,j,k)+144.*data_(i-2,j,k)
                -104.*data_(i-1,j,k)+35.data_(i,j,k))/12./dy;
    }
#endif

    double UniformMesh2D::Dy_zero_4(const int i, const int j, const int k) const
    {
        return (-data_(i,j+2,k)+8*data_(i,j+1,k)
                -8*data_(i,j-1,k)+data_(i,j-2,k))/12./dy;
    }

#if 1
    double UniformMesh2D::Dy_minus_4(const int i, const int j, const int k) const
    {
        return (3*data_(i,j,k)-4*data_(i,j-1,k)+data_(i,j-2,k))/2./dy;
    }

    double UniformMesh2D::Dy_plus_4(const int i, const int j, const int k) const
    {
        return (-data_(i,j+2,k)+4*data_(i,j+1,k)-3*data_(i,j,k))/2./dy;
    }
#endif

#if 0
    double UniformMesh2D::Dy_upwind_4(const int i, const int j, const int k,
                                      const int kv) const
    {
        double dym = Dy_minus_4(i,j,k);
        double dyp = Dy_plus_4(i,j,k);
        return (data_(i,j,kv) >= 0.) ?
            min(dym,0.)+max(dyp,0.) : max(dym,0.)+min(dyp,0.);
    }
#endif

    double UniformMesh2D::Dyy_zero_4(const int i, const int j, const int k) const
    {
        return (-data_(i,j+2,k)+16*data_(i,j+1,k)-30*data_(i,j,k)
                +16*data_(i,j-1,k)-data_(i,j-2,k))/12./dy/dy;
    }

#if 0
    double UniformMesh2D::Dyy_plus_4(const int i, const int j, const int k) const
    {
        return (11.*data_(i,j+4,k)-56.*data_(i,j+3,k)+144.*data_(i,j+2,k)
                -104.*data_(i,j+1,k)+35.data_(i,j,k))/12./dy;
    }

    double UniformMesh2D::Dyy_minus_4(const int i, const int j, const int k) const
    {
        return (11.*data_(i,j-4,k)-56.*data_(i,j-3,k)+144.*data_(i,j-2,k)
                -104.*data_(i,j-1,k)+35.data_(i,j,k))/12./dy;
    }
#endif

    double UniformMesh2D::Dxy_zero_4(const int i, const int j, const int k) const
    {
        return (-data_(i-2,j-2,k)+16*data_(i-1,j-1,k)
                +data_(i+2,j-2,k)-16*data_(i+1,j-1,k)
                +data_(i-2,j+2,k)-16*data_(i-1,j+1,k)
                -data_(i+2,j+2,k)+16*data_(i+1,j+1,k))/48/dx/dy;
    }

#if 0
    double UniformMesh2D::Dxy_plusplus_4(const int i, const int j, const int k) const
    {
        return (254*data_(i,j,k)-396*data_(i+1,j,k)+189*data_(i+2,j,k)-56*data_(i+3,j,k)
                +9*data_(i+4,j,k)-396*data_(i,j+1,k)+456*data_(i+1,j+1,k)
                -72*data_(i+3,j+1,k)+12*data_(i+4,j+1,k)+189*data_(i,j+2,k)
                -378*data_(i+2,j+2,k)+216*data_(i+3,j+2,k)-27*data_(i+4,j+2,k)
                -56*data_(i,j+3,k)-72*data_(i+1,j+3,k)+216*data_(i+2,j+3,k)
                -88*data_(i+3,j+3,k)+9*data_(i,j+4,k)+12*data_(i+1,j+4,k)
                -27*data_(i+2,j+4,k)+6*data_(i+4,j+4,k))/72/dx/dy;
    }

    double UniformMesh2D::Dxy_plusminus_4(const int i, const int j, const int k) const
    {
        return (-data_(i-2,j-2,k)+16*data_(i-1,j-1,k)
                +data_(i+2,j-2,k)-16*data_(i+1,j-1,k)
                +data_(i-2,j+2,k)-16*data_(i-1,j+1,k)
                -data_(i+2,j+2,k)+16*data_(i+1,j+1,k))/48/dx/dy;
    }

    double UniformMesh2D::Dxy_minusplus_4(const int i, const int j, const int k) const
    {
        return (-data_(i-2,j-2,k)+16*data_(i-1,j-1,k)
                +data_(i+2,j-2,k)-16*data_(i+1,j-1,k)
                +data_(i-2,j+2,k)-16*data_(i-1,j+1,k)
                -data_(i+2,j+2,k)+16*data_(i+1,j+1,k))/48/dx/dy;
    }

    double UniformMesh2D::Dxy_minusminus_4(const int i, const int j, const int k) const
    {
        return (-data_(i-2,j-2,k)+16*data_(i-1,j-1,k)
                +data_(i+2,j-2,k)-16*data_(i+1,j-1,k)
                +data_(i-2,j+2,k)-16*data_(i-1,j+1,k)
                -data_(i+2,j+2,k)+16*data_(i+1,j+1,k))/48/dx/dy;
    }

    double UniformMesh2D::Dx_plusplus_4(const int i, const int j, const double ifrac,
                                        const double jfrac, const int k) const
    {
        double w[10];
        w[0] = 1.;
        w[1] = ifrac;
        w[2] = jfrac;
        w[3] = w[1]*w[1];
        w[4] = w[1]*w[2];
        w[5] = w[2]*w[2];
        w[6] = w[1]*w[3];
        w[7] = w[2]*w[3];
        w[8] = w[1]*w[5];
        w[9] = w[2]*w[5];
        return ((-22*w[-1]+24+35*w[1]-6*w[2]-30*w[3]-15*w[4]+6*w[6]+6*w[7]+2*w[8])*data_(i,j,k)
                +(36*w[-1]-60-52*w[1]+18*w[2]+72*w[3]+18*w[4]-18*w[6]-12*w[7]-2*w[8])*data_(i+1,j,k)
                +(-18*w[-1]+48+21*w[1]-18*w[2]-54*w[3]-3*w[4]+18*w[6]+6*w[7])*data_(i+2,j,k)
                +(4*w[-1]-12-4*w[1]+6*w[2]+12*w[3]-6*w[6])*data_(i+3,j,k)
                +(-52*w[1]+36*w[3]+36*w[4]-6*w[6]-12*w[7]-6*w[8])*data_(i,j+1,k)
                +(72*w[1]-84*w[3]-42*w[4]+18*w[6]+24*w[7]+6*w[8])*data_(i+1,j+1,k)
                +(-24*w[1]+60*w[3]+6*w[4]-18*w[6]-12*w[7])*data_(i+2,j+1,k)
                +(4*w[1]-12*w[3]+6*w[6])*data_(i+3,j+1,k)
                +(21*w[1]-6*w[3]-27*w[4]+6*w[7]+6*w[8])*data_(i,j+2,k)
                +(-24*w[1]+12*w[3]+30*w[4]-12*w[7]-6*w[8])*data_(i+1,j+2,k)
                +(3*w[1]-6*w[3]-3*w[4]+6*w[7])*data_(i+2,j+2,k)
                +(-4*w[1]+6*w[4]-2*w[8])*data_(i,j+3,k)
                +(4*w[1]-6*w[4]+2*w[8])*data_(i+1,j+3,k))/12/dx;
    }
                        
    double UniformMesh2D::Dx_plusminus_4(const int i, const int j, const double ifrac,
                                         const double jfrac, const int k) const
    {
        double w[10];
        w[0] = 1.;
        w[1] = ifrac;
        w[2] = jfrac;
        w[3] = w[1]*w[1];
        w[4] = w[1]*w[2];
        w[5] = w[2]*w[2];
        w[6] = w[1]*w[3];
        w[7] = w[2]*w[3];
        w[8] = w[1]*w[5];
        w[9] = w[2]*w[5];
        return ()/12/dx;
    }


#endif
}
