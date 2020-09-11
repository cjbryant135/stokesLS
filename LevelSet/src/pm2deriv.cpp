//
//  pm2deriv.cpp
//  LevelSet
//
//  Created by David Chopp on 7/25/19.
//  Copyright © 2019 David Chopp. All rights reserved.
//

#include "polarmesh2d.h"
#include "pm2boundary.h"


namespace levelset {
    
    void PolarMesh2D::Dr_zero(const int k, const int kdr)
    {
        int i, j;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                data_(i,j,kdr) = dr[i-1]/dr[i]/(dr[i-1]+dr[i])*data_(i+1,j,k)
                                + (dr[i]-dr[i-1])/dr[i]/dr[i-1]*data_(i,j,k)
                                - dr[i]/dr[i-1]/(dr[i-1]+dr[i])*data_(i-1,j,k);
        bc->Apply(kdr);
    }

    void PolarMesh2D::Dr_minus(const int k, const int kdr)
    {
        int i, j;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                data_(i,j,kdr) = (data_(i,j,k)-data_(i-1,j,k))/dr[i-1];
        bc->Apply(kdr);
    }
    
    void PolarMesh2D::Dr_plus(const int k, const int kdr)
    {
        int i, j;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                data_(i,j,kdr) = (data_(i+1,j,k)-data_(i,j,k))/dr[i];
        bc->Apply(kdr);
    }
    
    void PolarMesh2D::Dr_upwind(const int v, const int k, const int kdr)
    {
        double drp, drm;
        int i, j;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j) {
                drp = (data_(i+1,j,k)-data_(i,j,k))/dr[i];
                drm = (data_(i,j,k)-data_(i-1,j,k))/dr[i-1];
                if (data_(i,j,v) >= 0.)
                    data_(i,j,kdr) = min(drm,0.)+max(drp,0.);
                else
                    data_(i,j,kdr) = max(drm,0.)+min(drp,0.);
            }
        bc->Apply(kdr);
    }
    
    void PolarMesh2D::Drr_zero(const int k, const int kdr)
    {
        int i, j;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
				//Original version. I think this is wrong
//                data_(i,j,kdr) = 2.*(data_(i+1,j,k)/dr[i-1]/(dr[i-1]-dr[i])
//                                     -data_(i,j,k)/dr[i-1]/dr[i]
//                                     +data_(i-1,j,k)/dr[i]/(dr[i-1]+dr[i]));
                //new version 
				 	data_(i,j,kdr) = 2.*(data_(i+1,j,k)/dr[i]/(dr[i-1]+dr[i])
                                    -data_(i,j,k)/dr[i-1]/dr[i]
                                    +data_(i-1,j,k)/dr[i-1]/(dr[i-1]+dr[i]));
		  bc->Apply(kdr);
    }
    
    void PolarMesh2D::Dth_zero(const int k, const int kdth)
    {
        int i, j;
        int jM, jP;
		  for (i=0; i<maxi; ++i) {
            for (j=0; j<maxj; ++j) {
                jM = j > 0 ? j-1 : maxj-1;
					 jP = j < maxj-1 ? j+1 : 0;
					 data_(i,j,kdth) = (data_(i,jP,k)-data_(i,jM,k))/2./dtheta;
        	   }
		  }
		  
		  bc->Apply(kdth);
    }
    
    void PolarMesh2D::Dth_minus(const int k, const int kdth)
    {
        int ii, j, jM;
        for (ii=0; ii<maxi; ++ii) {
            for (j=0; j<maxj; ++j) {
                jM = j > 0 ? j-1 : maxj-1;
					 data_(ii,j,kdth) = (data_(ii,j,k)-data_(ii,jM,k))/dtheta;
        		}
		  }
		  bc->Apply(kdth);
    }
    
    void PolarMesh2D::Dth_plus(const int k, const int kdth)
    {
        int ii, j, jP;
        for (ii=0; ii<maxi; ++ii) {
            for (j=0; j<maxj; ++j) {
					 jP = j < maxj-1 ? j+1 : 0;
                data_(ii,j,kdth) = (data_(ii,jP,k)-data_(ii,j,k))/dtheta;
        		}
		  }
		  bc->Apply(kdth);
    }
    
    void PolarMesh2D::Dth_upwind(const int v, const int k, const int kdth)
    {
        double dthp, dthm;
        int ii, j, jP, jM;
        for (ii=0; ii<maxi; ++ii)
            for (j=0; j<maxj; ++j) {
                jM = j > 0 ? j-1 : maxj-1;
					 jP = j < maxj-1 ? j+1 : 0;
					 dthp = (data_(ii,jP,k)-data_(ii,j,k))/dtheta;
                dthm = (data_(ii,j,k)-data_(ii,jM,k))/dtheta;
                if (data_(ii,j,v) >= 0.)
                    data_(ii,j,kdth) = min(dthm,0.)+max(dthp,0.);
                else
                    data_(ii,j,kdth) = max(dthm,0.)+min(dthp,0.);
            }
        bc->Apply(kdth);
    }
    
    void PolarMesh2D::Dthth_zero(const int k, const int kdth)
    {
        int ii, j, jM, jP;
        for (ii=0; ii<maxi; ++ii)
            for (j=0; j<maxj; ++j) {
                jM = j > 0 ? j-1 : maxj-1;
					 jP = j < maxj-1 ? j+1 : 0;
                data_(ii,j,kdth) = (data_(ii,jP,k)-2.*data_(ii,j,k)
                                   +data_(ii,jM,k))/dtheta/dtheta;
        		}
		  bc->Apply(kdth);
    }
    
    void PolarMesh2D::Drth_zero(const int k, const int kdrth)
    {
        int i, j, jM, jP;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j) {
//                data_(i,j,kdrth) = (dr[i-1]/dr[i]/(dr[i-1]+dr[i])*data_(i+1,j+1,k)
//                                    + (dr[i]-dr[i-1])/dr[i]/dr[i-1]*data_(i,j+1,k)
//                                    - dr[i]/dr[i-1]/(dr[i-1]+dr[i])*data_(i-1,j+1,k)
//                                    - dr[i-1]/dr[i]/(dr[i-1]+dr[i])*data_(i+1,j-1,k)
//                                    - (dr[i]-dr[i-1])/dr[i]/dr[i-1]*data_(i,j-1,k)
//                                    + dr[i]/dr[i-1]/(dr[i-1]+dr[i])*data_(i-1,-1,k))/2./dtheta;
        				//new version
                	jM = j > 0 ? j-1 : maxj-1;
					 	jP = j < maxj-1 ? j+1 : 0;
						data_(i,j,kdrth) = (data_(i+1,jP,k) + data_(i-1,jM,k) - 
													data_(i+1,jM,k) - data_(i-1,jP,k) )/2./dtheta/(dr[i-1]+dr[i]);
		  		}
		  bc->Apply(kdrth);
    }
    
    void PolarMesh2D::Dr_zero_4(const int k, const int kdr)
    {
        int i, j;
        double a[5];
        for (i=0; i<maxi; ++i) {
            a[0] = dr[i-1]*(dr[i-1]+dr[i-2])*dr[i]/dr[i+1]/(dr[i]+dr[i+1])/(dr[i-1]+dr[i-2]+dr[i]+dr[i+1]);
            a[1] = dr[i-1]*(dr[i-1]+dr[i-2])*(dr[i]+dr[i+1])/dr[i]/(dr[-1]+dr[i])/(dr[i-1]+dr[i-2]+dr[i+1])/dr[i+1];
            a[2] = ((2*dr[i-1]+dr[i-2])*dr[i]*(dr[i]+dr[i+1])-dr[i-1]*(dr[i-1]+dr[i-2])*(2*dr[i]+dr[i+1]))
                    /dr[i-1]/dr[i]/(dr[i-1]+dr[i-2])/(dr[i]+dr[i+1]);
            a[3] = -(dr[i-1]+dr[i-2])*dr[i]*(dr[i]+dr[i+1])/dr[i-1]/dr[i-2]/(dr[i-1]+dr[i])/(dr[i-1]+dr[i]+dr[i+1]);
            a[4] = dr[i-1]*dr[i]*(dr[i]+dr[i+1])/dr[i-2]/(dr[i-1]+dr[i-2])/(dr[i-1]+dr[i-2]+dr[i])/(dr[i-1]+dr[i-2]+dr[i]+dr[i+1]);
            for (j=0; j<maxj; ++j)
                data_(i,j,kdr) = a[0]*data_(i+2,j,k)
                                + a[1]*data_(i+1,j,k)
                                + a[2]*data_(i,j,k)
                                + a[3]*data_(i-1,j,k)
                                + a[4]*data_(i-2,j,k);
        }
        bc->Apply(kdr);
    }
    
    void PolarMesh2D::Dr_minus_4(const int k, const int kdr)
    {
        int i, j;
        double a[3];
        for (i=0; i<maxi; ++i) {
            a[0] = (2*dr[i-1]+dr[i-2])/dr[i-1]/(dr[i-1]+dr[i-2]);
            a[1] = -(dr[i-1]+dr[i-2])/dr[i-1]/dr[i-2];
            a[2] = dr[i-1]/dr[i-2]/(dr[i-1]+dr[i-2]);
            for (j=0; j<maxj; ++j)
                data_(i,j,kdr) = a[0]*data_(i,j,k)
                                  +a[1]*data_(i-1,j,k)
                                  +a[2]*data_(i-2,j,k);
        }
        bc->Apply(kdr);
    }
    
    void PolarMesh2D::Dr_plus_4(const int k, const int kdr)
    {
        int i, j;
        double a[3];
        for (i=0; i<maxi; ++i) {
            a[0] = -dr[i]/dr[i+1]/(dr[i]+dr[i+1]);
            a[1] = (dr[i]+dr[i+1])/dr[i]/dr[i+1];
            a[2] = -(2*dr[i]+dr[i+1])/dr[i]/(dr[i]+dr[i+1]);
            for (j=0; j<maxj; ++j)
                data_(i,j,kdr) = a[0]*data_(i+2,j,k)
                                  +a[1]*data_(i+1,j,k)
                                  +a[2]*data_(i,j,k);
        }
        bc->Apply(kdr);
    }
    
    void PolarMesh2D::Dr_upwind_4(const int v, const int kdrm, const int kdrp,
                                    const int kdr)
    {
        int i, j;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j) {
                if (data_(i,j,v) >= 0.)
                    data_(i,j,kdr) = min(data_(i,j,kdrm),0.)
                    +max(data_(i,j,kdrp),0.);
                else
                    data_(i,j,kdr) = max(data_(i,j,kdrm),0.)
                    +min(data_(i,j,kdrp),0.);
            }
        bc->Apply(kdr);
    }
    
    
    void PolarMesh2D::Drr_zero_4(const int k, const int kdr)
    {
        int i, j;
        double a[5];
        for (i=0; i<maxi; ++i) {
            a[0] = (2*(dr[i-1]*dr[i-1] + dr[i-1]*dr[i-2] - 2*dr[i-1]*dr[i] - dr[i-2]*dr[i]))/
                    (dr[i+1]*(dr[i] + dr[i+1])*(dr[i-1] + dr[i] + dr[i+1])*(dr[i-1] + dr[i-2] + dr[i] + dr[i+1]));
            a[1] = (2*(-dr[i-1]*dr[i-1] - dr[i-1]*dr[i-2] + 2*dr[i-1]*dr[i] + dr[i-2]*dr[i] + 2*dr[i-1]*dr[i+1] + dr[i-2]*dr[i+1]))/
                    (dr[i]*(dr[i-1] + dr[i])*(dr[i-1] + dr[i-2] + dr[i])*dr[i+1]);
            a[2] = (2*(dr[i-1]*dr[i-1] + dr[i-1]*dr[i-2] - 4*dr[i-1]*dr[i] - 2*dr[i-2]*dr[i] + dr[i]*dr[i] - 2*dr[i-1]*dr[i+1] -
                       dr[i-2]*dr[i+1] + dr[i]*dr[i+1]))/(dr[i-1]*(dr[i-1] + dr[i-2])*dr[i]*(dr[i] + dr[i+1]));
            a[3] = (2*(-(dr[i]*(dr[i] + dr[i+1])) + dr[i-1]*(2*dr[i] + dr[i+1]) + dr[i-2]*(2*dr[i] + dr[i+1])))/
                    (dr[i-1]*dr[i-2]*(dr[i-1] + dr[i])*(dr[i-1] + dr[i] + dr[i+1]));
            a[4] = (-2*(2*dr[i-1]*dr[i] - dr[i]*dr[i] + dr[i-1]*dr[i+1] - dr[i]*dr[i+1]))/
                    (dr[i-2]*(dr[i-1] + dr[i-2])*(dr[i-1] + dr[i-2] + dr[i])*(dr[i-1] + dr[i-2] + dr[i] + dr[i+1]));
            for (j=0; j<maxj; ++j)
                data_(i,j,kdr) = a[0]*data_(i+2,j,k)+a[1]*data_(i+1,j,k)+a[2]*data_(i,j,k)
                                +a[3]*data_(i-1,j,k)+a[4]*data_(i-2,j,k);
        }
        bc->Apply(kdr);
    }
    
    void PolarMesh2D::Dth_zero_4(const int k, const int kdth)
    {
        int i, j;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                data_(i,j,kdth) = (-data_(i,j+2,k)+8*data_(i,j+1,k)
                                  -8*data_(i,j-1,k)+data_(i,j-2,k))/12./dtheta;
        bc->Apply(kdth);
    }
    
    void PolarMesh2D::Dth_minus_4(const int k, const int kdth)
    {
        int i, j;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                data_(i,j,kdth) = (3*data_(i,j,k)-4*data_(i,j-1,k)
                                  +data_(i,j-2,k))/2./dtheta;
        bc->Apply(kdth);
    }
    
    void PolarMesh2D::Dth_plus_4(const int k, const int kdth)
    {
        int i, j;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                data_(i,j,kdth) = (-data_(i,j+2,k)+4*data_(i,j+1,k)
                                  -3*data_(i,j,k))/2./dtheta;
        bc->Apply(kdth);
    }
    
    void PolarMesh2D::Dth_upwind_4(const int v, const int kdthm, const int kdthp,
                                    const int kdth)
    {
        int i, j;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j) {
                if (data_(i,j,v) >= 0.)
                    data_(i,j,kdth) = min(data_(i,j,kdthm),0.)
                    +max(data_(i,j,kdthp),0.);
                else
                    data_(i,j,kdth) = max(data_(i,j,kdthm),0.)
                    +min(data_(i,j,kdthp),0.);
            }
        bc->Apply(kdth);
    }
    
    void PolarMesh2D::Dthth_zero_4(const int k, const int kdth)
    {
        int i, j;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                data_(i,j,kdth) = (-data_(i,j+2,k)+16*data_(i,j+1,k)
                                  -30*data_(i,j,k)+16*data_(i,j-1,k)
                                  -data_(i,j-2,k))/12./dtheta/dtheta;
        bc->Apply(kdth);
    }
    
    void PolarMesh2D::Drth_zero_4(const int k, const int kdrth)
    {
        int i, j;
        double a[5];
        for (i=0; i<maxi; ++i) {
            a[0] = dr[i-1]*(dr[i-1]+dr[i-2])*dr[i]/dr[i+1]/(dr[i]+dr[i+1])/(dr[i-1]+dr[i-2]+dr[i]+dr[i+1]);
            a[1] = dr[i-1]*(dr[i-1]+dr[i-2])*(dr[i]+dr[i+1])/dr[i]/(dr[-1]+dr[i])/(dr[i-1]+dr[i-2]+dr[i+1])/dr[i+1];
            a[2] = ((2*dr[i-1]+dr[i-2])*dr[i]*(dr[i]+dr[i+1])-dr[i-1]*(dr[i-1]+dr[i-2])*(2*dr[i]+dr[i+1]))
                        /dr[i-1]/dr[i]/(dr[i-1]+dr[i-2])/(dr[i]+dr[i+1]);
            a[3] = -(dr[i-1]+dr[i-2])*dr[i]*(dr[i]+dr[i+1])/dr[i-1]/dr[i-2]/(dr[i-1]+dr[i])/(dr[i-1]+dr[i]+dr[i+1]);
            a[4] = dr[i-1]*dr[i]*(dr[i]+dr[i+1])/dr[i-2]/(dr[i-1]+dr[i-2])/(dr[i-1]+dr[i-2]+dr[i])/(dr[i-1]+dr[i-2]+dr[i]+dr[i+1]);
            for (j=0; j<maxj; ++j)
                data_(i,j,kdrth) = (-(a[0]*data_(i+2,j+2,k) +a[1]*data_(i+1,j+2,k) +a[2]*data_(i,j+2,k) +a[3]*data_(i-1,j+2,k) +a[4]*data_(i-2,j+2,k))
                                +8*(a[0]*data_(i+2,j+1,k) +a[1]*data_(i+1,j+1,k) +a[2]*data_(i,j+1,k) +a[3]*data_(i-1,j+1,k) +a[4]*data_(i-2,j+1,k))
                                -8*(a[0]*data_(i+2,j-1,k) +a[1]*data_(i+1,j-1,k) +a[2]*data_(i,j-1,k) +a[3]*data_(i-1,j-1,k) +a[4]*data_(i-2,j-1,k))
                                +(a[0]*data_(i+2,j-2,k) +a[1]*data_(i+1,j-2,k) +a[2]*data_(i,j-2,k) +a[3]*data_(i-1,j-2,k) +a[4]*data_(i-2,j-2,k)))
                                /12./dtheta;
        }
        bc->Apply(kdrth);
    }

    void PolarMesh2D::NormUpwindGrad(const int v, const int kdrm, const int kdrp,
                        const int kdtm, const int kdtp, const int knorm)
    { 
        int i, j;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                if (data_(i,j,v) > 0.) 
                    data_(i,j,knorm) = sqrt(sqr(max(data_(i,j,kdrp),-data_(i,j,kdrm), 0.))
                                            +sqr(max(data_(i,j,kdtp),-data_(i,j,kdtm), 0.)/r[i]));
                else
                    data_(i,j,knorm) = sqrt(sqr(max(-data_(i,j,kdrp), data_(i,j,kdrm), 0.))
                                            +sqr(max(-data_(i,j,kdtp),data_(i,j,kdtm), 0.)/r[i]));
        bc->Apply(knorm);

    }
    
    void PolarMesh2D::NormUpwindGrad(const int v, const int k, const int knorm, const int dir)
    {
        int i, j, jP, jM;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j) {
                jM = j > 0 ? j-1 : maxj-1;
					 jP = j < maxj-1 ? j+1 : 0;
                if (data_(i,j,v)*dir > 0.)
                    data_(i,j,knorm) = sqrt(sqr(max((data_(i+1,j,k)-data_(i,j,k))/dr[i],
                                                    -(data_(i,j,k)-data_(i-1,j,k))/dr[i-1],0.))
                                            +sqr(max((data_(i,jP,k)-data_(i,j,k))/dtheta,
                                                     -(data_(i,j,k)-data_(i,jM,k))/dtheta,0.)/r[i]));
                else
                    data_(i,j,knorm) = sqrt(sqr(max(-(data_(i+1,j,k)-data_(i,j,k))/dr[i],
                                                    (data_(i,j,k)-data_(i-1,j,k))/dr[i-1],0.))
                                            +sqr(max(-(data_(i,jP,k)-data_(i,j,k))/dtheta,
                                                     (data_(i,j,k)-data_(i,jM,k))/dtheta,0.)/r[i]));
        		}
		  bc->Apply(knorm);
        
    }
    
    void PolarMesh2D::NormGrad(const int kdr, const int kdtheta, const int knorm)
    {
        int i, j;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                data_(i,j,knorm) = sqrt(sqr(data_(i,j,kdr))
                                        +sqr(data_(i,j,kdtheta)/r[i]));
        bc->Apply(knorm);
    }
    void PolarMesh2D::Curvature(const int kdr, const int kdtheta, const int kdrr,
                   const int kdrtheta, const int kdtheta2, const int kcurv)
    {
        int i, j;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                //original version
					 data_(i,j,kcurv) = (r[i]*data_(i,j,kdrr)*sqr(data_(i,j,kdtheta))
                                    +r[i]*data_(i,j,kdtheta2)*sqr(data_(i,j,kdr))
                                    -2.*r[i]*data_(i,j,kdrtheta)*data_(i,j,kdr)*data_(i,j,kdtheta)
                                    +data_(i,j,kdr)*(sqr(r[i]*data_(i,j,kdr))+2*sqr(data_(i,j,kdtheta))))
                /(pow(sqr(r[i]*data_(i,j,kdr))+sqr(data_(i,j,kdtheta)),1.5)+DBL_EPSILON);
                
        bc->Apply(kcurv);
    	//printf("HERE\n");
	 }
 
	  //single node methods (all one sided difference methods for finding curvature at the particle surface)
	  double PolarMesh2D::Drr_plus(const int i, const int j, const int k) 
	  {
	    (data_(i+2,j,k)+dr[i+1]*data_(i,j,k)/dr[i] 
		 	-(1.+dr[i+1]/dr[i])*data_(i+1,j,k) )
		 	/(dr[i]*dr[i+1]+ (sqr(dr[i])+sqr(dr[i+1]))/2.);
	  }
	  double PolarMesh2D::Drth_plus_zero(const int i, const int j, const int k) 
	  {
	  	 int jP, jM;
		 jM = j > 0 ? j-1 : maxj-1;
		 jP = j < maxj-1 ? j+1 : 0;
		 double a[3];
	 	 a[0] = -dr[i]/dr[i+1]/(dr[i]+dr[i+1]);
		 a[1] = (dr[i]+dr[i+1])/dr[i]/dr[i+1];
		 a[2] = -(2*dr[i]+dr[i+1])/dr[i]/(dr[i]+dr[i+1]);
	 
		 return (a[0]*(data_(i+2,jP,k)-data_(i+2,jM,k))
		 		 	+a[1]*(data_(i+1,jP,k)-data_(i+1,jM,k))
					+a[2]*(data_(i,jP,k)-data_(i,jM,k)))/2./dtheta;
		 //first order in r
		 //return (data_(i+1,jP,k)+data_(i,jM,k)-data_(i+1,jM,k)-data_(i,jP,k))/2./dtheta/dr[i];

	  }
	  
	  double PolarMesh2D::Curvature(const int i, const int j, const double Dr, const double Dth, const double Drr, 
		  					 const double Drth, const double Dthth) 
	  { 
		 return (r[i]*Drr*sqr(Dth)
		 			+r[i]*Dthth*sqr(Dr)
					-2.*r[i]*Drth*Dr*Dth
					+Dr*(sqr(r[i]*Dr)+2.*sqr(Dth)) )
				/(pow(sqr(r[i]*Dr)+sqr(Dth),1.5)+DBL_EPSILON);
	  }
	
	  double PolarMesh2D::Dr_plus_4(const int i, const int j, const int k)
	  {
        double a[3];
			a[0] = -dr[i]/dr[i+1]/(dr[i]+dr[i+1]);
			a[1] = (dr[i]+dr[i+1])/dr[i]/dr[i+1];
			a[2] = -(2*dr[i]+dr[i+1])/dr[i]/(dr[i]+dr[i+1]);
			return a[0]*data_(i+2,j,k)
									+a[1]*data_(i+1,j,k)
									+a[2]*data_(i,j,k);
	 }




}




