//
//  polarmesh2d.cpp
//  LevelSet
//
//  Created by David Chopp on 7/23/19.
//  Copyright © 2019 David Chopp. All rights reserved.
//

#include "polarmesh2d.h"
#include "numtrait.h"
#include "debug.h"
#include "pm2boundary.h"

namespace levelset {
    
    PolarMesh2D::PolarMesh2D(STEPSIZE, const int m, const int n,
                             const int sd, const int si,
                             const double deltar,
                             const double r0, const double stretchinc,
                             PM2_Boundary& thebc)
    : maxi(m), maxj(n), maxk(0), nextw(0), r(NULL), dr(NULL),
    thedata(NULL), theidata(NULL), bc(&thebc)
    {
        SetStepSize(m, n, sd, si, deltar, r0, thebc);
    }
    
    PolarMesh2D::PolarMesh2D(BOUNDS, const int m, const int n,
                                 const int sd, const int si,
                                 const double r0, const double r1,
                                 const double stretchinc,
                                 PM2_Boundary& thebc)
    : maxi(m), maxj(n), maxk(0), nextw(0), r(NULL), dr(NULL), thedata(NULL), bc(&thebc)
    {
        SetBoundsSize(m, n, sd, si, r0, r1, thebc);
    }
    
    PolarMesh2D::~PolarMesh2D(void)
    {
        if (thedata) delete[] thedata;
        if (theidata) delete[] theidata;
        if (fmm_coi) delete[] fmm_coi;
        if (tr) delete[] tr;
        if (tdr) delete[] tdr;
        for (int i=0; i<4; ++i)
            for (int j=0; j<4; ++j)
                if (fmm_coeff[i][j]) delete[] fmm_coeff[i][j];
    }
    
    void PolarMesh2D::SetStepSize(const int m, const int n,
                                  const int sd, const int si,
                                  const double deltar,
                                  const double r0,
                                  PM2_Boundary& thebc)
    {
        maxi = m;
        maxj = n;
        bc = &thebc;
        bc->SetMesh(this);
        tmaxi = bc->Width(BDRY2D_XLO)+bc->Width(BDRY2D_XHI)+maxi;
        tmaxj = bc->Width(BDRY2D_YLO)+bc->Width(BDRY2D_YHI)+maxj;
        BuildWorkGrid(sd);
        theidata = new int[tmaxi*tmaxj*si];
        memset(theidata, 0, tmaxi*tmaxj*si*sizeof(int));
        tr = new double[tmaxi];
        tdr = new double[tmaxi];
        //original version
//		  for (int i=0; i<tmaxi; ++i) {
//            tr[i] = r0+(i-1)*deltar;
//            tdr[i] = deltar;
//        }
        //new
		  for (int i=0; i<tmaxi; ++i) {
            tr[i] = r0+(i-bc->Width(BDRY2D_XLO))*deltar;
            tdr[i] = deltar;
        }
        r = &tr[bc->Width(BDRY2D_XLO)];
        dr = &tdr[bc->Width(BDRY2D_XLO)];
        dtheta = 2*M_PI/n;
        InitFMMCoeffs();
    }
    
    void PolarMesh2D::SetBoundsSize(const int m, const int n,
                                      const int sd, const int si,
                                      const double r0, const double r1,
                                      PM2_Boundary& thebc)
    {
        maxi = m;
        maxj = n;
        bc = &thebc;
        bc->SetMesh(this);
        tmaxi = bc->Width(BDRY2D_XLO)+bc->Width(BDRY2D_XHI)+maxi;
        tmaxj = bc->Width(BDRY2D_YLO)+bc->Width(BDRY2D_YHI)+maxj;
        BuildWorkGrid(sd);
        theidata = new int[tmaxi*tmaxj*si];
        memset(theidata, 0, tmaxi*tmaxj*si*sizeof(int));
        tr = new double[tmaxi];
        tdr = new double[tmaxi];
        double deltar = (r1-r0)/(m-1);
		  //original
//		  for(int i=0; i<tmaxi; ++i) {
//			   tr[i] = r0+(i-1)*deltar;
//				tdr[i] = deltar;
//		  } 
		  //new
        for (int i=0; i<tmaxi; ++i) {
            tr[i] = r0+(i-bc->Width(BDRY2D_XLO))*deltar;
            tdr[i] = deltar;
        }
        r = &tr[bc->Width(BDRY2D_XLO)];
        dr = &tdr[bc->Width(BDRY2D_XLO)];
        dtheta = 2*M_PI/n;
        InitFMMCoeffs();
    }
    
    void PolarMesh2D::InitFMMCoeffs(void)
    {
        // the fmm coeffs are set up for all cases, first index corresponds to the case #
        // case 0: (r cos(theta), r sin(theta)), ((r-dr) cos(theta), (r-dr) sin(theta)),
        //                                  ((r-dr) cos(theta+dtheta), (r-dr) sin(theta+dtheta))
        // case 1: (r cos(theta), r sin(theta)), ((r-dr) cos(theta+dtheta), (r-dr) sin(theta+dtheta)),
        //                                  (r cos(theta+dtheta), r sin(theta+dtheta))
        // case 2: (r cos(theta), r sin(theta)), ((r-dr) cos(theta), (r-dr) sin(theta)),
        //                                  (r cos(theta+dtheta), r sin(theta+dtheta))
        // case 3: (r cos(theta), r sin(theta)), (r cos(theta+dtheta), r sin(theta+dtheta)),
        //                                  ((r+dr) cos(theta), (r+dr) sin(theta))
        
        // fmm_coeff are the coefficients involved in computing the new u value, it has the form:
        // fmm_coeff[i][0] * u1 + fmm_coeff[i][1] * u2
        //                  + fmm_coeff[i][2] * sqrt(fmm_coeff[i][3]*F^2 + (u1-u2)^2)
        
        // fmm_val determines whether a computed u value should be accepted.
        // There are two conditions:
        // u >= fmm_val[i][0] * u1 + fmm_val[i][1] * u2
        // u >= fmm_val[i][2] * u1 + fmm_val[i][3] * u2
        
        // fmm_ext gives the velocity extension weights. There are overlaps, so watch indices
        // G = (fmm_ext[i][1] * (u-u1) + fmm_ext[i][2] * (u-u2)) * G1
        //     / ((fmm_ext[i][1]+fmm_ext[i][2]) * (u-u1) + (fmm_ext[i][0]+fmm_ext[i][2]) * (u-u2))
        //     + (fmm_ext[i][2] * (u-u1) + fmm_ext[i][0] * (u-u2)) * G2
        //     / ((fmm_ext[i][1]+fmm_ext[i][2]) * (u-u1) + (fmm_ext[i][0]+fmm_ext[i][2]) * (u-u2))
        
        for (int i=0; i<4; ++i) {
            for (int j=0; j<4; ++j) {
                fmm_coeff[i][j] = new double[tmaxi]; // u1, u2, outside sqrt, multiplying F^2
                fmm_val[i][j] = new double[tmaxi];
            }
            for (int j=0; j<3; ++j)
                fmm_ext[i][j] = new double[tmaxi];
        }
        // the fmm_coi is true if the fmm update must include the diagonal
        // in that case, use fmm_coeff[0] and fmm_coeff[1] for looking inward
        // otherwise use fmm_coeff[2]. Use fmm_coeff[3] for looking outward

		  //old (doesn't default to false?)
		  //fmm_coi = new bool[tmaxi];
		  //new
		  fmm_coi = new bool[tmaxi]();
        
        for (int i=0; i<maxi; ++i) {
            if(i==0) {
				 fmm_coi[i] = 0;

				} else {
					fmm_coi[i] = dr[i-1] < 2.*r[i]*(1.-cos(dtheta));
				}
				
				double x[3][2], d[2], P[2][2], PP, denom;
            // case 0:
            x[0][0] = r[i];                       x[0][1] = 0;
            x[1][0] = r[i]-dr[i-1];               x[1][1] = 0;
            x[2][0] = (r[i]-dr[i-1])*cos(dtheta); x[2][1] = (r[i]-dr[i-1])*sin(dtheta);
            d[0] = sqrt((x[0][0]-x[1][0])*(x[0][0]-x[1][0])+(x[0][1]-x[1][1])*(x[0][1]-x[1][1]));
            d[1] = sqrt((x[0][0]-x[2][0])*(x[0][0]-x[2][0])+(x[0][1]-x[2][1])*(x[0][1]-x[2][1]));
            P[0][0] = (x[0][0]-x[1][0])/d[0]; P[0][1] = (x[0][1]-x[1][1])/d[0];
            P[1][0] = (x[0][0]-x[2][0])/d[1]; P[1][1] = (x[0][1]-x[2][1])/d[1];
            PP = P[0][0]*P[1][0]+P[0][1]*P[1][1];
            denom = (d[0]*P[0][0]-d[1]*P[1][0])*(d[0]*P[0][0]-d[1]*P[1][0])
                        +(d[0]*P[0][1]-d[1]*P[1][1])*(d[0]*P[0][1]-d[1]*P[1][1]);
            fmm_val[0][0][i] = d[1]/(d[1]-d[0]*PP);
            fmm_val[0][1][i] = d[0]*PP/(d[0]*PP-d[1]);
            fmm_val[0][2][i] = d[1]*PP/(d[1]*PP-d[0]);
            fmm_val[0][3][i] = d[0]/(d[0]-d[1]*PP);
            fmm_coeff[0][0][i] = d[1]*(P[1][0]*(d[1]*P[1][0]-d[0]*P[0][0])+P[1][1]*(d[1]*P[1][1]-d[0]*P[0][1]))/denom;
            fmm_coeff[0][1][i] = -d[0]*(P[0][0]*(d[1]*P[1][0]-d[0]*P[0][0])+P[0][1]*(d[1]*P[1][1]-d[0]*P[0][1]))/denom;
            fmm_coeff[0][2][i] = -d[0]*d[1]*(P[0][0]*P[1][1]-P[0][1]*P[1][0])/denom;
            fmm_coeff[0][3][i] = denom;
            fmm_ext[0][0][i] = d[0]*d[0];
            fmm_ext[0][1][i] = d[1]*d[1];
            fmm_ext[0][2][i] = -d[0]*d[1]*(P[0][0]*P[1][0]+P[0][1]*P[1][1]);

            // case 1:
            x[0][0] = r[i];                       x[0][1] = 0;
            x[1][0] = (r[i]-dr[i-1])*cos(dtheta); x[1][1] = (r[i]-dr[i-1])*sin(dtheta);
            x[2][0] = r[i]*cos(dtheta);           x[2][1] = r[i]*sin(dtheta);
            d[0] = sqrt((x[0][0]-x[1][0])*(x[0][0]-x[1][0])+(x[0][1]-x[1][1])*(x[0][1]-x[1][1]));
            d[1] = sqrt((x[0][0]-x[2][0])*(x[0][0]-x[2][0])+(x[0][1]-x[2][1])*(x[0][1]-x[2][1]));
            P[0][0] = (x[0][0]-x[1][0])/d[0]; P[0][1] = (x[0][1]-x[1][1])/d[0];
            P[1][0] = (x[0][0]-x[2][0])/d[1]; P[1][1] = (x[0][1]-x[2][1])/d[1];
            PP = P[0][0]*P[1][0]+P[0][1]*P[1][1];
            denom = (d[0]*P[0][0]-d[1]*P[1][0])*(d[0]*P[0][0]-d[1]*P[1][0])
            +(d[0]*P[0][1]-d[1]*P[1][1])*(d[0]*P[0][1]-d[1]*P[1][1]);
            fmm_val[1][0][i] = d[1]/(d[1]-d[0]*PP);
            fmm_val[1][1][i] = d[0]*PP/(d[0]*PP-d[1]);
            fmm_val[1][2][i] = d[1]*PP/(d[1]*PP-d[0]);
            fmm_val[1][3][i] = d[0]/(d[0]-d[1]*PP);
            fmm_coeff[1][0][i] = d[1]*(P[1][0]*(d[1]*P[1][0]-d[0]*P[0][0])+P[1][1]*(d[1]*P[1][1]-d[0]*P[0][1]))/denom;
            fmm_coeff[1][1][i] = -d[0]*(P[0][0]*(d[1]*P[1][0]-d[0]*P[0][0])+P[0][1]*(d[1]*P[1][1]-d[0]*P[0][1]))/denom;
            fmm_coeff[1][2][i] = -d[0]*d[1]*(P[0][0]*P[1][1]-P[0][1]*P[1][0])/denom;
            fmm_coeff[1][3][i] = denom;
            fmm_ext[1][0][i] = d[0]*d[0];
            fmm_ext[1][1][i] = d[1]*d[1];
            fmm_ext[1][2][i] = -d[0]*d[1]*(P[0][0]*P[1][0]+P[0][1]*P[1][1]);

            // case 2:
            x[0][0] = r[i];                       x[0][1] = 0;
            x[1][0] = r[i]-dr[i-1];               x[1][1] = 0;
            x[2][0] = r[i]*cos(dtheta);           x[2][1] = r[i]*sin(dtheta);
            d[0] = sqrt((x[0][0]-x[1][0])*(x[0][0]-x[1][0])+(x[0][1]-x[1][1])*(x[0][1]-x[1][1]));
            d[1] = sqrt((x[0][0]-x[2][0])*(x[0][0]-x[2][0])+(x[0][1]-x[2][1])*(x[0][1]-x[2][1]));
            P[0][0] = (x[0][0]-x[1][0])/d[0]; P[0][1] = (x[0][1]-x[1][1])/d[0];
            P[1][0] = (x[0][0]-x[2][0])/d[1]; P[1][1] = (x[0][1]-x[2][1])/d[1];
            PP = P[0][0]*P[1][0]+P[0][1]*P[1][1];
            denom = (d[0]*P[0][0]-d[1]*P[1][0])*(d[0]*P[0][0]-d[1]*P[1][0])
            +(d[0]*P[0][1]-d[1]*P[1][1])*(d[0]*P[0][1]-d[1]*P[1][1]);
            fmm_val[2][0][i] = d[1]/(d[1]-d[0]*PP);
            fmm_val[2][1][i] = d[0]*PP/(d[0]*PP-d[1]);
            fmm_val[2][2][i] = d[1]*PP/(d[1]*PP-d[0]);
            fmm_val[2][3][i] = d[0]/(d[0]-d[1]*PP);
            fmm_coeff[2][0][i] = d[1]*(P[1][0]*(d[1]*P[1][0]-d[0]*P[0][0])+P[1][1]*(d[1]*P[1][1]-d[0]*P[0][1]))/denom;
            fmm_coeff[2][1][i] = -d[0]*(P[0][0]*(d[1]*P[1][0]-d[0]*P[0][0])+P[0][1]*(d[1]*P[1][1]-d[0]*P[0][1]))/denom;
            fmm_coeff[2][2][i] = -d[0]*d[1]*(P[0][0]*P[1][1]-P[0][1]*P[1][0])/denom;
            fmm_coeff[2][3][i] = denom;
            fmm_ext[2][0][i] = d[0]*d[0];
            fmm_ext[2][1][i] = d[1]*d[1];
            fmm_ext[2][2][i] = -d[0]*d[1]*(P[0][0]*P[1][0]+P[0][1]*P[1][1]);

            // case 3:
            x[0][0] = r[i];                       x[0][1] = 0;
            x[1][0] = r[i]*cos(dtheta);           x[1][1] = r[i]*sin(dtheta);
            x[2][0] = r[i]+dr[i];                 x[2][1] = 0;
            d[0] = sqrt((x[0][0]-x[1][0])*(x[0][0]-x[1][0])+(x[0][1]-x[1][1])*(x[0][1]-x[1][1]));
            d[1] = sqrt((x[0][0]-x[2][0])*(x[0][0]-x[2][0])+(x[0][1]-x[2][1])*(x[0][1]-x[2][1]));
            P[0][0] = (x[0][0]-x[1][0])/d[0]; P[0][1] = (x[0][1]-x[1][1])/d[0];
            P[1][0] = (x[0][0]-x[2][0])/d[1]; P[1][1] = (x[0][1]-x[2][1])/d[1];
            PP = P[0][0]*P[1][0]+P[0][1]*P[1][1];
            denom = (d[0]*P[0][0]-d[1]*P[1][0])*(d[0]*P[0][0]-d[1]*P[1][0])
            +(d[0]*P[0][1]-d[1]*P[1][1])*(d[0]*P[0][1]-d[1]*P[1][1]);
            fmm_val[3][0][i] = d[1]/(d[1]-d[0]*PP);
            fmm_val[3][1][i] = d[0]*PP/(d[0]*PP-d[1]);
            fmm_val[3][2][i] = d[1]*PP/(d[1]*PP-d[0]);
            fmm_val[3][3][i] = d[0]/(d[0]-d[1]*PP);
            fmm_coeff[3][0][i] = d[1]*(P[1][0]*(d[1]*P[1][0]-d[0]*P[0][0])+P[1][1]*(d[1]*P[1][1]-d[0]*P[0][1]))/denom;
            fmm_coeff[3][1][i] = -d[0]*(P[0][0]*(d[1]*P[1][0]-d[0]*P[0][0])+P[0][1]*(d[1]*P[1][1]-d[0]*P[0][1]))/denom;
            fmm_coeff[3][2][i] = -d[0]*d[1]*(P[0][0]*P[1][1]-P[0][1]*P[1][0])/denom;
            fmm_coeff[3][3][i] = denom;
            fmm_ext[3][0][i] = d[0]*d[0];
            fmm_ext[3][1][i] = d[1]*d[1];
            fmm_ext[3][2][i] = -d[0]*d[1]*(P[0][0]*P[1][0]+P[0][1]*P[1][1]);
        }
    }

    void PolarMesh2D::BuildWorkGrid(int num)
    {
        thedata = new double[tmaxi*tmaxj*num];
        datastart = Index(bc->Width(BDRY2D_XLO), bc->Width(BDRY2D_YLO));
        maxk = num;
        memset(thedata, 0, tmaxi*tmaxj*num*sizeof(double));
    }
    
	 void PolarMesh2D::CopyWorkGrid(int to, int from)
	 {
		  double* datato = &(thedata[Index(0,0,to)]);
		  double* datafrom = &(thedata[Index(0,0,from)]);
		  memcpy(datato, datafrom, tmaxi*tmaxj*sizeof(double));
	 }


    double PolarMesh2D::X(const int i, const int j) const
    {
        return r[i]*cos(j*dtheta);
    }
    
    double PolarMesh2D::Y(const int i, const int j) const
    {
        return r[i]*sin(j*dtheta);
    }
    
    double PolarMesh2D::X(const double i, const double j) const
    {
        return (r[int(i)]+(i-floor(i))*dr[int(i)])*cos(j*dtheta);
    }
    
    double PolarMesh2D::Y(const double i, const double j) const
    {
        return (r[int(i)]+(i-floor(i))*dr[int(i)])*sin(j*dtheta);
    }

    int PolarMesh2D::I(const double x, const double y) const
    {
        double rho=sqrt(x*x+y*y);
        int i;
        for (i=maxi-1+bc->Width(BDRY2D_XHI); i>=0 && r[i]>rho; --i) ;
        return i-bc->Width(BDRY2D_XLO);
    }
    
    double PolarMesh2D::Id(const double x, const double y) const
    {
        double rho=sqrt(x*x+y*y);
        int i;  
		  for (i=maxi-1+bc->Width(BDRY2D_XHI); i>=0 && r[i]>rho; --i) ;
        //return (i-bc->Width(BDRY2D_XLO))+(rho-r[i])/dr[i];
    	  return i+(rho-r[i])/dr[i];
	 }
    
    int PolarMesh2D::I(const double rho) const
    {
        int i;
        for (i=maxi-1+bc->Width(BDRY2D_XHI); i>=0 && r[i]>rho; --i) ;
        return i-bc->Width(BDRY2D_XLO);
    }
    
    double PolarMesh2D::Id(const double rho) const
    {
        int i;
        for (i=maxi-1+bc->Width(BDRY2D_XHI); i>=0 && r[i]>rho; --i) ;
        //return i-bc->Width(BDRY2D_XLO)+(rho-r[i])/dr[i];
    	  return i+(rho-r[i])/dr[i];
	 }
    
    void PolarMesh2D::LocToIndexXY(const double x, const double y, int& i, int& j,
                                   double& ifrac, double& jfrac, const char round) const
    {
        double id = Id(x, y);
        double jd = Jd(x, y);
//		  std::cout << "(id,jd) = " << id << ", " << jd << std::endl;
        if (round) {
            i = int(id+0.5);
            j = int(jd+0.5);
            ifrac = id - i;
            jfrac = jd - j;
        } else {
//            std::cout << "We did not round\n";
				i = int(id);
            j = int(jd);
            ifrac = id - i;
            jfrac = jd - j;
        }
    }

    void PolarMesh2D::LocToIndexRT(const double rho, const double theta, int& i, int& j,
                                   double& ifrac, double& jfrac, const char round) const
    {
        double id = Id(rho);
        double jd = Jd(theta);
        if (round) {
            i = int(id+0.5);
            j = int(jd+0.5);
            ifrac = id - i;
            jfrac = jd - j;
        } else {
            i = int(id);
            j = int(jd);
            ifrac = id - i;
            jfrac = jd - j;
        }
    }
    
    double PolarMesh2D::InterpXY(const double& x, const double& y, const int k) const
    {
        int i, j;
        double ifrac, jfrac;
        LocToIndexXY(x, y, i, j, ifrac, jfrac);
//		  std::cout << "(x,y,i,j,if,jf) = " << x << ", " << y << ", " << i << ", " << j << ", " << ifrac << ", " << jfrac << std::endl;
//		  std::cout << "phi = " << data_(i,j,k) << std::endl;
//		  std::cout << "phiE = " << data_(i+1,j,k) << std::endl;
//		  std::cout << "phiN = " << data_(i,j+1,k) << std::endl;
//		  std::cout << "phiNE = " << data_(i+1,j+1,k) << std::endl;
		  return data_(i,j,k)*(1.-ifrac)*(1.-jfrac) + data_(i+1,j,k)*ifrac*(1-jfrac)
            +data_(i,j+1,k)*(1.-ifrac)*jfrac + data_(i+1,j+1,k)*ifrac*jfrac;
    }

    double PolarMesh2D::InterpRT(const double& rho, const double& theta, const int k) const
    {
        int i, j;
        double ifrac, jfrac;
        LocToIndexRT(rho, theta, i, j, ifrac, jfrac);
        return data_(i,j,k)*(1.-ifrac)*(1.-jfrac) + data_(i+1,j,k)*ifrac*(1-jfrac)
        +data_(i,j+1,k)*(1.-ifrac)*jfrac + data_(i+1,j+1,k)*ifrac*jfrac;
    }
    
    double PolarMesh2D::Interp2XY(const double& x, const double& y, const int k) const
    {
        int i, j;
        double ifrac, jfrac;
        LocToIndexXY(x, y, i, j, ifrac, jfrac);
        return Interp2RT(i, j, k, ifrac, jfrac);
    }
    
    double PolarMesh2D::Interp2RT(const double& rho, const double& theta, const int k) const
    {
        int i, j;
        double ifrac, jfrac;
        LocToIndexRT(rho, theta, i, j, ifrac, jfrac);
        return Interp2RT(i, j, k, ifrac, jfrac);
    }
    
    double PolarMesh2D::Interp2RT(const int i, const int j, const int k, const double& ifrac,
                     const double& jfrac) const
    {
        GeneralBicubic p = Interp2(i, j, k);
        return p(ifrac*dr[i-1], jfrac*dtheta);
    }
    
    GeneralBicubic PolarMesh2D::Interp2(const int i, const int j, const int k) const
    {
        GeneralBicubic p;
        p.BuildwDeriv(data_(i-1,j-1,k),data_(i,j-1,k),data_(i+1,j-1,k),
                      data_(i+2,j-1,k),data_(i-1,j,k),data_(i,j,k),data_(i+1,j,k),
                      data_(i+2,j,k),data_(i-1,j+1,k),data_(i,j+1,k),
                      data_(i+1,j+1,k),data_(i+2,j+1,k),data_(i-1,j+2,k),
                      data_(i,j+2,k),data_(i+1,j+2,k),data_(i+2,j+2,k),
                      dr[i-1],dr[i],dr[i+1],dtheta,dtheta,dtheta);
        return p;
    }
}
