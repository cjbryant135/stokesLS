//
//  genbicubic.hpp
//  LevelSet
//
//  Created by David Chopp on 8/8/19.
//  Copyright © 2019 David Chopp. All rights reserved.
//

#ifndef genbicubic_hpp
#define genbicubic_hpp

namespace levelset {
    
    class GeneralBicubic
    {
    public:
        
        double a[4][4];
        double wx[3], wy[3];
        
        GeneralBicubic(const GeneralBicubic& p);
        
        GeneralBicubic(const double f11, const double f21, const double f12,
                const double f22, const double fx11, const double fx21,
                const double fx12, const double fx22, const double fy11,
                const double fy21, const double fy12, const double fy22,
                const double fxy11, const double fxy21, const double fxy12,
                const double fxy22, const double dx[3], const double dy[3])
        {
            for (int i=0; i<3; ++i) {
                wx[i] = dx[i]; wy[i] = dy[i];
            }
            Build(f11,f21,f12,f22,fx11,fx21,fx12,fx22,fy11,fy21,fy12,fy22,
                  fxy11,fxy21,fxy12,fxy22,dx[0], dx[1], dx[2], dy[0], dy[1], dy[2]);
        }
        
        GeneralBicubic(const double f11, const double f21, const double f12,
                       const double f22, const double fx11, const double fx21,
                       const double fx12, const double fx22, const double fy11,
                       const double fy21, const double fy12, const double fy22,
                       const double fxy11, const double fxy21, const double fxy12,
                       const double fxy22, const double dxm, const double dx0,
                       const double dxp, const double dym, const double dy0,
                       const double dyp)
        {
            wx[0] = dxm; wx[1] = dx0; wx[2] = dxp;
            wy[0] = dym; wy[1] = dy0; wy[2] = dyp;
            Build(f11,f21,f12,f22,fx11,fx21,fx12,fx22,fy11,fy21,fy12,fy22,
                  fxy11,fxy21,fxy12,fxy22,dxm,dx0,dxp,dym,dy0,dyp);
        }
        
        GeneralBicubic(void) {
            for (int i=0; i<4; ++i)
                for (int j=0; j<4; ++j) a[i][j] = 0.;
        }
        
        void Build(const double f11, const double f21, const double f12,
                   const double f22, const double fx11, const double fx21,
                   const double fx12, const double fx22, const double fy11,
                   const double fy21, const double fy12, const double fy22,
                   const double fxy11, const double fxy21, const double fxy12,
                   const double fxy22, const double dx0, const double dx1,
                   const double dx2, const double dy0, const double dy1,
                   const double dy2);
        
        void BuildwDeriv(const double f00, const double f10, const double f20,
                         const double f30, const double f01, const double f11,
                         const double f21, const double f31, const double f02,
                         const double f12, const double f22, const double f32,
                         const double f03, const double f13, const double f23,
                         const double f33, const double dx0, const double dx1,
                         const double dx2, const double dy0, const double dy1,
                         const double dy2);
               
        double operator()(const double x, const double y) const;
        double F(const double x, const double y) const
        {return operator()(x,y);}
        double Dx(const double x, const double y) const;
        double Dy(const double x, const double y) const;
        double Dxx(const double x, const double y) const;
        double Dxy(const double x, const double y) const;
        double Dyy(const double x, const double y) const;
        double Dxxx(const double x, const double y) const;
        double Dxxy(const double x, const double y) const;
        double Dxyy(const double x, const double y) const;
        double Dyyy(const double x, const double y) const;
        
        GeneralBicubic Dx(void);
        GeneralBicubic Dy(void);
        
        double LocalDist(const double x, const double y,
                         double& ax, double& ay, char& clean) const;
        double LocalDistPolar(const double r0, const double r, const double theta,
                         double& ar, double& atheta, char& clean) const;
        double XCrossing(const double x, const double y, double& ay, char& clean) const;
        double YCrossing(const double x, const double y, double& ax, char& clean) const;
    };
    
}
#endif /* genbicubic_hpp */

