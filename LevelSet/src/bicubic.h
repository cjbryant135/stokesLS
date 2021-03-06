/*************************************************
 bicubic.h
 
 $Header: bicubic.h,v 1.1 2000/05/31 11:08:23 chopp Exp $
 
 $Log:  bicubic.h,v $
 * Revision 1.1  2000/05/31  11:08:23  11:08:23  chopp (David Chopp)
 * Initial revision
 * 
 *************************************************/

#ifndef __BICUBIC_H__
#define __BICUBIC_H__

namespace levelset {
        
    class Bicubic 
    {
    public:
                        
        double a[4][4];
        double wx, wy;
                        
        Bicubic(const Bicubic& p);
                        
    Bicubic(const double f11, const double f21, const double f12,
            const double f22, const double fx11, const double fx21,
            const double fx12, const double fx22, const double fy11,
            const double fy21, const double fy12, const double fy22,
            const double fxy11, const double fxy21, const double fxy12,
            const double fxy22, const double dx, const double dy) :
        wx(dx), wy(dy)
        {Build(f11,f21,f12,f22,fx11,fx21,fx12,fx22,fy11,fy21,fy12,fy22,
               fxy11,fxy21,fxy12,fxy22,dx,dy);}
                        
    Bicubic(void) : wx(0.), wy(0.){
            for (int i=0; i<4; ++i)
                for (int j=0; j<4; ++j) a[i][j] = 0.;
        }
                        
        void Build(const double f11, const double f21, const double f12,
                   const double f22, const double fx11, const double fx21,
                   const double fx12, const double fx22, const double fy11,
                   const double fy21, const double fy12, const double fy22,
                   const double fxy11, const double fxy21, const double fxy12,
                   const double fxy22, const double dx, const double dy);
                        
        void BuildwDeriv(const double f00, const double f10, const double f20,
                         const double f30, const double f01, const double f11,
                         const double f21, const double f31, const double f02,
                         const double f12, const double f22, const double f32,
                         const double f03, const double f13, const double f23,
                         const double f33, const double dx, const double dy);
                        
        void BuildwDeriv(const double f00, const double f10, const double f20,
                         const double f30, const double f01, const double f11,
                         const double f21, const double f31, const double f02,
                         const double f12, const double f22, const double f32,
                         const double f03, const double f13, const double f23,
                         const double f33, const double dx, const double dy, const char* valid);
                        
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
                        
        Bicubic Dx(void);
        Bicubic Dy(void);
                        
        double LocalDist(const double x, const double y, 
                         double& ax, double& ay, char& clean) const;
        double LocalDistNewton(const double x, const double y, 
                               double& ax, double& ay, char& clean) const;
        double LocalDistByBoxes(const double x, const double y,
                                double& ax, double& ay) const;
        double XCrossing(const double x, const double y, double& ay, char& clean) const;
        double YCrossing(const double x, const double y, double& ax, char& clean) const;

        void FollowToBdry(const int dir, const double x, const double y, double& ax, double& ay) const;

        double Integral(void);
    };
        
}

#endif
