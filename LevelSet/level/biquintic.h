/*************************************************
    biquintic.h

    $Header: biquintic.h,v 1.1 2000/05/31 11:08:23 chopp Exp $

    $Log:       biquintic.h,v $
    * Revision 1.1  2000/05/31  11:08:23  11:08:23  chopp (David Chopp)
    * Initial revision
    * 
    *************************************************/

#ifndef __BIQUINTIC_H__
#define __BIQUINTIC_H__

namespace levelset {

    class Biquintic 
    {
    public:

        double a[6][6];
        double wx, wy;
   
        Biquintic(const Biquintic& p);
  
    Biquintic(const double* g, const double dx, const double dy) :
        wx(dx), wy(dy)
        {Build(g,dx,dy);}
   
    Biquintic(void) : wx(0.), wy(0.){
            for (int i=0; i<6; ++i)
                for (int j=0; j<6; ++j) a[i][j] = 0.;
        }
   
        void Build(const double* g, const double dx, const double dy);
   
        void BuildwDeriv(const double f[6][6], const double dx, const double dy);
   
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
        
        Biquintic Dx(void);
        Biquintic Dy(void);
    
        double LocalDist(const double x, const double y, 
                         double& ax, double& ay, char& clean) const;
    };
        
}

#endif
