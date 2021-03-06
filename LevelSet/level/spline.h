/*************************************************
    spline.h

    $Header: spline.h,v 1.1 99/02/04 14:38:44 chopp Exp $

    $Log:       spline.h,v $
    * Revision 1.1  99/02/04  14:38:44  14:38:44  chopp (David Chopp)
    * Initial revision
    * 
    *************************************************/

#ifndef __SPLINE_H__
#define __SPLINE_H__

#include <stdlib.h>

namespace levelset {
        
    class Spline {
    protected:

        int num;
        const int degree;
        double*** coeff;
        double* length;
        
        void Build(const int num = 0, const double* x_nodes = NULL, 
                   const double* y_nodes = NULL);

    public:

        Spline(const int n = 0, const double* x_nodes = NULL, 
               const double* y_nodes = NULL);
        
        ~Spline(void);
        
        void ChangeNodes(const int n, const double* x_nodes, 
                         const double* y_nodes);

        double Length(void) const;
        
        double x(const double s) const;
        double y(const double s) const;
        double dx(const double s) const;
        double dy(const double s) const;
        double d2x(const double s) const;
        double d2y(const double s) const;
        double d3x(const double s) const;
        double d3y(const double s) const;
        double d4x(const double s) const;
        double d4y(const double s) const;
        
        double x(const int i, const double s) const;
        double y(const int i, const double s) const;
        double dx(const int i, const double s) const;
        double dy(const int i, const double s) const;
        double d2x(const int i, const double s) const;
        double d2y(const int i, const double s) const;
        double d3x(const int i, const double s) const;
        double d3y(const int i, const double s) const;
        double d4x(const int i, const double s) const;
        double d4y(const int i, const double s) const;

        double Curvature(const int i, const double s) const;
        double Curvature(const double s) const;
        double Diffusion(const int i, const double s) const;
        double Diffusion(const double s) const;
    };

}
#endif // __SPLINE_H__
