#ifndef __MESH2D_H__
#define __MESH2D_H__

#include "compiler.h"
#include "initialfunc.h"
#include "interface.h"

namespace levelset {
        
/// List of neighboring gridpoints by deltai and deltaj
    const int     nabor[2][8] = {1, 1, 0, -1, -1, -1,  0,  1, 
                                 0, 1, 1,  1,  0, -1, -1, -1};

/// Abstract 2-dimensional mesh class, may be multi-valued
    class Mesh2D {

    public:

        /// Default Constructor
        Mesh2D(void) {;} 

        /** @name Initialization Routines */
        //@{
        /// Initialize the values of the mesh with function f
        virtual void SetValues(double (*f)(double,double), const int wi) =0;
        //@}
   
        /** @name Computational Routines */
        //@{
        /** @name Interface Information */
        //@{
        /** Retrieve a list of interface points. This routine must be
         * called before any of the below routines can be used. The
         * values returned by the other computational routines correspond
         * to the interface points computed in Interface
         */
        virtual void GetInterface(Interface& s, const int wi, const double value=0) =0;
        /// Retrieve list of unit normal vectors at interface points
        virtual void UnitNormal(Interface& s,
                                const int nx, const int ny, 
                                const int kdx, const int kdy) =0;
        /// Retrieve other data values at the interface points
        //virtual DataVec<double>& Values(const int slice) const =0;
        //@}
        /** @name Derivative Information 
         * The information returned corresponds to the values at each point
         * of the Interface. See Interface for more information
         * @see Mesh2D::Interface()
         */
        //@{
        /// Retrieve central difference x-derivative
        virtual void Dx_zero(Interface& s, const int n) =0;
        /// Retrieve left-sided x-derivative
        virtual void Dx_minus(Interface& s, const int n) =0;
        /// Retrieve right-sided x-derivative
        virtual void Dx_plus(Interface& s, const int n) =0;
        /// Retrieve upwind x-derivative
        virtual void Dx_upwind(Interface& s, const int n, const int d = 0) =0;
        /// Retrieve second x-derivative
        virtual void Dxx_zero(Interface& s, const int n) =0;
        /// Retrieve central difference y-derivative
        virtual void Dy_zero(Interface& s, const int n) =0;
        /// Retrieve down-sided y-derivative
        virtual void Dy_minus(Interface& s, const int n) =0;
        /// Retrieve top-sided y-derivative
        virtual void Dy_plus(Interface& s, const int n) =0;
        /// Retrieve upwind y-derivative
        virtual void Dy_upwind(Interface& s, const int n, const int d = 0) =0;
        /// Retrieve second y-derivative
        virtual void Dyy_zero(Interface& s, const int n) =0;
        //@}
        /** @name Motion Routines 
         * These routines may only be called after the Interface is located.
         * See Mesh2D::Interface for more information
         * @see Mesh2D::Interface()
         */
        //@{
        /// Advance the interface with normal velocity v and time step dt
        virtual void Advance(Interface& s, const double dt) =0;
        /// Advance the interface with velocity (vx, vy) and time step dt
//   void Advance(DataVec<double>& vx, DataVec<double>& vy, double dt) =0;
        //@}
        //@}

        /** @name Graphics Routines */
        //@{
#ifndef NO_GRAPHICS
        /// Plot the interface in the port
        virtual void Plot(PlotWindow2D& port, const double level_value = 0,
                          const char showboxes = 0x00) const =0;
#endif
        /// Plot the interface in text format to cout
        virtual void TPlot(const double level_value = 0,
                           const char showboxes = 0x00) const =0;
        //@}

    protected:
 
        char IsCrossing(const double v, const double y0, const double y1, double& f);
    };

}
#endif
