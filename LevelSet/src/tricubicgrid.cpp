/*************************************************
 tricubicgrid.cpp
 
 Revision 1.0  2017/11/28  11:05:46  11:05:46  chopp (David Chopp)
 Initial revision
 
*************************************************/

#include "tricubicgrid.h"
//#define DEBUG_TRICUBICGRID
#ifdef DEBUG_TRICUBICGRID
#include "plotwindow2d.h"
#endif

namespace levelset {

    TricubicGrid::TricubicGrid(const UniformMesh3D* themesh, const int l, const bool xperiodic, const bool yperiodic,
                               const bool zperiodic) 
        : mesh(themesh), grid(NULL), maxi(themesh->maxi), maxj(themesh->maxj), maxk(themesh->maxk),
          dx(themesh->dx), dy(themesh->dy), dz(themesh->dz), kval(l)
    {
        periodic[0] = xperiodic;
        periodic[1] = yperiodic;
        periodic[2] = zperiodic;
        grid = new Tricubic*[maxi*maxj*maxk];
        for (int i=0; i<maxi*maxj*maxk; ++i)
            grid[i] = NULL;
        zero[0] = themesh->zero[0];
        zero[1] = themesh->zero[1];
        zero[2] = themesh->zero[2];
    }

    TricubicGrid::~TricubicGrid(void)
    {
        if (grid) {
            for (int i=0; i<maxi*maxj*maxk; ++i)
                if (grid[i]) delete grid[i];
            delete[] grid;
        }
    }

    Tricubic* TricubicGrid::GetTricubic(const int ii, const int jj, const int kk)
    {
        int i, j, k;
        if (periodic[0]) i=(ii+maxi)%maxi;
        else i = ii;
        if (periodic[1]) j=(jj+maxj)%maxj;
        else j = jj;
        if (periodic[2]) k=(kk+maxk)%maxk;
        else k = kk;
        Tricubic* p = grid[Index(i,j,k)];
        if (p == NULL) {
            p = grid[Index(i,j,k)] = new Tricubic;
            double data[4][4][4];
            if (periodic[0]) {
                if (periodic[1]) {
                    if (periodic[2]) {
                        for (int l=0; l<4; ++l)
                            for (int m=0; m<4; ++m)
                                for (int n=0; n<4; ++n) 
                                    data[l][m][n] = mesh->data_((i+l-1+maxi)%maxi, (j+m-1+maxj)%maxj, 
                                                                (k+n-1+maxk)%maxk, kval);
                    } else {
                        for (int l=0; l<4; ++l)
                            for (int m=0; m<4; ++m)
                                for (int n=0; n<4; ++n) 
                                    data[l][m][n] = mesh->data_((i+l-1+maxi)%maxi, (j+m-1+maxj)%maxj, 
                                                                k+n-1, kval);
                    }
                } else {
                    if (periodic[2]) {
                        for (int l=0; l<4; ++l)
                            for (int m=0; m<4; ++m)
                                for (int n=0; n<4; ++n) 
                                    data[l][m][n] = mesh->data_((i+l-1+maxi)%maxi, j+m-1, 
                                                                (k+n-1+maxk)%maxk, kval);
                    } else {
                        for (int l=0; l<4; ++l)
                            for (int m=0; m<4; ++m)
                                for (int n=0; n<4; ++n) 
                                    data[l][m][n] = mesh->data_((i+l-1+maxi)%maxi, j+m-1, 
                                                                k+n-1, kval);
                    }
                }
            } else {
                if (periodic[1]) {
                    if (periodic[2]) {
                        for (int l=0; l<4; ++l)
                            for (int m=0; m<4; ++m)
                                for (int n=0; n<4; ++n) 
                                    data[l][m][n] = mesh->data_(i+l-1, (j+m-1+maxj)%maxj, 
                                                                (k+n-1+maxk)%maxk, kval);
                    } else {
                        for (int l=0; l<4; ++l)
                            for (int m=0; m<4; ++m)
                                for (int n=0; n<4; ++n) 
                                    data[l][m][n] = mesh->data_(i+l-1, (j+m-1+maxj)%maxj, 
                                                                k+n-1, kval);
                    }
                } else {
                    if (periodic[2]) {
                        for (int l=0; l<4; ++l)
                            for (int m=0; m<4; ++m)
                                for (int n=0; n<4; ++n) 
                                    data[l][m][n] = mesh->data_(i+l-1, j+m-1, 
                                                                (k+n-1+maxk)%maxk, kval);
                    } else {
                        for (int l=0; l<4; ++l)
                            for (int m=0; m<4; ++m)
                                for (int n=0; n<4; ++n) 
                                    data[l][m][n] = mesh->data_(i+l-1, j+m-1,  k+n-1, kval);
                    }
                }
            }

            p->BuildwDeriv(data, dx, dy, dz);
        }
        return p;
    }

    Tricubic* TricubicGrid::GetTricubic(const double x, const double y, const double z)
    {
        return GetTricubic(I(x), J(y), K(z));
    }

    double TricubicGrid::operator()(const double x, const double y, const double z) 
    {
        Tricubic* p = GetTricubic(x,y,z);
        return p->F(x-I(x)*dx,y-J(y)*dy,z-K(z)*dz);
    }
    
    double TricubicGrid::Dx(const double x, const double y, const double z) 
    {
        Tricubic* p = GetTricubic(x,y,z);
        return p->Dx(x-I(x)*dx,y-J(y)*dy,z-K(z)*dz);
    }
    
    double TricubicGrid::Dy(const double x, const double y, const double z) 
    {
        Tricubic* p = GetTricubic(x,y,z);
        return p->Dy(x-I(x)*dx,y-J(y)*dy,z-K(z)*dz);
    }
    
    double TricubicGrid::Dz(const double x, const double y, const double z) 
    {
        Tricubic* p = GetTricubic(x,y,z);
        return p->Dz(x-I(x)*dx,y-J(y)*dy,z-K(z)*dz);
    }
    
    double TricubicGrid::Dxx(const double x, const double y, const double z) 
    {
        Tricubic* p = GetTricubic(x,y,z);
        return p->Dxx(x-I(x)*dx,y-J(y)*dy,z-K(z)*dz);
    }
    
    double TricubicGrid::Dxy(const double x, const double y, const double z) 
    {
        Tricubic* p = GetTricubic(x,y,z);
        return p->Dxy(x-I(x)*dx,y-J(y)*dy,z-K(z)*dz);
    }
    
    double TricubicGrid::Dxz(const double x, const double y, const double z) 
    {
        Tricubic* p = GetTricubic(x,y,z);
        return p->Dxz(x-I(x)*dx,y-J(y)*dy,z-K(z)*dz);
    }
    
    double TricubicGrid::Dyy(const double x, const double y, const double z) 
    {
        Tricubic* p = GetTricubic(x,y,z);
        return p->Dyy(x-I(x)*dx,y-J(y)*dy,z-K(z)*dz);
    }
    
    double TricubicGrid::Dyz(const double x, const double y, const double z) 
    {
        Tricubic* p = GetTricubic(x,y,z);
        return p->Dyz(x-I(x)*dx,y-J(y)*dy,z-K(z)*dz);
    }
    
    double TricubicGrid::Dzz(const double x, const double y, const double z) 
    {
        Tricubic* p = GetTricubic(x,y,z);
        return p->Dzz(x-I(x)*dx,y-J(y)*dy,z-K(z)*dz);
    }
    
    double TricubicGrid::Dxxx(const double x, const double y, const double z) 
    {
        Tricubic* p = GetTricubic(x,y,z);
        return p->Dxxx(x-I(x)*dx,y-J(y)*dy,z-K(z)*dz);
    }
    
    double TricubicGrid::Dxxy(const double x, const double y, const double z) 
    {
        Tricubic* p = GetTricubic(x,y,z);
        return p->Dxxy(x-I(x)*dx,y-J(y)*dy,z-K(z)*dz);
    }
    
    double TricubicGrid::Dxxz(const double x, const double y, const double z) 
    {
        Tricubic* p = GetTricubic(x,y,z);
        return p->Dxxz(x-I(x)*dx,y-J(y)*dy,z-K(z)*dz);
    }
    
    double TricubicGrid::Dxyy(const double x, const double y, const double z) 
    {
        Tricubic* p = GetTricubic(x,y,z);
        return p->Dxyy(x-I(x)*dx,y-J(y)*dy,z-K(z)*dz);
    }
    
    double TricubicGrid::Dxyz(const double x, const double y, const double z) 
    {
        Tricubic* p = GetTricubic(x,y,z);
        return p->Dxyz(x-I(x)*dx,y-J(y)*dy,z-K(z)*dz);
    }
    
    double TricubicGrid::Dxzz(const double x, const double y, const double z) 
    {
        Tricubic* p = GetTricubic(x,y,z);
        return p->Dxzz(x-I(x)*dx,y-J(y)*dy,z-K(z)*dz);
    }
    
    double TricubicGrid::Dyyy(const double x, const double y, const double z) 
    {
        Tricubic* p = GetTricubic(x,y,z);
        return p->Dyyy(x-I(x)*dx,y-J(y)*dy,z-K(z)*dz);
    }

    double TricubicGrid::Dyyz(const double x, const double y, const double z) 
    {
        Tricubic* p = GetTricubic(x,y,z);
        return p->Dyyz(x-I(x)*dx,y-J(y)*dy,z-K(z)*dz);
    }

    double TricubicGrid::Dyzz(const double x, const double y, const double z) 
    {
        Tricubic* p = GetTricubic(x,y,z);
        return p->Dyzz(x-I(x)*dx,y-J(y)*dy,z-K(z)*dz);
    }

    double TricubicGrid::Dzzz(const double x, const double y, const double z) 
    {
        Tricubic* p = GetTricubic(x,y,z);
        return p->Dzzz(x-I(x)*dx,y-J(y)*dy,z-K(z)*dz);
    }

    void TricubicGrid::BestGuess(const double f000, const double f100, 
                                 const double f010, const double f110,
                                 const double f001, const double f101,
                                 const double f011, const double f111,
                                 const int boxscore, const int whichcorner, double* a) const
    {
        int score = boxscore;
        if (score < 0)
            score = (f000 >= 0. ? 1 : 0) + (f100 >= 0. ? 2 : 0) + (f010 >= 0. ? 4 : 0)
                + (f110 >= 0. ? 8 : 0) + (f001 >= 0. ? 16 : 0) + (f101 >= 0. ? 32 : 0)
                + (f011 >= 0. ? 64 : 0) + (f111 >= 0. ? 128 : 0);

        if (score > 0 && score < 255) {
            double alpha[4], beta[4], gamma[4];
            double del[3] = {dx, dy, dz};
            double f[2][2][2] = {f000, f001, f010, f011, f100, f101, f110, f111};
            int dir[3] = {1,1,1};
            int rot[3] = {0,0,0};
            int rotmap[3][8] = {{2, 3, 6, 7, 0, 1, 4, 5}, {4, 0, 6, 2, 5, 1, 7, 3}, 
                                {1, 3, 0, 2, 5, 7, 4, 6}};
            int corner = whichcorner;

            if (score > 127)
                score = 255-score;
            // reduce cases to core cases
            switch (score) {
                // case 1 no change
            case 2:
                score = 1; dir[0] = -1; break;
                // case 3 no change
            case 4:
                score = 1; dir[1] = -1; break;
            case 5:
                score = 3; rot[2] = 1; break;
                // case 6 no change
                // case 7 no change
            case 8:
                score = 1; dir[0] = -1; dir[1] = -1; break;
            case 9:
                score = 6; dir[1] = -1; break;
            case 10:
                score = 3; rot[2] = 1; dir[1] = -1; break;
            case 11:
                score = 31; dir[1] = -1; dir[2] = -1; break;
            case 12:
                score = 3; dir[1] = -1; break;
            case 13:
                score = 7; dir[1] = -1; break;
            case 14:
                score = 31; dir[2] = -1; break;
                // case 15 no change
            case 16:
                score = 1; dir[2] = -1; break;
            case 17:
                score = 3; rot[1] = 1; dir[2] = -1; break;
            case 18:
                score = 6; rot[0] = 1; dir[1] = -1; break;
            case 19:
                score = 7; rot[0] = 1; dir[1] = -1; break;
            case 20:
                score = 6; rot[1] = 1; dir[2] = -1; break;
            case 21:
                score = 7; rot[1] = 1; dir[2] = -1; break;
                // case 22 no change
                // case 23 no change
                // case 24 no change
                // case 25 no change
            case 26:
                score = 61; rot[2] = 1; dir[2] = -1; break;
                // case 27 no change
            case 28:
                score = 61; dir[0] = -1; dir[2] = -1; break;
            case 29:
                score = 27; rot[2] = 1; dir[0] = -1; break;
                // case 30 no change
                // case 31 no change
            case 32:
                score = 1; dir[0] = -1; dir[2] = -1; break;
            case 33:
                score = 6; rot[0] = 1; break;
            case 34:
                score = 3; rot[1] = 1; break;
            case 35:
                score = 31; rot[0] = 1; dir[2] = -1; break;
            case 36:
                score = 24; dir[1] = -1; dir[2] = -1; break;
            case 37:
                score = 61; rot[2] = 1; dir[1] = -1; dir[2] = -1; break;
            case 38:
                score = 25; dir[0] = -1; break;
            case 39:
                score = 27; dir[0] = -1; break;
            case 40:
                score = 6; rot[1] = 1; break;
            case 41:
                score = 22; dir[0] = -1; break;
            case 42:
                score = 7; rot[1] = 1; break;
            case 43:
                score = 23; dir[1] = -1; dir[2] = -1; break;
            case 44:
                score = 61; dir[2] = -1; break;
            case 45:
                score = 30; dir[0] = -1; break;
            case 46:
                score = 27; rot[2] = 1; dir[0] = -1; dir[2] = -1; break;
            case 47:
                score = 7; dir[1] = -1; dir[2] = -1; break;
            case 48:
                score = 3; dir[2] = -1; break;
            case 49:
                score = 7; rot[0] = 1; break;
            case 50:
                score = 31; rot[0] = 1; dir[1] = -1; dir[2] = -1; break;
            case 51:
                score = 15; rot[0] = 1; break;
            case 52:
                score = 61; dir[0] = -1; dir[1] = -1; break;
            case 53:
                score = 27; rot[1] = 1; dir[0] = -1; dir[2] = -1; break;
            case 54:
                score = 30; rot[0] = 1; dir[1] = -1; break;
            case 55:
                score = 31; rot[0] = 1; dir[1] = -1; break;
            case 56:
                score = 61; dir[1] = -1; break;
            case 57:
                score = 30; rot[0] = 1; dir[0] = -1; dir[1] = -1; break;
            case 58:
                score = 27; rot[1] = 1; dir[0] = -1; break;
            case 59:
                score = 7; rot[0] = 1; dir[2] = -1; break;
                // case 60 no change
                // case 61 no change
            case 62:
                score = 61; dir[0] = -1; break;
            case 63:
                score = 3; dir[1] = -1; dir[2] = -1; break;
            case 64:
                score = 1; dir[1] = -1; dir[2] = -1; break;
            case 65:
                score = 6; rot[1] = 1; dir[1] = -1; dir[2] = -1; break;
            case 66:
                score = 24; dir[1] = -1; break;
            case 67:
                score = 61; dir[0] = -1; dir[1] = -1; dir[2] = -1; break;
            case 68:
                score = 3; rot[1] = 1; dir[1] = -1; dir[2] = -1; break;
            case 69:
                score = 7; rot[1] = 1; dir[1] = -1; dir[2] = -1; break;
            case 70:
                score = 25; dir[1] = -1; break;
            case 71:
                score = 27; rot[2] = 1; break;
            case 72:
                score = 6; rot[0] = 1; dir[1] = -1; dir[2] = -1; break;
            case 73:
                score = 22; dir[1] = -1; break;
            case 74:
                score = 61; rot[2] = 1; dir[0] = -1; dir[2] = -1; break;
            case 75:
                score = 30; dir[1] = -1; break;
            case 76:
                score = 7; rot[0] = 1; dir[1] = -1; dir[2] = -1; break;
            case 77:
                score = 23; dir[1] = -1; break;
            case 78:
                score = 27; dir[2] = -1; break;
            case 79:
                score = 31; dir[1] = -1; break;
            case 80:
                score = 3; rot[2] = 1; dir[2] = -1; break;
            case 81:
                score = 31; rot[1] = 1; dir[1] = -1; break;
            case 82:
                score = 61; rot[2] = 1; dir[1] = -1; break;
            case 83:
                score = 27; rot[1] = 1; dir[2] = -1; break;
            case 84:
                score = 31; rot[1] = 1; break;
            case 85:
                score = 15; rot[1] = 1; break;
            case 86:
                score = 30; rot[1] = 1; break;
            case 87:
                score = 31; rot[1] = 1; dir[2] = -1; break;
            case 88:
                score = 61; rot[2] = 1; dir[0] = -1; dir[1] = -1; break;
            case 89:
                score = 30; rot[1] = 1; dir[1] = -1; break;
            case 90:
                score = 60; rot[2] = 1; break;
            case 91:
                score = 61; rot[2] = 1; dir[0] = -1; break;
            case 92:
                score = 27; rot[1] = 1; break;
            case 93:
                score = 31; rot[1] = 1; dir[1] = -1; dir[2] = -1; break;
            case 94:
                score = 61; rot[2] = 1; break;
            case 95:
                score = 3; rot[2] = 1; dir[1] = -1; dir[2] = -1; break;
            case 96:
                score = 6; dir[2] = -1; break;
            case 97:
                score = 22; dir[2] = -1; break;
            case 98:
                score = 25; dir[0] = -1; dir[2] = -1; break;
            case 99:
                score = 30; rot[0] = 1; break;
            case 100:
                score = 25; dir[1] = -1; dir[2] = -1; break;
            case 101:
                score = 30; rot[1] = 1; dir[0] = -1; break;
            case 102:
                score = 60; rot[1] = 1; break;
            case 103:
                score = 25; dir[0] = -1; dir[1] = -1; dir[2] = -1; break;
            case 104:
                score = 22; dir[0] = -1; dir[1] = -1; dir[2] = -1; break;
                // case 105 no change
            case 106:
                score = 30; rot[1] = 1; dir[0] = -1; dir[1] = -1; break;
            case 107:
                score = 22; dir[1] = -1; dir[2] = -1; break;
            case 108:
                score = 30; rot[0] = 1; dir[0] = -1; break;
            case 109:
                score = 22; dir[0] = -1; dir[2] = -1; break;
            case 110:
                score = 25; dir[2] = -1; break;
            case 111:
                score = 6; dir[1] = -1; dir[2] = -1; break;
            case 112:
                score = 7; dir[2] = -1; break;
            case 113:
                score = 23; dir[2] = -1; break;
            case 114:
                score = 27; dir[0] = -1; dir[2] = -1; break;
            case 115:
                score = 31; rot[0] = 1; break;
            case 116:
                score = 27; rot[2] = 1; dir[2] = -1; break;
            case 117:
                score = 7; rot[1] = 1; dir[1] = -1; break;
            case 118:
                score = 25; dir[0] = -1; dir[1] = -1; break;
            case 119:
                score = 3; rot[1] = 1; dir[1] = -1; break;
            case 120:
                score = 30; dir[0] = -1; dir[1] = -1; break;
            case 121:
                score = 22; dir[0] = -1; dir[1] = -1; break;
            case 122:
                score = 61; rot[2] = 1; dir[0] = -1; dir[1] = -1; dir[2] = -1; break;
            case 123:
                score = 6; rot[0] = 1; dir[2] = -1; break;
            case 124:
                score = 61; dir[1] = -1; dir[2] = -1; break;
            case 125:
                score = 6; rot[1] = 1; dir[1] = -1; break;
            case 126:
                score = 24; dir[2] = -1; break;
            case 127:
                score = 1; dir[0] = -1; dir[1] = -1; dir[2] = -1; break;
            }
            for (int r=0; r<rot[0]; ++r) {
                double temp = f[0][0][0];
                f[0][0][0] = f[0][0][1]; f[0][0][1] = f[0][1][1]; f[0][1][1] = f[0][1][0]; f[0][1][0] = temp;
                temp = f[1][0][0];
                f[1][0][0] = f[1][0][1]; f[1][0][1] = f[1][1][1]; f[1][1][1] = f[1][1][0]; f[1][1][0] = temp;
                swap(del[1], del[2]);
                corner = rotmap[0][corner];
            }
            for (int r=0; r<rot[1]; ++r) {
                double temp = f[0][0][0];
                f[0][0][0] = f[1][0][0]; f[1][0][0] = f[1][0][1]; f[1][0][1] = f[0][0][1]; f[0][0][1] = temp;
                temp = f[0][1][0];
                f[0][1][0] = f[1][1][0]; f[1][1][0] = f[1][1][1]; f[1][1][1] = f[0][1][1]; f[0][1][1] = temp;
                swap(del[0], del[2]);
                corner = rotmap[1][corner];
            }   
            for (int r=0; r<rot[2]; ++r) {
                double temp = f[0][0][0];
                f[0][0][0] = f[0][1][0]; f[0][1][0] = f[1][1][0]; f[1][1][0] = f[1][0][0]; f[1][0][0] = temp;
                temp = f[0][0][1];
                f[0][0][1] = f[0][1][1]; f[0][1][1] = f[1][1][1]; f[1][1][1] = f[1][0][1]; f[1][0][1] = temp;
                swap(del[0], del[1]);
                corner = rotmap[2][corner];
            }   
            if (dir[0] < 0) {
                for (int j=0; j<=1; ++j)
                    for (int k=0; k<=1; ++k)
                        swap(f[0][j][k], f[1][j][k]);
                corner = corner ^ 1;
            }
            if (dir[1] < 0) {
                for (int i=0; i<=1; ++i)
                    for (int k=0; k<=1; ++k)
                        swap(f[i][0][k], f[i][1][k]);
                corner = corner ^ 2;
            }
            if (dir[2] < 0) {
                for (int i=0; i<=1; ++i)
                    for (int j=0; j<=1; ++j)
                        swap(f[i][j][0], f[i][j][1]);
                corner = corner ^ 4;
            }
        
            switch (score) {
            case 1:     // single corner
                alpha[0] = -f[0][0][0]/(f[1][0][0]-f[0][0][0]);
                beta[0] = -f[0][0][0]/(f[0][1][0]-f[0][0][0]);
                gamma[0] = -f[0][0][0]/(f[0][0][1]-f[0][0][0]);
                switch(corner) {
                case 0:
                case 7:
                    a[0] = alpha[0]/3.; a[1] = beta[0]/3.; a[2] = gamma[0]/3.; break;
                case 1:
                    a[0] = alpha[0]; a[1] = 0.; a[2] = 0.; break;
                case 2:
                    a[0] = 0.; a[1] = beta[0]; a[2] = 0.; break;
                case 3:
                    a[0] = alpha[0]/2.; a[1] = beta[0]/2.; a[2] = 0.; break;
                case 4:
                    a[0] = 0.; a[1] = 0.; a[2] = gamma[0]; break;
                case 5:
                    a[0] = alpha[0]/2.; a[1] = 0.; a[2] = gamma[0]/2.; break;
                case 6:
                    a[0] = 0.; a[1] = beta[0]/2.; a[2] = gamma[0]/2.; break;
                } break;
            case 3:   // single edge
                beta[0] = -f[0][0][0]/(f[0][1][0]-f[0][0][0]);
                beta[1] = -f[1][0][0]/(f[1][1][0]-f[1][0][0]);
                gamma[0] = -f[0][0][0]/(f[0][0][1]-f[0][0][0]);
                gamma[1] = -f[1][0][0]/(f[1][0][1]-f[1][0][0]);
                switch(corner) {
                case 0:
                case 6:
                    a[0] = 0.; a[1] = beta[0]/2.; a[2] = gamma[0]/2.; break;
                case 1:
                case 7:
                    a[0] = 1.; a[1] = beta[1]/2.; a[2] = gamma[1]/2.; break;
                case 2:
                    a[0] = 0.; a[1] = beta[0]; a[2] = 0.; break;
                case 3:
                    a[0] = 1.; a[1] = beta[1]; a[2] = 0.; break;
                case 4:
                    a[0] = 0.; a[1] = 0.; a[2] = gamma[0]; break;
                case 5:
                    a[0] = 1.; a[1] = 0.; a[2] = gamma[1]; break;
                } break;          
            case 6:   // two corners, same facet
                alpha[0] = -f[0][0][0]/(f[1][0][0]-f[0][0][0]);
                alpha[1] = -f[0][1][0]/(f[1][1][0]-f[0][1][0]);
                beta[0] = -f[0][0][0]/(f[0][1][0]-f[0][0][0]);
                beta[1] = -f[1][0][0]/(f[1][1][0]-f[1][0][0]);
                gamma[0] = -f[1][0][0]/(f[1][0][1]-f[1][0][0]);
                gamma[1] = -f[0][1][0]/(f[0][1][1]-f[0][1][0]);
                switch(corner) {
                case 0:
                    if (alpha[0] < beta[0]) {
                        a[0] = alpha[0]; a[1] = 0.; a[2] = 0.;
                    } else {
                        a[0] = 0.; a[1] = beta[0]; a[2] = 0.;
                    } break;
                case 1:
                    a[0] = (2+alpha[0])/3.; a[1] = beta[1]/3.; a[2] = gamma[0]/3.; break;
                case 2:
                    a[0] = alpha[1]/3.; a[1] = (2+beta[0])/3.; a[2] = gamma[1]/3.; break;
                case 3:
                    if (alpha[1] < beta[1]) {
                        a[0] = alpha[1]; a[1] = 1.; a[2] = 0.;
                    } else {
                        a[0] = 1.; a[1] = beta[1]; a[2] = 0.;
                    } break;
                case 4:
                    if (alpha[0]-gamma[0] < beta[0]-gamma[1]) {
                        a[0] = (1+alpha[0])/2.; a[1] = 0.; a[2] = gamma[0]/2.;
                    } else {
                        a[0] = 0.; a[1] = (1+beta[0])/2.; a[2] = gamma[1]/2.;
                    } break;
                case 5:
                    a[0] = 1.; a[1] = 0.; a[2] = gamma[0]; break;
                case 6:
                    a[0] = 0.; a[1] = 1.; a[2] = gamma[1]; break;
                case 7:
                    if (alpha[1] + gamma[1] < beta[1] + gamma[0]) {
                        a[0] = 1.; a[1] = beta[1]/2.; a[2] = gamma[0]/2.;
                    } else {
                        a[0] = alpha[1]/2.; a[1] = 1.; a[2] = gamma[1]/2.;
                    } break;
                }
            case 7:
                alpha[0] = -f[0][1][0]/(f[1][1][0]-f[0][1][0]);
                beta[0] = -f[1][0][0]/(f[1][1][0]-f[1][0][0]);
                gamma[0] = -f[0][0][0]/(f[0][0][1]-f[0][0][0]);
                gamma[1] = -f[1][0][0]/(f[1][0][1]-f[1][0][0]);
                gamma[2] = -f[0][1][0]/(f[0][1][1]-f[0][1][0]);
                switch(corner) {
                case 0:
                case 4:
                    a[0] = 0.; a[1] = 0.; a[2] = gamma[0]; break;
                case 1:
                    a[0] = 1.; a[1] = beta[0]/2.; a[2] = gamma[1]/2.; break;
                case 2:
                    a[0] = alpha[0]/2.; a[1] = 1.; a[2] = gamma[2]/2.; break;
                case 3:
                    a[0] = (1+alpha[0])/2.; a[1] = (1+beta[0])/2.; a[2] = 0.; break;
                case 5:
                    a[0] = 1.; a[1] = 0.; a[2] = gamma[1]; break;
                case 6:
                    a[0] = 0.; a[1] = 1.; a[3] = gamma[2]; break;
                case 7:
                    a[0] = (2+alpha[0])/5.; a[1] = (2+beta[0])/5.; a[2] = (gamma[0]+gamma[1]+gamma[2])/5.; break;
                } break;
            case 15:
                gamma[0] = -f[0][0][0]/(f[0][0][1]-f[0][0][0]);
                gamma[1] = -f[1][0][0]/(f[1][0][1]-f[1][0][0]);
                gamma[2] = -f[0][1][0]/(f[0][1][1]-f[0][1][0]);
                gamma[3] = -f[1][1][0]/(f[1][1][1]-f[1][1][0]);
                switch(corner) {
                case 0:
                case 4:
                    a[0] = 0.; a[1] = 0.; a[2] = gamma[0]; break;
                case 1:
                case 5:
                    a[0] = 1.; a[1] = 0.; a[2] = gamma[1]; break;
                case 2:
                case 6:
                    a[0] = 0.; a[1] = 1.; a[2] = gamma[2]; break;
                case 3:
                case 7:
                    a[0] = 1.; a[1] = 1.; a[2] = gamma[3]; break;
                } break;
            case 22:
                alpha[0] = -f[0][0][0]/(f[1][0][0]-f[0][0][0]);
                alpha[1] = -f[0][1][0]/(f[1][1][0]-f[0][1][0]);
                alpha[2] = -f[0][0][1]/(f[1][0][1]-f[0][0][1]);
                beta[0] = -f[0][0][0]/(f[0][1][0]-f[0][0][0]);
                beta[1] = -f[1][0][0]/(f[1][1][0]-f[1][0][0]);
                beta[2] = -f[0][0][1]/(f[0][1][1]-f[0][0][1]);
                gamma[0] = -f[0][0][0]/(f[0][0][1]-f[0][0][0]);
                gamma[1] = -f[1][0][0]/(f[1][0][1]-f[1][0][0]);
                gamma[2] = -f[0][1][0]/(f[0][1][1]-f[0][1][0]);
                switch(corner) {
                case 0:
                    if (alpha[0] < beta[0]) {
                        if (alpha[0] < gamma[0]) {
                            a[0] = alpha[0]; a[1] = 0.; a[2] = 0.;
                        } else {
                            a[0] = 0.; a[1] = 0.; a[2] = gamma[0];
                        }
                    } else {
                        if (beta[0] < gamma[0]) {
                            a[0] = 0; a[1] = beta[0]; a[2] = 0.;
                        } else {
                            a[0] = 0.; a[1] = 0.; a[2] = gamma[0];
                        }
                    } break;
                case 1:
                    a[0] = (2+alpha[0])/3.; a[1] = beta[1]/3.; a[2] = gamma[1]/3.; break;
                case 2:
                    a[0] = alpha[1]/3.; a[1] = (2+beta[0])/3.; a[2] = gamma[2]/3.; break;
                case 3:
                    if (beta[1] < alpha[1]) {
                        a[0] = 1.; a[1] = beta[1]; a[2] = 0;
                    } else {
                        a[0] = alpha[1]; a[1] = 1.; a[2] = 0.;
                    } break;
                case 4:
                    a[0] = alpha[2]/3.; a[1] = beta[2]/3.; a[2] = (2+gamma[0])/3.; break;
                case 5:
                    if (gamma[1] < alpha[2]) {
                        a[0] = alpha[2]; a[1] = 0.; a[2] = 1.;
                    } else {
                        a[0] = 1.; a[1] = 0.; a[2] = gamma[1];
                    } break;
                case 6:
                    if (gamma[2] < beta[2]) {
                        a[0] = 0.; a[1] = beta[2]; a[2] = 1.; 
                    } else {
                        a[0] = 0.; a[1] = 1.; a[2] = gamma[2];
                    } break;
                case 7:
                    if (beta[1] < alpha[1]) {
                        a[0] = alpha[1]; a[1] = 1.; a[2] = 0.;
                    } else {
                        a[0] = 1.; a[1] = beta[1]; a[2] = 0.;
                    } break;
                } break;
            case 23:
                alpha[1] = -f[0][1][0]/(f[1][1][0]-f[0][1][0]);
                alpha[2] = -f[0][0][1]/(f[1][0][1]-f[0][0][1]);
                beta[1] = -f[1][0][0]/(f[1][1][0]-f[1][0][0]);
                beta[2] = -f[0][0][1]/(f[0][1][1]-f[0][0][1]);
                gamma[1] = -f[1][0][0]/(f[1][0][1]-f[1][0][0]);
                gamma[2] = -f[0][1][0]/(f[0][1][1]-f[0][1][0]);
                switch(corner) {
                case 0:
                case 7:
                    a[0] = (2+alpha[1]+alpha[2])/6.; a[1] = (2+beta[1]+beta[2])/6.; a[2] = (2+gamma[1]+gamma[2])/6.; break;
                case 1:
                    a[0] = 1.; a[1] = beta[1]/2.; a[2] = gamma[1]/2.; break;
                case 2:
                    a[0] = alpha[1]/2.; a[1] = 1.; a[2] = gamma[2]/2.; break;
                case 3:
                    a[0] = (1+alpha[1])/2.; a[1] = (1+beta[1])/2.; a[2] = 0.; break;
                case 4:
                    a[0] = alpha[2]/2.; a[1] = beta[2]/2.; a[2] = 1.; break;
                case 5:
                    a[0] = (1+alpha[2])/2.; a[1] = 0.; a[2] = (1+gamma[1])/2.; break;
                case 6:
                    a[0] = 0.; a[1] = (1+beta[2])/2.; a[2] = (1+gamma[2])/2.; break;
                } break;
            case 24:
                alpha[1] = -f[0][1][0]/(f[1][1][0]-f[0][1][0]);
                alpha[2] = -f[0][0][1]/(f[1][0][1]-f[0][0][1]);
                beta[1] = -f[1][0][0]/(f[1][1][0]-f[1][0][0]);
                beta[2] = -f[0][0][1]/(f[0][1][1]-f[0][0][1]);
                gamma[0] = -f[0][0][0]/(f[0][0][1]-f[0][0][0]);
                gamma[3] = -f[1][1][0]/(f[1][1][1]-f[1][1][0]);
                switch(corner) {
                case 0:
                    a[0] = 0.; a[1] = 0.; a[2] = gamma[0]; break;
                case 1:
                    a[0] = 1.; a[1] = beta[1]; a[2] = 0.; break;
                case 2:
                    a[0] = alpha[1]; a[1] = 1.; a[2] = 0.; break;
                case 3:
                    a[0] = (2+alpha[1])/3.; a[1] = (2+beta[1])/3.; a[2] = gamma[3]/3.; break;
                case 4:
                    a[0] = alpha[2]/3.; a[1] = beta[2]/3.; a[2] = (2+gamma[0])/3.; break;
                case 5:
                    a[0] = alpha[2]; a[1] = 0.; a[2] = 1.; break;
                case 6:
                    a[0] = 0.; a[1] = beta[2]; a[2] = 1.; break;
                case 7:
                    a[0] = 1.; a[1] = 1.; a[2] = gamma[3]; break;
                } break;
            case 25:
                alpha[0] = -f[0][0][0]/(f[1][0][0]-f[0][0][0]);
                alpha[1] = -f[0][1][0]/(f[1][1][0]-f[0][1][0]);
                alpha[2] = -f[0][0][1]/(f[1][0][1]-f[0][0][1]);
                beta[0] = -f[0][0][0]/(f[0][1][0]-f[0][0][0]);
                beta[1] = -f[1][0][0]/(f[1][1][0]-f[1][0][0]);
                beta[2] = -f[0][0][1]/(f[0][1][1]-f[0][0][1]);
                gamma[3] = -f[1][1][0]/(f[1][1][1]-f[1][1][0]);
                switch(corner) {
                case 0:
                    a[0] = alpha[0]/2.; a[1] = beta[0]/2.; a[2] = 0.; break;
                case 1:
                    if (1-alpha[0] < beta[1]) {
                        a[0] = alpha[0]; a[1] = 0.; a[2] = 0.; 
                    } else {
                        a[0] = 1.; a[1] = beta[1]; a[2] = 0.;
                    } break;
                case 2:
                    if (alpha[1] < 1-beta[0]) {
                        a[0] = alpha[1]; a[1] = 1.; a[2] = 0.; 
                    } else {
                        a[0] = 0.; a[1] = beta[0]; a[2] = 0.;
                    } break;
                case 3:
                    a[0] = (2+alpha[1])/3.; a[1] = (3+beta[1])/3.; a[2] = gamma[3]/3.; break;
                case 4:
                    a[0] = alpha[2]/2.; a[1] = beta[2]/2.; a[2] = 1.; break;
                case 5:
                    a[0] = alpha[2]; a[1] = 0.; a[2] = 1.; break;
                case 6:
                    a[0] = 0.; a[1] = beta[2]; a[2] = 1.; break;
                case 7:
                    a[0] = 1.; a[1] = 1.; a[2] = gamma[3]; break;
                } break;
            case 27:
                alpha[1] = -f[0][1][0]/(f[1][1][0]-f[0][1][0]);
                alpha[2] = -f[0][0][1]/(f[1][0][1]-f[0][0][1]);
                beta[0] = -f[0][0][0]/(f[0][1][0]-f[0][0][0]);
                beta[2] = -f[0][0][1]/(f[0][1][1]-f[0][0][1]);
                gamma[1] = -f[1][0][0]/(f[1][0][1]-f[1][0][0]);
                gamma[3] = -f[1][1][0]/(f[1][1][1]-f[1][1][0]);
                switch(corner) {
                case 0:
                    a[0] = 0.; a[1] = beta[0]; a[2] = 0.; break;
                case 1:
                    a[0] = 1.; a[1] = 0.; a[2] = gamma[1]; break;
                case 2:
                    a[0] = alpha[1]/2.; a[1] = (1+beta[0])/2.; a[2] = 0.; break;
                case 3:
                    a[0] = (1+alpha[1])/2.; a[1] = 1.; a[2] = gamma[3]/2.; break;
                case 4:
                    a[0] = alpha[2]/2.; a[1] = beta[2]/2.; gamma[3] = 1.; break;
                case 5:
                    a[0] = (1+alpha[2])/2.; a[1] = 0.; a[2] = (1.+gamma[1])/2.; break;
                case 6:
                    a[0] = 0.; a[1] = beta[2]; a[2] = 1.; break;
                case 7:
                    a[0] = 1.; a[1] = 1.; a[2] = gamma[3]; break;
                } break;
            case 30:
                alpha[0] = -f[0][0][0]/(f[1][0][0]-f[0][0][0]);
                alpha[2] = -f[0][0][1]/(f[1][0][1]-f[0][0][1]);
                beta[0] = -f[0][0][0]/(f[0][1][0]-f[0][0][0]);
                beta[2] = -f[0][0][1]/(f[0][1][1]-f[0][0][1]);
                gamma[0] = -f[0][0][0]/(f[0][0][1]-f[0][0][0]);
                gamma[1] = -f[1][0][0]/(f[1][0][1]-f[1][0][0]);
                gamma[2] = -f[0][1][0]/(f[0][1][1]-f[0][1][0]);
                gamma[3] = -f[1][1][0]/(f[1][1][1]-f[1][1][0]);
                switch(corner) {
                case 0:
                    a[0] = alpha[0]/2.; a[1] = beta[0]/2.; a[2] = 0.; break;
                case 1:
                    a[0] = (1+alpha[0])/2.; a[1] = 0.; a[2] = gamma[1]/2.; break;
                case 2:
                    a[0] = 0.; a[1] = (1+beta[0])/2.; a[2] = gamma[2]/2.; break;
                case 3:
                    a[0] = 1.; a[1] = 1.; a[2] = gamma[3]; break;
                case 4:
                    a[0] = alpha[2]/3.; a[1] = beta[2]/3.; a[2] = (2+gamma[0])/3.; break;
                case 5:
                    if (alpha[0] > gamma[1]) {
                        a[0] = alpha[0]; a[1] = 0.; a[2] = 1.; 
                    } else {
                        a[0] = 1.; a[1] = 0.; a[2] = gamma[1];
                    } break;
                case 6:
                    a[0] = 0.; a[1] = 1.; a[2] = gamma[2]; break;
                case 7:
                    a[0] = 1.; a[1] = 1.; a[2] = gamma[3]; break;
                } break;
            case 31:
                alpha[2] = -f[0][0][1]/(f[1][0][1]-f[0][0][1]);
                beta[2] = -f[0][0][1]/(f[0][1][1]-f[0][0][1]);
                gamma[1] = -f[1][0][0]/(f[1][0][1]-f[1][0][0]);
                gamma[2] = -f[0][1][0]/(f[0][1][1]-f[0][1][0]);
                gamma[3] = -f[1][1][0]/(f[1][1][1]-f[1][1][0]);
                switch(corner) {
                case 0:
                    a[0] = (2+alpha[2])/5.; a[1] = (2+beta[2])/5.; a[2] = (2+gamma[1]+gamma[2]+gamma[3])/5.; break;
                case 1:
                    a[0] = 1.; a[1] = 0.; a[2] = gamma[1]; break;
                case 2:
                    a[0] = 0.; a[1] = 1.; a[2] = gamma[2]; break;
                case 3:
                    a[0] = 1.; a[1] = 1.; a[2] = gamma[3]; break;
                case 4:
                    a[0] = alpha[2]/2.; a[1] = beta[2]/2.; a[2] = 1.; break;
                case 5:
                    a[0] = (1+alpha[2])/2.; a[1] = 0.; a[2] = (1+gamma[1])/2.; break;
                case 6:
                    a[0] = 0.; a[1] = (1+beta[2])/2.; a[2] = (1+gamma[2])/2.; break;
                case 7:
                    a[0] = 1.; a[1] = 1; a[2] = gamma[3]; break;
                } break;
            case 60:
                beta[0] = -f[0][0][0]/(f[0][1][0]-f[0][0][0]);
                beta[1] = -f[1][0][0]/(f[1][1][0]-f[1][0][0]);
                beta[2] = -f[0][0][1]/(f[0][1][1]-f[0][0][1]);
                beta[3] = -f[1][0][1]/(f[1][1][1]-f[1][0][1]);
                gamma[0] = -f[0][0][0]/(f[0][0][1]-f[0][0][0]);
                gamma[1] = -f[1][0][0]/(f[1][0][1]-f[1][0][0]);
                gamma[2] = -f[0][1][0]/(f[0][1][1]-f[0][1][0]);
                gamma[3] = -f[1][1][0]/(f[1][1][1]-f[1][1][0]);
                switch(corner) {
                case 0:
                    if (beta[0] < gamma[0]) {
                        a[0] = 0.; a[1] = beta[0]; a[2] = 0.;
                    } else {
                        a[0] = 0.; a[1] = 0.; a[2] = gamma[0]; 
                    } break;
                case 1:
                    if (beta[1] < gamma[1]) {
                        a[0] = 1.; a[1] = beta[1]; a[2] = 0.;
                    } else {
                        a[0] = 1.; a[1] = 0.; a[2] = gamma[1]; 
                    } break;
                case 2:
                    a[0] = 0.; a[1] = (1+beta[0])/2.; a[2] = gamma[2]/2.; break;
                case 3:
                    a[0] = 1.; a[1] = (1+beta[1])/2.; a[2] = gamma[3]/2.; break;
                case 4:
                    a[0] = 0.; a[1] = beta[2]/2.; a[2] = (1+gamma[0])/2.; break;
                case 5:
                    a[0] = 1.; a[1] = beta[3]/2.; a[2] = (1+gamma[1])/2.; break;
                case 6:
                    if (gamma[2] < beta[2]) {
                        a[0] = 0.; a[1] = beta[2]; a[2] = 1.; 
                    } else {
                        a[0] = 0.; a[1] = 1.; a[2] = gamma[2]; 
                    } break;
                case 7:
                    if (gamma[3] < beta[3]) {
                        a[0] = 1.; a[1] = beta[3]; a[2] = 1.;
                    } else {
                        a[0] = 1.; a[1] = 1.; a[2] = gamma[3];
                    } break;
                } break;
            case 61:
                alpha[0] = -f[0][0][0]/(f[1][0][0]-f[0][0][0]);
                beta[1] = -f[1][0][0]/(f[1][1][0]-f[1][0][0]);
                beta[2] = -f[0][0][1]/(f[0][1][1]-f[0][0][1]);
                beta[3] = -f[1][0][1]/(f[1][1][1]-f[1][0][1]);
                gamma[1] = -f[1][0][0]/(f[1][0][1]-f[1][0][0]);
                gamma[2] = -f[0][1][0]/(f[0][1][1]-f[0][1][0]);
                gamma[3] = -f[1][1][0]/(f[1][1][1]-f[1][1][0]);
                switch(corner) {
                case 0:
                    a[0] = alpha[0]; a[1] = 0.; a[2] = 0.; break;
                case 1:
                    a[0] = (2+alpha[0])/3.; a[1] = beta[1]/3.; a[2] = gamma[1]/3.; break;
                case 2:
                    a[0] = 0.; a[1] = 1.; a[2] = gamma[2]; break;
                case 3:
                    a[0] = 1.; a[1] = 1.; a[2] = gamma[3]; break;
                case 4:
                    a[0] = 0.; a[1] = beta[2]; a[2] = 1.; break;
                case 5:
                    if (1-gamma[1] < beta[3]) {
                        a[0] = 1.; a[1] = 0.; a[2] = gamma[1]; 
                    } else {
                        a[0] = 1.; a[1] = beta[3]; a[2] = 1.;
                    } break;
                case 6:
                    a[0] = 0.; a[1] = (1+beta[2])/2.; a[2] = (1+gamma[2])/2.; break;
                case 7:
                    a[0] = 1.; a[1] = (1+beta[3])/2.; a[2] = (1+gamma[3])/2.; break;
                } break;
            case 105:
                alpha[0] = -f[0][0][0]/(f[1][0][0]-f[0][0][0]);
                alpha[1] = -f[0][1][0]/(f[1][1][0]-f[0][1][0]);
                alpha[2] = -f[0][0][1]/(f[1][0][1]-f[0][0][1]);
                alpha[3] = -f[0][1][1]/(f[1][1][1]-f[0][1][1]);
                beta[0] = -f[0][0][0]/(f[0][1][0]-f[0][0][0]);
                beta[1] = -f[1][0][0]/(f[1][1][0]-f[1][0][0]);
                beta[2] = -f[0][0][1]/(f[0][1][1]-f[0][0][1]);
                beta[3] = -f[1][0][1]/(f[1][1][1]-f[1][0][1]);
                gamma[0] = -f[0][0][0]/(f[0][0][1]-f[0][0][0]);
                gamma[1] = -f[1][0][0]/(f[1][0][1]-f[1][0][0]);
                gamma[2] = -f[0][1][0]/(f[0][1][1]-f[0][1][0]);
                gamma[3] = -f[1][1][0]/(f[1][1][1]-f[1][1][0]);
                switch(corner) {
                case 0:
                    a[0] = alpha[0]/3.; a[1] = beta[0]/3.; a[2] = gamma[0]/3.; break;
                case 1:
                    a[0] = (2+alpha[0])/3.; a[1] = beta[1]/3.; a[2] = gamma[1]/3.; break;
                case 2:
                    a[0] = alpha[1]/3.; a[1] = (1+beta[0])/3.; a[2] = gamma[2]/3.; break;
                case 3:
                    a[0] = (2+alpha[1])/3.; a[1] = (1+beta[1])/3.; a[2] = gamma[3]/3.; break;
                case 4:
                    a[0] = alpha[2]/3.; a[1] = beta[2]/3.; a[2] = (2+gamma[0])/3.; break;
                case 5:
                    a[0] = (2+alpha[2])/3.; a[1] = beta[3]/3.; a[2] = (2+gamma[1])/3.; break;
                case 6:
                    a[0] = alpha[3]/3.; a[1] = (2+beta[2])/3.; a[2] = (2+gamma[2])/3.; break;
                case 7:
                    a[0] = (2+alpha[3])/3.; a[1] = (2+beta[3])/3.; a[2] = (2+gamma[3])/3.; break;
                } break;
            }
            if (dir[2] < 0) {
                a[2] = 1.-a[2];
            }
            if (dir[1] < 0) {
                a[1] = 1.-a[1];
            }
            if (dir[0] < 0) {
                a[0] = 1.-a[0];
            }
            for (int r=0; r<rot[2]; ++r) {
                double temp = a[0];
                a[0] = a[1];
                a[1] = 1-temp;
            }   
            for (int r=0; r<rot[1]; ++r) {
                double temp = a[0];
                a[0] = 1-a[2];
                a[2] = a[0];
            }   
            for (int r=0; r<rot[0]; ++r) {
                double temp = a[1];
                a[1] = a[2];
                a[2] = 1-a[1];
            }

            for (int d=0; d<3; ++d) 
                a[d] *= del[d];
        }
    }

    
    double TricubicGrid::LocalDist(const double x, const double y, const double z,
                                   double& ax, double& ay, double& az, char& clean)
    {
        ax = x;
        ay = y;
        double tol = 1.e-3;
        double p, px, py, pz, tax, tay, taz;
        double delta[6] = {1000.,1000.,1000.,1000.,1000.,1000.};
        int itcount = 0;
        double resid = 1000.;
        double oldresid = 2000.;
                
        while (resid > tol*dx*dy && (itcount < 20 || resid < oldresid)) {
                        
            // Find zero
                        
            p = F(ax,ay,az);
            px = Dx(ax,ay,az);
            py = Dy(ax,ay,az);
            pz = Dz(ax,ay,az);
            delta[0] = -p*px/(px*px+py*py+pz*pz);
            delta[1] = -p*py/(px*px+py*py+pz*pz);
            delta[2] = -p*pz/(px*px+py*py+pz*pz);
            tax = ax+delta[0];
            tay = ay+delta[1];
            taz = az+delta[2];
                        
            // Find min distance
                        
            p = F(tax,tay,taz);
            px = Dx(tax,tay,taz);
            py = Dy(tax,tay,taz);
            double denom = px*px+py*py+pz*pz;
            delta[3] = ((py*py+pz*pz)*(x-tax)-px*(py*(y-tay)+pz*(z-taz)))/denom;
            delta[4] = ((px*px+pz*pz)*(y-tay)-py*(px*(x-tax)+pz*(z-taz)))/denom;
            delta[5] = ((px*px+py*py)*(z-taz)-pz*(px*(x-tax)+py*(y-tay)))/denom;
                        
            ax += delta[0]+delta[2];
            ay += delta[1]+delta[3];
            az += delta[2]+delta[5];
            ++itcount;
            oldresid = resid;
            resid = sqrt(sqr(delta[0])+sqr(delta[1])+sqr(delta[2])+sqr(delta[3]));
        }
                
        //    clean = itcount < 20;
        clean = resid < tol*dx*dy && ax >= 0. && ax <= dx && ay >= 0. && ay <= dy
            && az >= 0. && az <= dz;
        
        return sqrt(sqr(x-ax)+sqr(y-ay)+sqr(z-az));
    }

    double TricubicGrid::LocalDistNewton(const double x, const double y, const double z, 
                                         double& ax, double& ay, double& az, char& clean)
    {
        double tol = 1.e-3;
        double p, px, py, pz, pxx, pxy, pxz, pyy, pyz, pzz;
        double Jinv[4][4];
        double det;
        double delta[4];
        int itcount = 0;
        double resid = 1000.;
        double oldresid = 2000.;
        double lambda = 0.;
        double Phi[4];

        bool doublecheck = false;
        while (resid > tol*dx*dy && (itcount < 20 || resid < oldresid)) {
                        
            // Find zero
                        
            p = F(ax,ay,az);
            px = Dx(ax,ay,az);
            py = Dy(ax,ay,az);
            pz = Dz(ax,ay,az);
            pxx = lambda*Dxx(ax,ay,az);
            pxy = lambda*Dxy(ax,ay,az);
            pxz = lambda*Dxz(ax,ay,az);
            pyy = lambda*Dyy(ax,ay,az);
            pyz = lambda*Dyz(ax,ay,az);
            pzz = lambda*Dzz(ax,ay,az);
            det = (-pyy*pzz + pyz*pyz - 2*pyy - 2*pzz - 4)*px*px 
                + (2*pxy*pzz - 2*pxz*pyz + 4*pxy)*px*py 
                + (-2*pxy*pyz + 2*pxz*pyy + 4*pxz)*px*pz 
                + (-pxx*pzz + pxz*pxz - 2*pxx - 2*pzz - 4)*py*py 
                + (2*pxx*pyz - 2*pxy*pxz + 4*pyz)*pz*py 
                + (-pxx*pyy + pxy*pxy - 2*pxx - 2*pyy - 4)*pz*pz;
            Jinv[0][0] = (- py*py*pzz + 2* py*pyz*pz - pyy*pz*pz - 2*py*py - 2*pz*pz)/det;
            Jinv[0][1] = ( px*py*pzz - px*pyz*pz + pxy*pz*pz - pxz*py*pz + 2*px*py)/det;
            Jinv[0][2] = (- px*py*pyz + px*pyy*pz - pxy*py*pz + pxz*py*py + 2*px*pz)/det;
            Jinv[0][3] = (- px*pyy*pzz + px*pyz*pyz + pxy*py*pzz - pxy*pyz*pz - pxz*py*pyz 
                          + pxz*pyy*pz - 2*px*pyy - 2*px*pzz + 2*pxy*py + 2*pxz*pz - 4*px)/det;
            Jinv[1][0] = Jinv[0][1];
            Jinv[1][1] = (- px*px*pzz + 2*px*pxz*pz - pxx*pz*pz - 2*px*px - 2*pz*pz)/det;
            Jinv[1][2] = ( px*px*pyz - px*pxy*pz - px*pxz*py + pxx*py*pz + 2*py*pz)/det;
            Jinv[1][3] = ( px*pxy*pzz - px*pxz*pyz - pxx*py*pzz + pxx*pyz*pz - pxy*pxz*pz 
                           + pxz*pxz*py + 2*px*pxy - 2*pxx*py - 2*py*pzz + 2*pyz*pz - 4*py)/det;
            Jinv[2][0] = Jinv[0][2];
            Jinv[2][1] = Jinv[1][2];
            Jinv[2][2] = (- px*px*pyy + 2*px*pxy*py - pxx*py*py - 2*px*px - 2*py*py)/det;
            Jinv[2][3] = (- px*pxy*pyz + px*pxz*pyy + pxx*py*pyz - pxx*pyy*pz + pxy*pxy*pz 
                          - pxy*pxz*py + 2*px*pxz - 2*pxx*pz + 2*py*pyz - 2*pyy*pz - 4*pz)/det;
            Jinv[3][0] = Jinv[0][3];
            Jinv[3][1] = Jinv[1][3];
            Jinv[3][2] = Jinv[2][3];
            Jinv[3][3] = (pxx*pyy*pzz - pxx*pyz*pyz - pxy*pxy*pzz + 2*pxy*pxz*pyz - pxz*pxz*pyy 
                          + 2*pxx*pyy + 2*pxx*pzz - 2*pxy*pxy - 2*pxz*pxz + 2*pyy*pzz - 2*pyz*pyz 
                          + 4*pxx + 4*pyy + 4*pzz + 8)/det;
      
            Phi[0] = 2*(ax-x)+lambda*px;
            Phi[1] = 2*(ay-y)+lambda*py;
            Phi[2] = 2*(az-z)+lambda*pz;
            Phi[3] = p;
            for (int i=0; i<4; ++i) {
                delta[i] = 0.;
                for (int j=0; j<4; ++j)
                    delta[i] += -Jinv[i][j]*Phi[j];
            }
            if (fabs(delta[0]) > dx/2.) {
                for (int i=0; i<4; ++i) 
                    delta[i] *= dx/2./fabs(delta[0]);
            }
            if (fabs(delta[1]) > dy/2.) {
                for (int i=0; i<4; ++i) 
                    delta[i] *= dy/2./fabs(delta[1]);
            }
            if (fabs(delta[2]) > dz/2.) {
                for (int i=0; i<4; ++i) 
                    delta[i] *= dz/2./fabs(delta[2]);
            }
                                                
            ax += delta[0];
            ay += delta[1];
            az += delta[2];
            lambda += delta[3];
            ++itcount;
            oldresid = resid;
            resid = sqrt(sqr(delta[0])+sqr(delta[1])+sqr(delta[2])+sqr(delta[3]));
            if (!doublecheck && itcount >= 20 && resid > tol*dx*dy && oldresid < resid) {
                doublecheck = true;
                ax = x;
                ay = y;
                az = z;
                itcount = 0;
            }
        }
                
        //    clean = itcount < 20;
        clean = resid < tol*dx*dy;
        
        return sqrt(sqr(x-ax)+sqr(y-ay)+sqr(z-az));
    }

#if 0
    void TricubicGrid::FollowToBdry(const int dir, const double x, const double y, double& ax, double &ay)
    {
        Tricubic* p = GetTricubic(x,y);
        int i = I(x);
        int j = J(y);
        double lx = x-X(i);
        double ly = y-Y(j);
        p->FollowToBdry(dir, lx, ly, ax, ay);
    }
#endif
}
  
  
