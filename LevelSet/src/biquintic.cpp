/*************************************************
    biquintic.cpp

    $Header: biquintic.cpp,v 1.1 2000/05/31 11:05:46 chopp Exp $

    $Log:       biquintic.cpp,v $
Revision 1.1  2000/05/31  11:05:46  11:05:46  chopp (David Chopp)
Initial revision

*************************************************/

#include "biquintic.h"
#include <math.h>
#include <iostream>
#include <fstream>

namespace levelset {

    Biquintic::Biquintic(const Biquintic& p)
    {
        for (int i=0; i<4; ++i)
            for (int j=0; j<4; ++j)
                a[i][j] = p.a[i][j];
        wx = p.wx;
        wy = p.wy;
    }


    void Biquintic::Build(const double* g, const double dx, const double dy)
    {    
        // g[0] = f[0][0], g[1] = f[1][0], g[2] = f[0][1], g[3] = f[1][1]
        // g[4] = fx[0][0], g[5] = fx[1][0], g[6] = fx[0][1], g[7] = fx[1][1]
        // g[8] = fy[0][0], g[9] = fy[1][0], g[10] = fy[0][1], g[11] = fy[1][1]
        // g[12] = fxy[0][0], g[13] = fxy[1][0], g[14] = fxy[0][1], g[15] = fxy[1][1]
        // g[16] = fxx[0][0], g[17] = fxx[1][0], g[18] = fxx[0][1], g[19] = fxx[1][1]
        // g[20] = fyy[0][0], g[21] = fyy[1][0], g[22] = fyy[0][1], g[23] = fyy[1][1]
        // g[24] = fxxy[0][0], g[25] = fxxy[1][0], g[26] = fxxy[0][1], g[27] = fxxy[1][1]
        // g[28] = fxyy[0][0], g[29] = fxyy[1][0], g[30] = fxyy[0][1], g[31] = fxyy[1][1]
        // g[32] = fxxyy[0][0], g[33] = fxxyy[1][0], g[34] = fxxyy[0][1], g[35] = fxxyy[1][1]
        
        a[0][0] = g[0];
        a[0][1] = g[8];
        a[0][2] = g[20]/2.;
        a[0][3] = -(20 * g[0] + 12 * g[8] * dy + 8 * dy * g[10] 
                    - 20 * g[2] + 3 * dy * dy * g[20] - dy * dy * g[22]) 
            /dy/dy/dy / 2.;
        a[0][4] = (16 * g[8] * dy - 2 * dy * dy * g[22] + 30 * g[0] 
                   + 14 * dy * g[10] - 30 * g[2] + 3 * dy * dy * g[20])
            /dy/dy/dy/dy / 2.;
        a[0][5] = - (12 * g[0] + 6 * g[8] * dy + 6 * dy * g[10] 
                     - 12 * g[2] + dy * dy * g[20] - dy * dy * g[22])
            /dy/dy/dy/dy/dy / 2.;
        
        a[1][0] = g[4];
        a[1][1] = g[12];
        a[1][2] = g[28] / 2.;
        a[1][3] = - (20 * g[4] - dy * dy * g[30] + 12 * g[12] * dy 
                     - 20 * g[6] + 3 * dy * dy * g[28] + 8 * g[14] * dy)
            /dy/dy/dy / 2.;
        a[1][4] =  (30 * g[4] - 2 * dy * dy * g[30] + 16 * g[12] * dy 
                    - 30 * g[6] + 3 * dy * dy * g[28] + 14 * g[14] * dy) 
            /dy/dy/dy/dy/ 2.;
        a[1][5] = - (12 * g[4] - dy * dy * g[30] + 6 * g[12] * dy 
                     - 12 * g[6] + dy * dy * g[28] + 6 * g[14] * dy) 
            /dy/dy/dy/dy/dy / 2.;
        
        a[2][0] = g[16] / 2.;
        a[2][1] = g[24] / 2.;
        a[2][2] = g[32] / 4.;
        a[2][3] = - (3 * dy * dy * g[32] - 20 * g[18] + 20 * g[16] 
                     + 12 * dy * g[24] - dy * dy * g[34] + 8 * dy * g[26]) 
            /dy/dy/dy / 4.;
        a[2][4] =  (3 * dy * dy * g[32] - 30 * g[18] + 30 * g[16] 
                    + 16 * dy * g[24] - 2 * dy * dy * g[34] + 14 * dy * g[26])
            /dy/dy/dy/dy/ 4.;
        a[2][5] = - (dy * dy * g[32] - 12 * g[18] + 12 * g[16] 
                     + 6 * dy * g[24] - dy * dy * g[34] + 6 * dy * g[26]) 
            /dy/dy/dy/dy/dy/ 4.;
        
        a[3][0] = - (3 * dx * dx * g[16] + 20 * g[0] - dx * dx * g[17] 
                     + 12 * g[4] * dx - 20 * g[1] + 8 * dx * g[5]) 
            /dx/dx/dx / 2.;
        a[3][1] = - (20 * g[8] - dx * dx * g[25] + 12 * g[12] * dx 
                     - 20 * g[9] + 3 * dx * dx * g[24] + 8 * dx * g[13])
            /dx/dx/dx / 2.;
        a[3][2] = - (3 * dx * dx * g[32] - 20 * g[21] + 20 * g[20] 
                     + 12 * dx * g[28] - dx * dx * g[33] + 8 * dx * g[29]) 
            /dx/dx/dx / 4.;
        a[3][3] =  (400 * g[0] - 400 * g[1] - 400 * g[2] + 400 * g[3] 
                    + 144 * g[12] * dx * dy + 60 * dx * dx * g[16] 
                    + 20 * dx * dx * g[19] - 160 * dy * g[11] 
                    + 9 * dx * dx * dy * dy * g[32] + 36 * dx * dy * dy * g[28] 
                    + 240 * g[4] * dx - 160 * dx * g[7] + dx * dx * dy * dy * g[35] 
                    + 36 * dx * dx * dy * g[24] - 3 * dx * dx * dy * dy * g[34] 
                    + 160 * dy * g[10] + 96 * dx * dy * g[13] + 60 * dy * dy * g[20] 
                    - 60 * dy * dy * g[21] - 12 * dx * dy * dy * g[30] 
                    + 24 * dx * dx * dy * g[26] - 8 * dx * dx * dy * g[27] 
                    + 240 * g[8] * dy - 20 * dy * dy * g[22] + 20 * dy * dy * g[23] 
                    + 160 * dx * g[5] - 20 * dx * dx * g[17] 
                    - 3 * dx * dx * dy * dy * g[33] - 60 * dx * dx * g[18]
                    - 240 * dy * g[9] - 12 * dx * dx * dy * g[25] 
                    + 96 * dx * dy * g[14] + 24 * dx * dy * dy * g[29] 
                    - 8 * dx * dy * dy * g[31] + 64 * dx * dy * g[15]
                    - 240 * dx * g[6]) /dx/dx/dx/dy/dy/dy / 4.;
        a[3][4] = - (600 * g[0] - 600 * g[1] - 600 * g[2] + 600 * g[3] 
                     + 192 * g[12] * dx * dy + 90 * dx * dx * g[16] 
                     + 30 * dx * dx * g[19] - 280 * dy * g[11] 
                     + 9 * dx * dx * dy * dy * g[32] + 36 * dx * dy * dy * g[28] 
                     + 360 * g[4] * dx - 240 * dx * g[7] 
                     + 2 * dx * dx * dy * dy * g[35] + 48 * dx * dx * dy * g[24] 
                     - 6 * dx * dx * dy * dy * g[34] + 280 * dy * g[10] 
                     + 128 * dx * dy * g[13] + 60 * dy * dy * g[20] 
                     - 60 * dy * dy * g[21] - 24 * dx * dy * dy * g[30] 
                     + 42 * dx * dx * dy * g[26] - 14 * dx * dx * dy * g[27] 
                     + 320 * g[8] * dy - 40 * dy * dy * g[22] + 40 * dy * dy * g[23] 
                     + 240 * dx * g[5] - 30 * dx * dx * g[17] 
                     - 3 * dx * dx * dy * dy * g[33] - 90 * dx * dx * g[18] 
                     - 320 * dy * g[9] - 16 * dx * dx * dy * g[25] 
                     + 168 * dx * dy * g[14] + 24 * dx * dy * dy * g[29] 
                     - 16 * dx * dy * dy * g[31] + 112 * dx * dy * g[15] 
                     - 360 * dx * g[6]) /dx/dx/dx/dy/dy/dy/dy / 4.;
        a[3][5] =  (240 * g[0] - 240 * g[1] - 240 * g[2] + 240 * g[3] 
                    + 72 * g[12] * dx * dy + 36 * dx * dx * g[16] 
                    + 12 * dx * dx * g[19] - 120 * dy * g[11] 
                    + 3 * dx * dx * dy * dy * g[32] + 12 * dx * dy * dy * g[28] 
                    + 144 * g[4] * dx - 96 * dx * g[7] + dx * dx * dy * dy * g[35] 
                    + 18 * dx * dx * dy * g[24] - 3 * dx * dx * dy * dy * g[34] 
                    + 120 * dy * g[10] + 48 * dx * dy * g[13] + 20 * dy * dy * g[20] 
                    - 20 * dy * dy * g[21] - 12 * dx * dy * dy * g[30]
                    + 18 * dx * dx * dy * g[26] - 6 * dx * dx * dy * g[27]
                    + 120 * g[8] * dy - 20 * dy * dy * g[22] + 20 * dy * dy * g[23] 
                    + 96 * dx * g[5] - 12 * dx * dx * g[17] 
                    - dx * dx * dy * dy * g[33] - 36 * dx * dx * g[18] 
                    - 120 * dy * g[9] - 6 * dx * dx * dy * g[25]
                    + 72 * dx * dy * g[14] + 8 * dx * dy * dy * g[29] 
                    - 8 * dx * dy * dy * g[31] + 48 * dx * dy * g[15]
                    - 144 * dx * g[6]) /dx/dx/dx/dy/dy/dy/dy/dy / 4.;
        
        a[4][0] =  (16 * g[4] * dx + 30 * g[0] - 2 * dx * dx * g[17] - 30 * g[1] 
                    + 3 * dx * dx * g[16] + 14 * dx * g[5]) /dx/dx/dx/dx / 2.;
        a[4][1] =  (30 * g[8] - 2 * dx * dx * g[25] + 16 * g[12] * dx 
                    - 30 * g[9] + 3 * dx * dx * g[24] + 14 * dx * g[13]) 
            /dx/dx/dx/dx / 2.;
        a[4][2] =  (3 * dx * dx * g[32] - 30 * g[21] + 30 * g[20] 
                    + 16 * dx * g[28] - 2 * dx * dx * g[33] + 14 * dx * g[29]) 
            /dx/dx/dx/dx / 4.;
        a[4][3] = - (600 * g[0] - 600 * g[1] - 600 * g[2] + 600 * g[3] 
                     + 192 * g[12] * dx * dy + 60 * dx * dx * g[16]
                     + 40 * dx * dx * g[19] - 240 * dy * g[11] 
                     + 9 * dx * dx * dy * dy * g[32] + 48 * dx * dy * dy * g[28] 
                     + 320 * g[4] * dx - 280 * dx * g[7] 
                     + 2 * dx * dx * dy * dy * g[35] + 36 * dx * dx * dy * g[24] 
                     - 3 * dx * dx * dy * dy * g[34] + 240 * dy * g[10] 
                     + 168 * dx * dy * g[13] + 90 * dy * dy * g[20] 
                     - 90 * dy * dy * g[21] - 16 * dx * dy * dy * g[30] 
                     + 24 * dx * dx * dy * g[26] - 16 * dx * dx * dy * g[27] 
                     + 360 * g[8] * dy - 30 * dy * dy * g[22] + 30 * dy * dy * g[23] 
                     + 280 * dx * g[5] - 40 * dx * dx * g[17] 
                     - 6 * dx * dx * dy * dy * g[33] - 60 * dx * dx * g[18] 
                     - 360 * dy * g[9] - 24 * dx * dx * dy * g[25] 
                     + 128 * dx * dy * g[14] + 42 * dx * dy * dy * g[29] 
                     - 14 * dx * dy * dy * g[31] + 112 * dx * dy * g[15] 
                     - 320 * dx * g[6]) /dx/dx/dx/dx/dy/dy/dy / 4.;
        a[4][4] =  (900 * g[0] - 900 * g[1] - 900 * g[2] + 900 * g[3] 
                    + 256 * g[12] * dx * dy + 90 * dx * dx * g[16] 
                    + 60 * dx * dx * g[19] - 420 * dy * g[11] 
                    + 9 * dx * dx * dy * dy * g[32] + 48 * dx * dy * dy * g[28] 
                    + 480 * g[4] * dx - 420 * dx * g[7] 
                    + 4 * dx * dx * dy * dy * g[35] + 48 * dx * dx * dy * g[24] 
                    - 6 * dx * dx * dy * dy * g[34] + 420 * dy * g[10] 
                    + 224 * dx * dy * g[13] + 90 * dy * dy * g[20] 
                    - 90 * dy * dy * g[21] - 32 * dx * dy * dy * g[30] 
                    + 42 * dx * dx * dy * g[26] - 28 * dx * dx * dy * g[27] 
                    + 480 * g[8] * dy - 60 * dy * dy * g[22] + 60 * dy * dy * g[23] 
                    + 420 * dx * g[5] - 60 * dx * dx * g[17] 
                    - 6 * dx * dx * dy * dy * g[33] - 90 * dx * dx * g[18]
                    - 480 * dy * g[9] - 32 * dx * dx * dy * g[25] 
                    + 224 * dx * dy * g[14] + 42 * dx * dy * dy * g[29] 
                    - 28 * dx * dy * dy * g[31] + 196 * dx * dy * g[15] 
                    - 480 * dx * g[6])/dx/dx/dx/dx/dy/dy/dy/dy / 4.;
        a[4][5] = - (360 * g[0] - 360 * g[1] - 360 * g[2] + 360 * g[3] 
                     + 96 * g[12] * dx * dy + 36 * dx * dx * g[16] 
                     + 24 * dx * dx * g[19] - 180 * dy * g[11] 
                     + 3 * dx * dx * dy * dy * g[32] + 16 * dx * dy * dy * g[28] 
                     + 192 * g[4] * dx - 168 * dx * g[7] 
                     + 2 * dx * dx * dy * dy * g[35] + 18 * dx * dx * dy * g[24]
                     - 3 * dx * dx * dy * dy * g[34] + 180 * dy * g[10]
                     + 84 * dx * dy * g[13] + 30 * dy * dy * g[20]
                     - 30 * dy * dy * g[21] - 16 * dx * dy * dy * g[30]
                     + 18 * dx * dx * dy * g[26] - 12 * dx * dx * dy * g[27] 
                     + 180 * g[8] * dy - 30 * dy * dy * g[22] + 30 * dy * dy * g[23]
                     + 168 * dx * g[5] - 24 * dx * dx * g[17] 
                     - 2 * dx * dx * dy * dy * g[33] - 36 * dx * dx * g[18] 
                     - 180 * dy * g[9] - 12 * dx * dx * dy * g[25] 
                     + 96 * dx * dy * g[14] + 14 * dx * dy * dy * g[29] 
                     - 14 * dx * dy * dy * g[31] + 84 * dx * dy * g[15]
                     - 192 * dx * g[6])/dx/dx/dx/dx/dy/dy/dy/dy/dy / 4.;
        
        a[5][0] = - (12 * g[0] - dx * dx * g[17] + 6 * g[4] * dx - 12 * g[1] 
                     + dx * dx * g[16] + 6 * dx * g[5])/dx/dx/dx/dx/dx / 2.;
        a[5][1] = - (12 * g[8] - dx * dx * g[25] + 6 * g[12] * dx - 12 * g[9] 
                     + dx * dx * g[24] + 6 * dx * g[13]) /dx/dx/dx/dx/dx / 2.;
        a[5][2] = - (dx * dx * g[32] - 12 * g[21] + 12 * g[20] + 6 * dx * g[28] 
                     - dx * dx * g[33] + 6 * dx * g[29])/dx/dx/dx/dx/dx / 4.;
        a[5][3] =  (240 * g[0] - 240 * g[1] - 240 * g[2] + 240 * g[3] 
                    + 72 * g[12] * dx * dy + 20 * dx * dx * g[16] 
                    + 20 * dx * dx * g[19] - 96 * dy * g[11] 
                    + 3 * dx * dx * dy * dy * g[32] + 18 * dx * dy * dy * g[28]
                    + 120 * g[4] * dx - 120 * dx * g[7] + dx * dx * dy * dy * g[35]
                    + 12 * dx * dx * dy * g[24] - dx * dx * dy * dy * g[34]
                    + 96 * dy * g[10] + 72 * dx * dy * g[13] + 36 * dy * dy * g[20] 
                    - 36 * dy * dy * g[21] - 6 * dx * dy * dy * g[30] 
                    + 8 * dx * dx * dy * g[26] - 8 * dx * dx * dy * g[27] 
                    + 144 * g[8] * dy - 12 * dy * dy * g[22] + 12 * dy * dy * g[23] 
                    + 120 * dx * g[5] - 20 * dx * dx * g[17] 
                    - 3 * dx * dx * dy * dy * g[33] - 20 * dx * dx * g[18] 
                    - 144 * dy * g[9] - 12 * dx * dx * dy * g[25] 
                    + 48 * dx * dy * g[14] + 18 * dx * dy * dy * g[29] 
                    - 6 * dx * dy * dy * g[31] + 48 * dx * dy * g[15] 
                    - 120 * dx * g[6])/dx/dx/dx/dx/dx/dy/dy/dy / 4.;
        a[5][4] = - (360 * g[0] - 360 * g[1] - 360 * g[2] + 360 * g[3] 
                     + 96 * g[12] * dx * dy + 30 * dx * dx * g[16] 
                     + 30 * dx * dx * g[19] - 168 * dy * g[11] 
                     + 3 * dx * dx * dy * dy * g[32] + 18 * dx * dy * dy * g[28] 
                     + 180 * g[4] * dx - 180 * dx * g[7] 
                     + 2 * dx * dx * dy * dy * g[35] + 16 * dx * dx * dy * g[24] 
                     - 2 * dx * dx * dy * dy * g[34] + 168 * dy * g[10] 
                     + 96 * dx * dy * g[13] + 36 * dy * dy * g[20] 
                     - 36 * dy * dy * g[21] - 12 * dx * dy * dy * g[30] 
                     + 14 * dx * dx * dy * g[26] - 14 * dx * dx * dy * g[27]
                     + 192 * g[8] * dy - 24 * dy * dy * g[22] + 24 * dy * dy * g[23]
                     + 180 * dx * g[5] - 30 * dx * dx * g[17] 
                     - 3 * dx * dx * dy * dy * g[33] - 30 * dx * dx * g[18] 
                     - 192 * dy * g[9] - 16 * dx * dx * dy * g[25]
                     + 84 * dx * dy * g[14] + 18 * dx * dy * dy * g[29] 
                     - 12 * dx * dy * dy * g[31] + 84 * dx * dy * g[15] 
                     - 180 * dx * g[6]) /dx/dx/dx/dx/dx/dy/dy/dy/dy / 4.;
        a[5][5] =  (144 * g[0] - 144 * g[1] - 144 * g[2] + 144 * g[3] 
                    + 36 * g[12] * dx * dy + 12 * dx * dx * g[16] 
                    + 12 * dx * dx * g[19] - 72 * dy * g[11]
                    + dx * dx * dy * dy * g[32] + 6 * dx * dy * dy * g[28] 
                    + 72 * g[4] * dx - 72 * dx * g[7] + dx * dx * dy * dy * g[35] 
                    + 6 * dx * dx * dy * g[24] - dx * dx * dy * dy * g[34] 
                    + 72 * dy * g[10] + 36 * dx * dy * g[13] + 12 * dy * dy * g[20] 
                    - 12 * dy * dy * g[21] - 6 * dx * dy * dy * g[30] 
                    + 6 * dx * dx * dy * g[26] - 6 * dx * dx * dy * g[27]
                    + 72 * g[8] * dy - 12 * dy * dy * g[22] + 12 * dy * dy * g[23] 
                    + 72 * dx * g[5] - 12 * dx * dx * g[17] 
                    - dx * dx * dy * dy * g[33] - 12 * dx * dx * g[18] 
                    - 72 * dy * g[9] - 6 * dx * dx * dy * g[25] 
                    + 36 * dx * dy * g[14] + 6 * dx * dy * dy * g[29] 
                    - 6 * dx * dy * dy * g[31] + 36 * dx * dy * g[15] 
                    - 72 * dx * g[6])/dx/dx/dx/dx/dx/dy/dy/dy/dy/dy / 4.;
        
        wx = dx;
        wy = dy;
    }

    void Biquintic::BuildwDeriv(const double f[6][6], const double dx, const double dy)
    {
        double g[36];
        g[0] = f[2][2];
        g[1] = f[3][2];
        g[2] = f[2][3];
        g[3] = f[3][3];
        
        g[4] = (f[0][2]-8*f[1][2]+8*f[3][2]-f[4][2])/12/dx;
        g[5] = (f[1][2]-8*f[2][2]+8*f[4][2]-f[5][2])/12/dx;
        g[6] = (f[0][3]-8*f[1][3]+8*f[3][3]-f[4][3])/12/dx;
        g[7] = (f[1][3]-8*f[2][3]+8*f[4][3]-f[5][3])/12/dx;
        
        g[8] = (f[2][0]-8*f[2][1]+8*f[2][3]-f[2][4])/12/dy;
        g[9] = (f[3][0]-8*f[3][1]+8*f[3][3]-f[3][4])/12/dy;
        g[10] = (f[2][1]-8*f[2][2]+8*f[2][4]-f[2][5])/12/dy;
        g[11] = (f[3][1]-8*f[3][2]+8*f[3][4]-f[3][5])/12/dy;
        
        g[12] = (f[0][0]-8*f[1][0]+8*f[3][0]-f[4][0]
                 -8*f[0][1]+64*f[1][1]-64*f[3][1]+8*f[4][1]
                 +8*f[0][3]-64*f[1][3]+64*f[3][3]-8*f[4][3]
                 -f[0][4]+8*f[1][4]-8*f[3][4]+f[4][4])/144./dx/dy;
        g[13] = (f[1][0]-8*f[2][0]+8*f[4][0]-f[5][0]
                 -8*f[1][1]+64*f[2][1]-64*f[4][1]+8*f[5][1]
                 +8*f[1][3]-64*f[2][3]+64*f[4][3]-8*f[5][3]
                 -f[1][4]+8*f[2][4]-8*f[4][4]+f[5][4])/144./dx/dy;
        g[14] = (f[0][1]-8*f[1][1]+8*f[3][1]-f[4][1]
                 -8*f[0][2]+64*f[1][2]-64*f[3][2]+8*f[4][2]
                 +8*f[0][4]-64*f[1][4]+64*f[3][4]-8*f[4][4]
                 -f[0][5]+8*f[1][5]-8*f[3][5]+f[4][5])/144./dx/dy;
        g[15] = (f[1][1]-8*f[2][1]+8*f[4][1]-f[5][1]
                 -8*f[1][2]+64*f[2][2]-64*f[4][2]+8*f[5][2]
                 +8*f[1][4]-64*f[2][4]+64*f[4][4]-8*f[5][4]
                 -f[1][5]+8*f[2][5]-8*f[4][5]+f[5][5])/144./dx/dy;
        
        g[16] = (-f[0][2]+16*f[1][2]-30*f[2][2]+16*f[3][2]-f[4][2])/12/dx/dx;
        g[17] = (-f[1][2]+16*f[2][2]-30*f[3][2]+16*f[4][2]-f[5][2])/12/dx/dx;
        g[18] = (-f[0][3]+16*f[1][3]-30*f[2][3]+16*f[3][3]-f[4][3])/12/dx/dx;
        g[19] = (-f[1][3]+16*f[2][3]-30*f[3][3]+16*f[4][3]-f[5][3])/12/dx/dx;
        
        g[20] = (-f[2][0]+16*f[2][1]-30*f[2][2]+16*f[2][3]-f[2][4])/12/dy/dy;
        g[21] = (-f[3][0]+16*f[3][1]-30*f[3][2]+16*f[3][3]-f[3][4])/12/dy/dy;
        g[22] = (-f[2][1]+16*f[2][2]-30*f[2][3]+16*f[2][4]-f[2][5])/12/dy/dy;
        g[23] = (-f[3][1]+16*f[3][2]-30*f[3][3]+16*f[3][4]-f[3][5])/12/dy/dy;
        
        g[24] = (-f[0][0]+16*f[1][0]-30*f[2][0]+16*f[3][0]-f[4][0]
                 +8*f[0][1]-128*f[1][1]+240*f[2][1]-128*f[3][1]+8*f[4][1]
                 -8*f[0][3]+128*f[1][3]-240*f[2][3]+128*f[3][3]-8*f[4][3]
                 +f[0][4]-16*f[1][4]+30*f[2][4]-16*f[3][4]+f[4][4])/144./dx/dx/dy;      
        g[25] = (-f[1][0]+16*f[2][0]-30*f[3][0]+16*f[4][0]-f[5][0]
                 +8*f[1][1]-128*f[2][1]+240*f[3][1]-128*f[4][1]+8*f[5][1]
                 -8*f[1][3]+128*f[2][3]-240*f[3][3]+128*f[4][3]-8*f[5][3]
                 +f[1][4]-16*f[2][4]+30*f[3][4]-16*f[4][4]+f[5][4])/144./dx/dx/dy;      
        g[26] = (-f[0][1]+16*f[1][1]-30*f[2][1]+16*f[3][1]-f[4][1]
                 +8*f[0][2]-128*f[1][2]+240*f[2][2]-128*f[3][2]+8*f[4][2]
                 -8*f[0][4]+128*f[1][4]-240*f[2][4]+128*f[3][4]-8*f[4][4]
                 +f[0][5]-16*f[1][5]+30*f[2][5]-16*f[3][5]+f[4][5])/144./dx/dx/dy;      
        g[27] = (-f[1][1]+16*f[2][1]-30*f[3][1]+16*f[4][1]-f[5][1]
                 +8*f[1][2]-128*f[2][2]+240*f[3][2]-128*f[4][2]+8*f[5][2]
                 -8*f[1][4]+128*f[2][4]-240*f[3][4]+128*f[4][4]-8*f[5][4]
                 +f[1][5]-16*f[2][5]+30*f[3][5]-16*f[4][5]+f[5][5])/144./dx/dx/dy;      
        
        g[28] = (-f[0][0]+16*f[0][1]-30*f[0][2]+16*f[0][3]-f[0][4]
                 +8*f[1][0]-128*f[1][1]+240*f[1][2]-128*f[1][3]+8*f[1][4]
                 -8*f[3][0]+128*f[3][1]-240*f[3][2]+128*f[3][3]-8*f[3][4]
                 +f[4][0]-16*f[4][1]+30*f[4][2]-16*f[4][3]+f[4][4])/144./dx/dy/dy;      
        g[29] = (-f[1][0]+16*f[1][1]-30*f[1][2]+16*f[1][3]-f[1][4]
                 +8*f[2][0]-128*f[2][1]+240*f[2][2]-128*f[2][3]+8*f[2][4]
                 -8*f[4][0]+128*f[4][1]-240*f[4][2]+128*f[4][3]-8*f[4][4]
                 +f[5][0]-16*f[5][1]+30*f[5][2]-16*f[5][3]+f[5][4])/144./dx/dy/dy;      
        g[30] = (-f[0][1]+16*f[0][2]-30*f[0][3]+16*f[0][4]-f[0][5]
                 +8*f[1][1]-128*f[1][2]+240*f[1][3]-128*f[1][4]+8*f[1][5]
                 -8*f[3][1]+128*f[3][2]-240*f[3][3]+128*f[3][4]-8*f[3][5]
                 +f[4][1]-16*f[4][2]+30*f[4][3]-16*f[4][4]+f[4][5])/144./dx/dy/dy;      
        g[31] = (-f[1][1]+16*f[1][2]-30*f[1][3]+16*f[1][4]-f[1][5]
                 +8*f[2][1]-128*f[2][2]+240*f[2][3]-128*f[2][4]+8*f[2][5]
                 -8*f[4][1]+128*f[4][2]-240*f[4][3]+128*f[4][4]-8*f[4][5]
                 +f[5][1]-16*f[5][2]+30*f[5][3]-16*f[5][4]+f[5][5])/144./dx/dy/dy;      
                                
        g[32] = (f[0][0]-16*f[1][0]+30*f[2][0]-16*f[3][0]+f[4][0]
                 -16*f[0][1]+256*f[1][1]-480*f[2][1]+256*f[3][1]-16*f[4][1]
                 +30*f[0][2]-480*f[1][2]+900*f[2][2]-480*f[3][2]+30*f[4][2]
                 -16*f[0][3]+256*f[1][3]-480*f[2][3]+256*f[3][3]-16*f[4][3]
                 +f[0][4]-16*f[1][4]+30*f[2][4]-16*f[3][4]+f[4][4])/144./dx/dx/dy/dy;
        g[33] = (f[1][0]-16*f[2][0]+30*f[3][0]-16*f[4][0]+f[5][0]
                 -16*f[1][1]+256*f[2][1]-480*f[3][1]+256*f[4][1]-16*f[5][1]
                 +30*f[1][2]-480*f[2][2]+900*f[3][2]-480*f[4][2]+30*f[5][2]
                 -16*f[1][3]+256*f[2][3]-480*f[3][3]+256*f[4][3]-16*f[5][3]
                 +f[1][4]-16*f[2][4]+30*f[3][4]-16*f[4][4]+f[5][4])/144./dx/dx/dy/dy;
        g[34] = (f[0][1]-16*f[1][1]+30*f[2][1]-16*f[3][1]+f[4][1]
                 -16*f[0][2]+256*f[1][2]-480*f[2][2]+256*f[3][2]-16*f[4][2]
                 +30*f[0][3]-480*f[1][3]+900*f[2][3]-480*f[3][3]+30*f[4][3]
                 -16*f[0][4]+256*f[1][4]-480*f[2][4]+256*f[3][4]-16*f[4][4]
                 +f[0][5]-16*f[1][5]+30*f[2][5]-16*f[3][5]+f[4][5])/144./dx/dx/dy/dy;
        g[35] = (f[1][1]-16*f[2][1]+30*f[3][1]-16*f[4][1]+f[5][1]
                 -16*f[1][2]+256*f[2][2]-480*f[3][2]+256*f[4][2]-16*f[5][2]
                 +30*f[1][3]-480*f[2][3]+900*f[3][3]-480*f[4][3]+30*f[5][3]
                 -16*f[1][4]+256*f[2][4]-480*f[3][4]+256*f[4][4]-16*f[5][4]
                 +f[1][5]-16*f[2][5]+30*f[3][5]-16*f[4][5]+f[5][5])/144./dx/dx/dy/dy;
        
        
        Build(g, dx, dy);
    }

    double Biquintic::operator()(const double x, const double y) const
    {
        double ans = 0;
    
        for (int i=0; i<6; ++i) 
            for (int j=0; j<6; ++j)
                ans += a[i][j]*pow(x,i)*pow(y,j);
    
        return ans;
    }

    double Biquintic::Dx(const double x, const double y) const
    {
        double ans = 0;
    
        for (int i=1; i<6; ++i) 
            for (int j=0; j<6; ++j)
                ans += i*a[i][j]*pow(x,i-1)*pow(y,j);
    
        return ans;
    }

    double Biquintic::Dy(const double x, const double y) const
    {
        double ans = 0;
    
        for (int i=0; i<6; ++i) 
            for (int j=1; j<6; ++j)
                ans += j*a[i][j]*pow(x,i)*pow(y,j-1);
    
        return ans;
    }

    double Biquintic::Dxx(const double x, const double y) const
    {
        double ans = 0;
    
        for (int i=2; i<6; ++i) 
            for (int j=0; j<6; ++j)
                ans += i*(i-1)*a[i][j]*pow(x,i-2)*pow(y,j);
    
        return ans;
    }

    double Biquintic::Dxy(const double x, const double y) const
    {
        double ans = 0;
    
        for (int i=1; i<6; ++i) 
            for (int j=1; j<6; ++j)
                ans += i*j*a[i][j]*pow(x,i-1)*pow(y,j-1);
    
        return ans;
    }

    double Biquintic::Dyy(const double x, const double y) const
    {
        double ans = 0;
    
        for (int i=0; i<6; ++i) 
            for (int j=2; j<6; ++j)
                ans += j*(j-1)*a[i][j]*pow(x,i)*pow(y,j-2);
    
        return ans;
    }

    double Biquintic::Dxxx(const double x, const double y) const
    {
        double ans = 0;
   
        for (int i=3; i<6; ++i)
            for (int j=0; j<6; ++j)
                ans += 6*a[i][j]*pow(x,i-3)*pow(y,j);
    
        return ans;
    }

    double Biquintic::Dxxy(const double x, const double y) const
    {
        double ans = 0;
    
        for (int i=2; i<6; ++i) 
            for (int j=1; j<6; ++j)
                ans += i*(i-1)*j*a[i][j]*pow(x,i-2)*pow(y,j-1);
    
        return ans;
    }

    double Biquintic::Dxyy(const double x, const double y) const
    {
        double ans = 0;
    
        for (int i=1; i<6; ++i) 
            for (int j=2; j<6; ++j)
                ans += i*j*(j-1)*a[i][j]*pow(x,i-1)*pow(y,j-2);
    
        return ans;
    }

    double Biquintic::Dyyy(const double x, const double y) const
    {
        double ans = 0;
    
        for (int i=0; i<6; ++i)
            for (int j=3; j<6; ++j) 
                ans += 6*a[i][3]*pow(x,i)*pow(y,j-3);
    
        return ans;
    }

    Biquintic Biquintic::Dx(void)
    {
        Biquintic d;
    
        for (int j=0; j<6; ++j) {
            d.a[5][j] = 0.;
            for (int i=1; i<6; ++i)
                d.a[i-1][j] = i*a[i][j];
        }
    
        return d;
    }

    Biquintic Biquintic::Dy(void)
    {
        Biquintic d;
    
        for (int i=0; i<6; ++i) {
            d.a[i][5] = 0.;
            for (int j=1; j<6; ++j)
                d.a[i][j-1] = j*a[i][j];
        }
    
        return d;
    }

    inline double sqr(const double x) {return x*x;}

    double Biquintic::LocalDist(const double x, const double y, 
                                double& ax, double& ay, char& clean) const
    {
        ax = x;
        ay = y;
        double tol = 1.e-3;
        double p, px, py, tax, tay;
        double delta[4] = {1000.,1000.,1000.,1000.};
        int itcount = 0;

        while (sqrt(sqr(delta[0])+sqr(delta[1])+sqr(delta[2])+sqr(delta[3])) 
               > tol*wx*wy && itcount < 20) {

            // Find zero

            p = F(ax,ay);
            px = Dx(ax,ay);
            py = Dy(ax,ay);
            delta[0] = -p*px/(px*px+py*py);
            delta[1] = -p*py/(px*px+py*py);
            tax = ax+delta[0];
            tay = ay+delta[1];

            // Find min distance

            p = F(tax,tay);
            px = Dx(tax,tay);
            py = Dy(tax,tay);
            double denom = px*px+py*py;
            delta[2] = py*(px*(tay-y)-py*(tax-x))/denom;
            delta[3] = -px*(px*(tay-y)-py*(tax-x))/denom;

            ax += delta[0]+delta[2];
            ay += delta[1]+delta[3];
            ++itcount;
        }

        clean = itcount < 20;
   
        return sqrt(sqr(x-ax)+sqr(y-ay));
    }
        
}
