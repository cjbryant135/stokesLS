//
//  pm2advect.cpp
//  LevelSet
//
//  Created by David Chopp on 7/25/19.
//  Copyright © 2019 David Chopp. All rights reserved.
//

#include "numtrait.h"
#include "datavec.h"
#include "polarmesh2d.h"
#include "pm2heap.h"
#include "pm2boundary.h"
#include "utility.h"
#include "genbicubic.h"

//DEBUG
#include<fstream>
namespace levelset {
    /*
		Original version
	 */
    void PolarMesh2D::ExtendVelocity(const int ktemp, const int kvel, const int ikus,
                                     const int ikhi, const int imask)
    {
        NumericTrait<double> t;
        
        // Initialize update status using mask
        
        int i, j;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                if (idata_(i,j,imask)) {
                    idata_(i,j,ikus) = 2;
                    data_(i,j,ktemp) = 0.;
                } else {
                    idata_(i,j,ikus) = 0;
                    data_(i,j,ktemp) = t.max_value;
                }
        
        // Now do fast marching to complete extension
        
        PolarMesh2D_Heap theHeap(this);
        
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                if (idata_(i,j,ikus) == 2) {
                    if (i > 0) {
                        UpdateNeighbor(theHeap, i-1, j, ktemp, ikus, ikhi);
                        if (fmm_coi[i]) {
                            if (j > 0)
                                UpdateNeighbor(theHeap, i-1, j-1, ktemp, ikus, ikhi);
                            else
                                UpdateNeighbor(theHeap, i-1, maxj-1, ktemp, ikus, ikhi);
                            if (j < maxj-1)
                                UpdateNeighbor(theHeap, i-1, j+1, ktemp, ikus, ikhi);
                            else
                                UpdateNeighbor(theHeap, i-1, 0, ktemp, ikus, ikhi);
                        }
                    }
                    if (i < maxi-1)
                        UpdateNeighbor(theHeap, i+1, j, ktemp, ikus, ikhi);
                    if (j > 0)
                        UpdateNeighbor(theHeap, i, j-1, ktemp, ikus, ikhi);
                    else
                        UpdateNeighbor(theHeap, i, maxj-1, ktemp, ikus, ikhi);
                    if (j < maxj-1)
                        UpdateNeighbor(theHeap, i, j+1, ktemp, ikus, ikhi);
                    else
                        UpdateNeighbor(theHeap, i, 0, ktemp, ikus, ikhi);
                }
        
        {
            int i, j;
            //         PlotWindow2D exdisp(true);
            //         exdisp.SetRange(-1., maxi, -1., maxj);
				while (!theHeap.Empty()) {
                theHeap.Extract(i,j,ikhi); 
					 //            exdisp.Dot(i,j);
                idata_(i,j,ikus) = 2;
                if (!idata_(i,j,imask))
                    ExtendVals(i, j, ktemp, ikus, ikhi, &kvel, 1);
                
                if (i > 0) {
                    UpdateNeighbor(theHeap, i-1, j, ktemp, ikus, ikhi);
                    if (fmm_coi[i]) {
                        if (j > 0)
                            UpdateNeighbor(theHeap, i-1, j-1, ktemp, ikus, ikhi);
                        else
                            UpdateNeighbor(theHeap, i-1, maxj-1, ktemp, ikus, ikhi);
                        if (j < maxj-1)
                            UpdateNeighbor(theHeap, i-1, j+1, ktemp, ikus, ikhi);
                        else
                            UpdateNeighbor(theHeap, i-1, 0, ktemp, ikus, ikhi);
                    }
                }
                if (i < maxi-1)
                    UpdateNeighbor(theHeap, i+1, j, ktemp, ikus, ikhi);
                if (j > 0)
                    UpdateNeighbor(theHeap, i, j-1, ktemp, ikus, ikhi);
                else
                    UpdateNeighbor(theHeap, i, maxj-1, ktemp, ikus, ikhi);
                if (j < maxj-1)
                    UpdateNeighbor(theHeap, i, j+1, ktemp, ikus, ikhi);
                else
                    UpdateNeighbor(theHeap, i, 0, ktemp, ikus, ikhi);
            }
        }

    }
    /*
		overloaded version that allows for velocity to be intialized at non-zero phi values
	 	Version used by the chimera code
	 */
	 void PolarMesh2D::ExtendVelocity(const int kphi, const int ktemp, const int kvel, const int ikus,
                                     const int ikhi, const int imask)
    {
        NumericTrait<double> t;
        
        // Initialize update status using mask
        
        int i, j;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                if (idata_(i,j,imask)) {
                    idata_(i,j,ikus) = 2;
                    //data_(i,j,ktemp) = 0.;
						  data_(i,j,ktemp) = data_(i,j,kphi);
                } else {
                    idata_(i,j,ikus) = 0;
                    data_(i,j,ktemp) = data_(i,j,kphi) >= 0. ? t.max_value : t.min_value;
                }
        
        // Now do fast marching to complete extension
        
        PolarMesh2D_Heap theHeap(this);
//		  std::cout << "Initializing the tentative set\n";        
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                if (idata_(i,j,ikus) == 2) {
                    if (i > 0) {
                        UpdateNeighbor(theHeap, i-1, j, ktemp, ikus, ikhi);
                        if (fmm_coi[i]) {
                            if (j > 0)
                                UpdateNeighbor(theHeap, i-1, j-1, ktemp, ikus, ikhi);
                            else
                                UpdateNeighbor(theHeap, i-1, maxj-1, ktemp, ikus, ikhi);
                            if (j < maxj-1)
                                UpdateNeighbor(theHeap, i-1, j+1, ktemp, ikus, ikhi);
                            else
                                UpdateNeighbor(theHeap, i-1, 0, ktemp, ikus, ikhi);
                        }
                    }
                    if (i < maxi-1)
                        UpdateNeighbor(theHeap, i+1, j, ktemp, ikus, ikhi);
                    if (j > 0)
                        UpdateNeighbor(theHeap, i, j-1, ktemp, ikus, ikhi);
                    else
                        UpdateNeighbor(theHeap, i, maxj-1, ktemp, ikus, ikhi);
                    if (j < maxj-1)
                        UpdateNeighbor(theHeap, i, j+1, ktemp, ikus, ikhi);
                    else
                        UpdateNeighbor(theHeap, i, 0, ktemp, ikus, ikhi);
                }
        
        {
            int i, j;
            //         PlotWindow2D exdisp(true);
            //         exdisp.SetRange(-1., maxi, -1., maxj);
//            std::cout << "Doing fast marching \n";
				//DEBUG
//				std::ofstream heapOut;
//				heapOut.open("heapOrder.csv");
				while (!theHeap.Empty()) {
                theHeap.Extract(i,j,ikhi);
                //DEBUG
//					 heapOut << i << "," << j << std::endl;
//					 if(fmm_coi[i]) {
//						std::cout << "Diagonal calculation being used\n";
//					 }
					 //            exdisp.Dot(i,j);
                idata_(i,j,ikus) = 2;
                if (!idata_(i,j,imask))
						  ExtendVals(i, j, ktemp, ikus, ikhi, &kvel, 1);
                if (i > 0) {
                    UpdateNeighbor(theHeap, i-1, j, ktemp, ikus, ikhi);
                    if (fmm_coi[i]) {
                        if (j > 0)
                            UpdateNeighbor(theHeap, i-1, j-1, ktemp, ikus, ikhi);
                        else
                            UpdateNeighbor(theHeap, i-1, maxj-1, ktemp, ikus, ikhi);
                        if (j < maxj-1)
                            UpdateNeighbor(theHeap, i-1, j+1, ktemp, ikus, ikhi);
                        else
                            UpdateNeighbor(theHeap, i-1, 0, ktemp, ikus, ikhi);
                    }
                }
                if (i < maxi-1)
                    UpdateNeighbor(theHeap, i+1, j, ktemp, ikus, ikhi);
                if (j > 0)
                    UpdateNeighbor(theHeap, i, j-1, ktemp, ikus, ikhi);
                else
                    UpdateNeighbor(theHeap, i, maxj-1, ktemp, ikus, ikhi);
                if (j < maxj-1)
                    UpdateNeighbor(theHeap, i, j+1, ktemp, ikus, ikhi);
                else
                    UpdateNeighbor(theHeap, i, 0, ktemp, ikus, ikhi);
            }
        	//DEBUG
			//heapOut.close();
		  }

    }
    
    double PolarMesh2D::ComputeTVal(const int i, const int j, const int k,
                                      const int ikus)
    {
        int sr[2];
        int sth[2][2];

		  int jM, jP; 
		  jM = j > 0 ? j-1 : maxj-1;
		  jP = j < maxj-1 ? j+1 : 0;
        
        sr[0] = idata_(i-1,j,ikus) == 2 ? 1 : 0;
        sr[1] = idata_(i+1,j,ikus) == 2 ? 1 : 0;
        sth[0][0] = idata_(i,jM,ikus) == 2 ? 1 : 0;
        sth[0][1] = idata_(i,jP,ikus) == 2 ? 1 : 0;
        sth[1][0] = (idata_(i-1,jM,ikus) == 2 && fmm_coi[i]) ? 1 : 0;
        sth[1][1] = (idata_(i-1,jP,ikus) == 2 && fmm_coi[i]) ? 1 : 0;
        
        double temptval;
        double tval = DBL_MAX;
        int sgn = data_(i,j,k) >= 0. ? 1 : -1;
        
        // In clockwise
        
        if (fmm_coi[i]) {
            if (sr[0] && sth[1][0]) {
                temptval = fmm_coeff[0][0][i]*data_(i-1,j,k) + fmm_coeff[0][1][i]*data_(i-1,jM,k)
                            +sgn*fmm_coeff[0][2][i]*sqrt(fmm_coeff[0][3][i]-sqr(data_(i-1,j,k)-data_(i-1,jM,k)));
                if (sgn*temptval > sgn*(fmm_val[0][0][i]*data_(i-1,j,k)+fmm_val[0][1][i]*data_(i-1,jM,k))
                    && sgn*temptval > sgn*(fmm_val[0][2][i]*data_(i-1,j,k)+fmm_val[0][3][i]*data_(i-1,jM,k))
                    && abs(temptval) < abs(tval))
                    tval = temptval;
            }
            
            if (sth[1][0] && sth[0][0]) {
                temptval = fmm_coeff[1][0][i]*data_(i-1,jM,k)+fmm_coeff[1][1][i]*data_(i,jM,k)
                            +sgn*fmm_coeff[1][2][i]*sqrt(fmm_coeff[1][3][i]-sqr(data_(i-1,jM,k)-data_(i,jM,k)));
                if (sgn*temptval > sgn*(fmm_val[1][0][i]*data_(i-1,jM,k)+fmm_val[1][1][i]*data_(i,jM,k))
                    && sgn*temptval > sgn*(fmm_val[1][2][i]*data_(i-1,jM,k)+fmm_val[1][3][i]*data_(i,jM,k))
                    && abs(temptval) < abs(tval))
                    tval = temptval;
            }

            if (sr[0] && sth[1][1]) {
                temptval = fmm_coeff[0][0][i]*data_(i-1,j,k) + fmm_coeff[0][1][i]*data_(i-1,jP,k)
                +sgn*fmm_coeff[0][2][i]*sqrt(fmm_coeff[0][3][i]-sqr(data_(i-1,j,k)-data_(i-1,jP,k)));
                if (sgn*temptval > sgn*(fmm_val[0][0][i]*data_(i-1,j,k)+fmm_val[0][1][i]*data_(i-1,jP,k))
                    && sgn*temptval > sgn*(fmm_val[0][2][i]*data_(i-1,j,k)+fmm_val[0][3][i]*data_(i-1,jP,k))
                    && abs(temptval) < abs(tval))
                    tval = temptval;
            }
            
            if (sth[1][0] && sth[0][1]) {
                temptval = fmm_coeff[1][0][i]*data_(i-1,jP,k)+fmm_coeff[1][1][i]*data_(i,jP,k)
                +sgn*fmm_coeff[1][2][i]*sqrt(fmm_coeff[1][3][i]-sqr(data_(i-1,jP,k)-data_(i,jP,k)));
                if (sgn*temptval > sgn*(fmm_val[1][0][i]*data_(i-1,jP,k)+fmm_val[1][1][i]*data_(i,jP,k))
                    && sgn*temptval > sgn*(fmm_val[1][2][i]*data_(i-1,jP,k)+fmm_val[1][3][i]*data_(i,jP,k))
                    && abs(temptval) < abs(tval))
                    tval = temptval;
            }
        } else {
            if (sr[0] && sth[0][0]) {
                temptval = fmm_coeff[2][0][i]*data_(i-1,j,k) + fmm_coeff[2][1][i]*data_(i,jM,k)
                +sgn*fmm_coeff[2][2][i]*sqrt(fmm_coeff[2][3][i]-sqr(data_(i-1,j,k)-data_(i,jM,k)));
                if (sgn*temptval > sgn*(fmm_val[2][0][i]*data_(i-1,j,k)+fmm_val[2][1][i]*data_(i,jM,k))
                    && sgn*temptval > sgn*(fmm_val[2][2][i]*data_(i-1,j,k)+fmm_val[2][3][i]*data_(i,jM,k))
                    && abs(temptval) < abs(tval))
                    tval = temptval;
            }

            if (sr[0] && sth[0][1]) {
                temptval = fmm_coeff[2][0][i]*data_(i-1,j,k) + fmm_coeff[2][1][i]*data_(i,jP,k)
                +sgn*fmm_coeff[2][2][i]*sqrt(fmm_coeff[2][3][i]-sqr(data_(i-1,j,k)-data_(i,jP,k)));
                if (sgn*temptval > sgn*(fmm_val[2][0][i]*data_(i-1,j,k)+fmm_val[2][1][i]*data_(i,jP,k))
                    && sgn*temptval > sgn*(fmm_val[2][2][i]*data_(i-1,j,k)+fmm_val[2][3][i]*data_(i,jP,k))
                    && abs(temptval) < abs(tval))
                    tval = temptval;
            }
        }
        
        if (sr[1] && sth[0][0]) {
            temptval = fmm_coeff[3][0][i]*data_(i,jM,k) + fmm_coeff[3][1][i]*data_(i+1,j,k)
            +sgn*fmm_coeff[3][2][i]*sqrt(fmm_coeff[3][3][i]-sqr(data_(i,jM,k)-data_(i+1,j,k)));
            if (sgn*temptval > sgn*(fmm_val[3][0][i]*data_(i,jM,k)+fmm_val[3][1][i]*data_(i+1,j,k))
                && sgn*temptval > sgn*(fmm_val[3][2][i]*data_(i,jM,k)+fmm_val[3][3][i]*data_(i+1,j,k))
                && abs(temptval) < abs(tval))
                tval = temptval;
        }

        if (sr[1] && sth[0][1]) {
            temptval = fmm_coeff[3][0][i]*data_(i,jP,k) + fmm_coeff[3][1][i]*data_(i+1,j,k)
            +sgn*fmm_coeff[3][2][i]*sqrt(fmm_coeff[3][3][i]-sqr(data_(i,jP,k)-data_(i+1,j,k)));
            if (sgn*temptval > sgn*(fmm_val[3][0][i]*data_(i,jP,k)+fmm_val[3][1][i]*data_(i+1,j,k))
                && sgn*temptval > sgn*(fmm_val[3][2][i]*data_(i,jP,k)+fmm_val[3][3][i]*data_(i+1,j,k))
                && abs(temptval) < abs(tval))
                tval = temptval;
        }

        if (tval == DBL_MAX) {
            if (sr[0]) {
                temptval = data_(i-1,j,k)+sgn*dr[i-1];
                if (abs(temptval) < abs(tval))
                    tval = temptval;
            }
            
            if (sr[1]) {
                temptval = data_(i+1,j,k)+sgn*dr[i];
                if (abs(temptval) < abs(tval))
                    tval = temptval;
            }
            
            if (sth[0][0]) {
                temptval = data_(i,jM,k)+sgn*r[i]*sqrt(2*(1-cos(dtheta)));
                if (abs(temptval) < abs(tval))
                    tval = temptval;
            }
            
            if (sth[0][1]) {
                temptval = data_(i,jP,k)+sgn*r[i]*sqrt(2*(1-cos(dtheta)));
                if (abs(temptval) < abs(tval))
                    tval = temptval;
            }
        }
       	
//		  if(tval < 0.) {
//			std::cout << "Negative tval computed\n";
//		  }
		  
        return data_(i,j,k) = tval;
    }

    
    void PolarMesh2D::UpdateNeighbor(PolarMesh2D_Heap& heap, const int i,
                                       const int j, const int k, const int ikus,
                                       const int ikhi)
    {
        switch(idata_(i,j,ikus)) {
            case 0:
                heap.Insert(i,j,fabs(ComputeTVal(i,j,k,ikus)),ikhi);
                idata_(i,j,ikus) = 1;
                break;
            case 1:
                heap.Change(idata_(i,j,ikhi),fabs(ComputeTVal(i,j,k,ikus)),ikhi);
                break;
            case 2:
                break;
            default:
                std::cerr << "Error: Gridpoint " << i << "," << j << " not initialized.\n";
                break;
        }
    }
    
    void PolarMesh2D::ExtendVals(const int i, const int j, const int k,
                                   const int ikus, const int ikhi,
                                   const int* p, const int pl)
    {
        int sr[2];
        int sth[2][2];
        char done = false;
       	
		  int jM, jP;
		  jM = j > 0 ? j-1 : maxj-1;
		  jP = j < maxj-1 ? j+1 : 0;
		  
		  //is i-1 neighbor accepted?
        sr[0] = (i > 0 && idata_(i-1,j,ikus) == 2) ? 1 : 0;
        //is i+1 neighbor accepted
		  sr[1] = (i < maxi-1 && idata_(i+1,j,ikus) == 2) ? 1 : 0;
        //is j-1 neighbor accepted
		  sth[0][0] = idata_(i,jM,ikus) == 2 ? 1 : 0;
        //is j+1 neighbor accepted
		  sth[0][1] = idata_(i,jP,ikus) == 2 ? 1 : 0;
        sth[1][0] = (i > 0 && idata_(i-1,jM,ikus) == 2 && fmm_coi[i]) ? 1 : 0;
        sth[1][1] = (i > 0 && idata_(i-1,jP,ikus) == 2 && fmm_coi[i]) ? 1 : 0;

        if (fmm_coi[i]) {
            if (sr[0] && sth[1][0]
                && fabs(data_(i,j,k)) > fabs(fmm_val[0][0][i]*data_(i-1,j,k)+fmm_val[0][1][i]*data_(i-1,jM,k))
                && fabs(data_(i,j,k)) > fabs(fmm_val[0][2][i]*data_(i-1,j,k)+fmm_val[0][3][i]*data_(i-1,jM,k))) {
                for (int l=0; l<pl; ++l)
                    data_(i,j,p[l]) = (data_(i-1,j,p[l])*(fmm_ext[0][1][i]*(data_(i,j,k)-data_(i-1,j,k))
                                                         +fmm_ext[0][2][i]*(data_(i,j,k)-data_(i-1,jM,k)))
                                       +data_(i-1,jM,p[l])*(fmm_ext[0][2][i]*(data_(i,j,k)-data_(i-1,j,k))
                                                             +fmm_ext[0][0][i]*(data_(i,j,k)-data_(i-1,jM,k))))
                                    / ((fmm_ext[0][1][i]+fmm_ext[0][2][i])*(data_(i,j,k)-data_(i-1,j,k))
                                       +(fmm_ext[0][0][i]+fmm_ext[0][2][i])*(data_(i,j,k)-data_(i-1,jM,k)));
                done = true;
            }
            
            if (sth[1][0] && sth[0][0]
                && fabs(data_(i,j,k)) > fabs(fmm_val[1][0][i]*data_(i-1,jM,k)+fmm_val[1][1][i]*data_(i,jM,k))
                && fabs(data_(i,j,k)) > fabs(fmm_val[1][2][i]*data_(i-1,jM,k)+fmm_val[1][3][i]*data_(i,jM,k))) {
                for (int l=0; l<pl; ++l)
                    data_(i,j,p[l]) = (data_(i-1,jM,p[l])*(fmm_ext[1][1][i]*(data_(i,j,k)-data_(i-1,jM,k))
                                                          +fmm_ext[1][2][i]*(data_(i,j,k)-data_(i,jM,k)))
                                       +data_(i,jM,p[l])*(fmm_ext[1][2][i]*(data_(i,j,k)-data_(i-1,jM,k))
                                                             +fmm_ext[1][0][i]*(data_(i,j,k)-data_(i,jM,k))))
                                    / ((fmm_ext[1][1][i]+fmm_ext[1][2][i])*(data_(i,j,k)-data_(i-1,jM,k))
                                       +(fmm_ext[1][0][i]+fmm_ext[1][2][i])*(data_(i,j,k)-data_(i,jM,k)));
                done = true;
            }
            
            if (sr[0] && sth[1][1]
                && fabs(data_(i,j,k)) > fabs(fmm_val[0][0][i]*data_(i-1,j,k)+fmm_val[0][1][i]*data_(i-1,jP,k))
                && fabs(data_(i,j,k)) > fabs(fmm_val[0][2][i]*data_(i-1,j,k)+fmm_val[0][3][i]*data_(i-1,jP,k))) {
                    for (int l=0; l<pl; ++l)
                        data_(i,j,p[l]) = (data_(i-1,j,p[l])*(fmm_ext[0][1][i]*(data_(i,j,k)-data_(i-1,j,k))
                                                                +fmm_ext[0][2][i]*(data_(i,j,k)-data_(i-1,jP,k)))
                                           +data_(i-1,jP,p[l])*(fmm_ext[0][2][i]*(data_(i,j,k)-data_(i-1,j,k))
                                                               +fmm_ext[0][0][i]*(data_(i,j,k)-data_(i-1,jP,k))))
                                        / ((fmm_ext[0][1][i]+fmm_ext[0][2][i])*(data_(i,j,k)-data_(i-1,j,k))
                                           +(fmm_ext[0][0][i]+fmm_ext[0][2][i])*(data_(i,j,k)-data_(i-1,jP,k)));
                done = true;
            }
            
            if (sth[1][0] && sth[0][1]
                && fabs(data_(i,j,k)) > fabs(fmm_val[1][0][i]*data_(i-1,jP,k)+fmm_val[1][1][i]*data_(i,jP,k))
                && fabs(data_(i,j,k)) > fabs(fmm_val[1][2][i]*data_(i-1,jP,k)+fmm_val[1][3][i]*data_(i,jP,k))) {
                for (int l=0; l<pl; ++l)
                    data_(i,j,p[l]) = (data_(i-1,jP,p[l])*(fmm_ext[1][1][i]*(data_(i,j,k)-data_(i-1,jP,k))
                                                          +fmm_ext[1][2][i]*(data_(i,j,k)-data_(i,jP,k)))
                                       +data_(i,jP,p[l])*(fmm_ext[1][2][i]*(data_(i,j,k)-data_(i-1,jP,k))
                                                           +fmm_ext[1][0][i]*(data_(i,j,k)-data_(i,jP,k))))
                                        / ((fmm_ext[1][1][i]+fmm_ext[1][2][i])*(data_(i,j,k)-data_(i-1,jP,k))
                                           +(fmm_ext[1][0][i]+fmm_ext[1][2][i])*(data_(i,j,k)-data_(i,jP,k)));
                done = true;
            }
            
        } else {
            
            if (sr[0] && sth[0][0]
                && (fabs(data_(i,j,k)) > fabs(fmm_val[2][0][i]*data_(i-1,j,k)+fmm_val[2][1][i]*data_(i,jM,k)) ||
						fabs(data_(i,j,k) - (fmm_val[2][0][i]*data_(i-1,j,k)+fmm_val[2][1][i]*data_(i,jM,k))) < 1.e-12
					 )
                && (fabs(data_(i,j,k)) > fabs(fmm_val[2][2][i]*data_(i-1,j,k)+fmm_val[2][3][i]*data_(i,jM,k)) ||
					   fabs(data_(i,j,k) - (fmm_val[2][2][i]*data_(i-1,j,k)+fmm_val[2][3][i]*data_(i,jM,k))) < 1.e-12
					 )
					 ) {
                for (int l=0; l<pl; ++l)
                    data_(i,j,p[l]) = (data_(i-1,j,p[l])*(fmm_ext[2][1][i]*(data_(i,j,k)-data_(i-1,j,k))
                                                          +fmm_ext[2][2][i]*(data_(i,j,k)-data_(i,jM,k)))
                                       +data_(i,jM,p[l])*(fmm_ext[2][2][i]*(data_(i,j,k)-data_(i-1,j,k))
                                                           +fmm_ext[2][0][i]*(data_(i,j,k)-data_(i,jM,k))))
                    / ((fmm_ext[2][1][i]+fmm_ext[2][2][i])*(data_(i,j,k)-data_(i-1,j,k))
                       +(fmm_ext[2][0][i]+fmm_ext[2][2][i])*(data_(i,j,k)-data_(i,jM,k)));
                done = true;
            }
            
            if (sr[0] && sth[0][1]
                && (fabs(data_(i,j,k)) > fabs(fmm_val[2][0][i]*data_(i-1,j,k)+fmm_val[2][1][i]*data_(i,jP,k)) ||
                	fabs(data_(i,j,k) - (fmm_val[2][0][i]*data_(i-1,j,k)+fmm_val[2][1][i]*data_(i,jP,k))) < 1.e-12 )
					 && (fabs(data_(i,j,k)) > fabs(fmm_val[2][2][i]*data_(i-1,j,k)+fmm_val[2][3][i]*data_(i,jP,k)) ||
					   fabs(data_(i,j,k) - (fmm_val[2][2][i]*data_(i-1,j,k)+fmm_val[2][3][i]*data_(i,jP,k)) ) < 1.e-12 )
					 ) {
                
					 for (int l=0; l<pl; ++l)
                    data_(i,j,p[l]) = (data_(i-1,j,p[l])*(fmm_ext[2][1][i]*(data_(i,j,k)-data_(i-1,j,k))
                                                          +fmm_ext[2][2][i]*(data_(i,j,k)-data_(i,jP,k)))
                                       +data_(i,jP,p[l])*(fmm_ext[2][2][i]*(data_(i,j,k)-data_(i-1,j,k))
                                                           +fmm_ext[2][0][i]*(data_(i,j,k)-data_(i,jP,k))))
                    / ((fmm_ext[2][1][i]+fmm_ext[2][2][i])*(data_(i,j,k)-data_(i-1,j,k))
                       +(fmm_ext[2][0][i]+fmm_ext[2][2][i])*(data_(i,j,k)-data_(i,jP,k)));
                done = true;
            }
        }
        
        if (sr[1] && sth[0][0]
            && (fabs(data_(i,j,k)) > fabs(fmm_val[3][0][i]*data_(i,jM,k)+fmm_val[3][1][i]*data_(i+1,j,k)) ||
					fabs(data_(i,j,k) - (fmm_val[3][0][i]*data_(i,jM,k)+fmm_val[3][1][i]*data_(i+1,j,k)) ) < 1.e-12
					) 
				&& (fabs(data_(i,j,k)) > fabs(fmm_val[3][2][i]*data_(i,jM,k)+fmm_val[3][3][i]*data_(i+1,j,k)) ||
				  fabs(data_(i,j,k) - (fmm_val[3][2][i]*data_(i,jM,k)+fmm_val[3][3][i]*data_(i+1,j,k)) ) < 1.e-12
				  )
				) {
				for (int l=0; l<pl; ++l)
                data_(i,j,p[l]) = (data_(i,jM,p[l])*(fmm_ext[3][1][i]*(data_(i,j,k)-data_(i,jM,k))
                                                      +fmm_ext[3][2][i]*(data_(i,j,k)-data_(i+1,j,k)))
                                   +data_(i+1,j,p[l])*(fmm_ext[3][2][i]*(data_(i,j,k)-data_(i,jM,k))
                                                       +fmm_ext[3][0][i]*(data_(i,j,k)-data_(i+1,j,k))))
                / ((fmm_ext[3][1][i]+fmm_ext[3][2][i])*(data_(i,j,k)-data_(i,jM,k))
                   +(fmm_ext[3][0][i]+fmm_ext[3][2][i])*(data_(i,j,k)-data_(i+1,j,k)));
            done = true;
        }
        
        if (sr[1] && sth[0][1]
            && (fabs(data_(i,j,k)) > fabs(fmm_val[3][0][i]*data_(i,jP,k)+fmm_val[3][1][i]*data_(i+1,j,k)) ||
            	fabs(data_(i,j,k) - (fmm_val[3][0][i]*data_(i,jP,k)+fmm_val[3][1][i]*data_(i+1,j,k))) < 1.e-12
				)
				&& (fabs(data_(i,j,k)) > fabs(fmm_val[3][2][i]*data_(i,jP,k)+fmm_val[3][3][i]*data_(i+1,j,k)) ||
					fabs(data_(i,j,k) - (fmm_val[3][2][i]*data_(i,jP,k)+fmm_val[3][3][i]*data_(i+1,j,k))) < 1.e-12
				)
				) {
				for (int l=0; l<pl; ++l)
                data_(i,j,p[l]) = (data_(i,jP,p[l])*(fmm_ext[3][1][i]*(data_(i,j,k)-data_(i,jP,k))
                                                      +fmm_ext[3][2][i]*(data_(i,j,k)-data_(i+1,j,k)))
                                   +data_(i+1,j,p[l])*(fmm_ext[3][2][i]*(data_(i,j,k)-data_(i,jP,k))
                                                       +fmm_ext[3][0][i]*(data_(i,j,k)-data_(i+1,j,k))))
                / ((fmm_ext[3][1][i]+fmm_ext[3][2][i])*(data_(i,j,k)-data_(i,jP,k))
                   +(fmm_ext[3][0][i]+fmm_ext[3][2][i])*(data_(i,j,k)-data_(i+1,j,k)));
            done = true;
        }
        
        if (!done) {
            if (sr[0]) {
                for (int l=0; l<pl; ++l)
                    data_(i,j,p[l]) = data_(i-1,j,p[l]);
            }
            
            if (sr[1]) {
                for (int l=0; l<pl; ++l)
                    data_(i,j,p[l]) = data_(i+1,j,p[l]);
            }
            
            if (sth[0][0]) {
                for (int l=0; l<pl; ++l)
                    data_(i,j,p[l]) = data_(i,jM,p[l]);
            }
            
            if (sth[0][1]) {
                for (int l=0; l<pl; ++l)
                    data_(i,j,p[l]) = data_(i,jP,p[l]);
            }
        }
    }
    
    void PolarMesh2D::Reinitialize(const int func, const int phi, const int ikus,
                                     const int ikhi, const int dir)
    {
        NumericTrait<double> t;
        
        // Initialize the accepted and tentative points
        
        // Initialize the ikus
        
        int i, j;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j) {
                idata_(i,j,ikus) = 0;
                idata_(i,j,ikhi) = -1;
                data_(i,j,phi) = data_(i,j,func) >= 0. ? t.max_value : t.min_value;
            }
        
        int cross;
        for (i=0; i<maxi-1; ++i)
            for (j=0; j<maxj-1; ++j) {
                cross = data_(i,j,func) > 0. ? 1 : 0;
                cross += data_(i+1,j,func) > 0 ? 2 : 0;
                cross += data_(i,j+1,func) > 0 ? 4 : 0;
                cross += data_(i+1,j+1,func) > 0 ? 8 : 0;
                if (cross > 0 && cross < 15)
                    InitLevelSet(i,j,func,phi,ikus);
            }
        
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                if (idata_(i,j,ikus) == 1) idata_(i,j,ikus) = 0;
        
        // Now do fast marching to complete extension
        
        PolarMesh2D_Heap theHeap(this);
        
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j) {
                if (idata_(i,j,ikus)==2) {
                    if (i > 0) UpdateNeighbor(theHeap, i-1, j, phi, ikus, ikhi);
                    if (j > 0) UpdateNeighbor(theHeap, i, j-1, phi, ikus, ikhi);
                    if (i < maxi-1) UpdateNeighbor(theHeap, i+1, j, phi, ikus, ikhi);
                    if (j < maxj-1) UpdateNeighbor(theHeap, i, j+1, phi, ikus, ikhi);
                }
            }
        
        {
            int i, j;
            while (!theHeap.Empty()) {
                theHeap.Extract(i,j,ikhi);
                idata_(i,j,ikus) = 2;
                if (i > 0) UpdateNeighbor(theHeap, i-1, j, phi, ikus, ikhi);
                if (j > 0) UpdateNeighbor(theHeap, i, j-1, phi, ikus, ikhi);
                if (i < maxi-1) UpdateNeighbor(theHeap, i+1, j, phi, ikus, ikhi);
                if (j < maxj-1) UpdateNeighbor(theHeap, i, j+1, phi, ikus, ikhi);
            }
        }
        
        switch (dir) {
            case 0:
                for (i=0; i<maxi; ++i)
                    for (j=0; j<maxj; ++j)
                        data_(i,j,func) = data_(i,j,phi);
                break;
            case 1:
                for (i=0; i<maxi; ++i)
                    for (j=0; j<maxj; ++j)
                        if (data_(i,j,phi) > 0) data_(i,j,func) = data_(i,j,phi);
                break;
            case -1:
                for (i=0; i<maxi; ++i)
                    for (j=0; j<maxj; ++j)
                        if (data_(i,j,phi) < 0) data_(i,j,func) = data_(i,j,phi);
                break;
        }
        
        bc->Apply(func);
        
    }

    void PolarMesh2D::InitLevelSet(const int i, const int j, const int func,
                                     const int phi, const int ikus)
    {
        double c[6] = {dr[i-1]/dr[i]/(dr[i]+dr[i-1]),
                        (dr[i]-dr[i-1])/dr[i]/dr[i-1],
                        -dr[i]/(dr[i-1]*(dr[i-1]+dr[i])),
                        dr[i]/dr[i+1]/(dr[i+1]+dr[i]),
                        (dr[i+1]-dr[i])/dr[i+1]/dr[i],
                        -dr[i+1]/(dr[i]*(dr[i]+dr[i+1]))
        };
        GeneralBicubic p(data_(i,j,func), data_(i+1,j,func), data_(i,j+1,func), data_(i+1,j+1,func),
                  c[0]*data_(i+1,j,func)+c[1]*data_(i,j,func)+c[2]*data_(i-1,j,func),
                  c[3]*data_(i+2,j,func)+c[4]*data_(i+1,j,func)+c[5]*data_(i,j,func),
                  c[0]*data_(i+1,j+1,func)+c[1]*data_(i,j+1,func)+c[2]*data_(i-1,j+1,func),
                  c[3]*data_(i+2,j+1,func)+c[4]*data_(i+1,j+1,func)+c[5]*data_(i,j+1,func),
                  (data_(i,j+1,func)-data_(i,j-1,func))/2/dtheta,
                  (data_(i+1,j+1,func)-data_(i+1,j-1,func))/2/dtheta,
                  (data_(i,j+2,func)-data_(i,j,func))/2/dtheta,
                  (data_(i+1,j+2,func)-data_(i+1,j,func))/2/dtheta,
                  ((c[0]*data_(i+1,j+1,func)+c[1]*data_(i,j+1,func)+c[2]*data_(i-1,j+1,func))
                    -(c[0]*data_(i+1,j-1,func)+c[1]*data_(i,j-1,func)+c[2]*data_(i-1,j-1,func)))/2/dtheta,
                  ((c[3]*data_(i+1,j+1,func)+c[4]*data_(i,j+1,func)+c[5]*data_(i-1,j+1,func))
                   -(c[3]*data_(i+1,j-1,func)+c[4]*data_(i,j-1,func)+c[5]*data_(i-1,j-1,func)))/2/dtheta,
                  ((c[0]*data_(i+1,j+2,func)+c[1]*data_(i,j+2,func)+c[2]*data_(i-1,j+2,func))
                   -(c[0]*data_(i+1,j,func)+c[1]*data_(i,j,func)+c[2]*data_(i-1,j,func)))/2/dtheta,
                  ((c[3]*data_(i+1,j+2,func)+c[4]*data_(i,j+2,func)+c[5]*data_(i-1,j+2,func))
                   -(c[3]*data_(i+1,j,func)+c[4]*data_(i,j,func)+c[5]*data_(i-1,j,func)))/2/dtheta,
                  dr[i-1], dr[i], dr[i+1], dtheta, dtheta, dtheta);
        
        double ar, atheta;
        double dist[2][2];
        char clean[2][2];
        dist[0][0] = p.LocalDistPolar(r[i], 0., 0., ar, atheta, clean[0][0]);
        dist[1][0] = p.LocalDistPolar(r[i], dr[i], 0., ar, atheta, clean[1][0]);
        dist[0][1] = p.LocalDistPolar(r[i], 0., dtheta, ar, atheta, clean[0][1]);
        dist[1][1] = p.LocalDistPolar(r[i], dr[i], dtheta, ar, atheta, clean[1][1]);
        
        for (int ii=0; ii<2; ++ii)
            for (int jj=0; jj<2; ++jj) {
                if (clean[ii][jj]) {
                    if (idata_(i+ii,j+jj,ikus) < 2) {
                        data_(i+ii,j+jj,phi) = copysign(dist[ii][jj],
                                                        data_(i+ii,j+jj,func));
                        idata_(i+ii,j+jj,ikus) = 2;
                    }
                    else if (abs(data_(i+ii,j+jj,phi)) > dist[ii][jj]) {
                        data_(i+ii,j+jj,phi) =copysign(dist[ii][jj],
                                                       data_(i+ii,j+jj,func));
                        idata_(i+ii,j+jj,ikus) = 2;
                    }
                } else if (abs(data_(i+ii,j+jj,phi)) > dist[ii][jj]) {
                    data_(i+ii,j+jj,phi) =copysign(dist[ii][jj],data_(i+ii,j+jj,func));
                    idata_(i+ii,j+jj,ikus) = 1;
                }
            }
    }
    
    void PolarMesh2D::Advance(const int k, const int kv, const int knorm, const double dt)
    {
        int i, j;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                data_(i,j,k) += dt*data_(i,j,knorm)*data_(i,j,kv);
        bc->Apply(k);
    }
    
    void PolarMesh2D::Advance(const int k, const int kv, const double dt)
    {
        int i, j;
        for (i=0; i<maxi; ++i)
            for (j=0; j<maxj; ++j)
                data_(i,j,k) += dt*data_(i,j,kv);
        bc->Apply(k);
    }

}
