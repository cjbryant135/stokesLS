#include "uniformmesh2d.h"

// Solve simple Poisson equation Del^2 u = f

#include <vector>
#include <map>
#include <mkl.h>

namespace levelset {
        
    class SparseMat 
    {
        std::vector<std::map<int,double> > ja;
        int count;
        int N;

        void CompressedRows(int*& i, int*& j, double*& a) const;
   
    public:
   
        SparseMat(const int n) : ja(n), count(0), N(n) {} 
   
        void Insert(const int i, const int j, const double x);

        std::vector<double> Solve(std::vector<double> b);
   
        friend std::vector<double> operator*(const SparseMat& M, const std::vector<double> x);
    };
    void SparseMat::Insert(const int r, const int c, const double x)
    {
        ja[r].insert(std::make_pair(c,x));
        ++count;
    }

    void SparseMat::CompressedRows(int*& i, int*& j, double*& a) const
    {
        i = new int[N+1];
        j = new int[count];
        a = new double[count];
        int k = 0;
        for (int r=0; r<N; ++r) {
            i[r] = r == 0 ? 1 : i[r-1]+ja[r-1].size();
            for (std::map<int,double>::const_iterator l=ja[r].begin(); 
                 l != ja[r].end(); ++l) {
                j[k] = l->first+1;
                a[k++] = l->second;
            }
        }
        i[N] = k+1;
    }

    extern "C" 
    {
//   void pardisoinit_(void* pt, int* mtype, int* iparm);
        void pardiso_(void* pt, int* maxfct, int* mnum, int* mtype, int* phase, 
                      int* n, double* a, int* ia, int* ja, int* perm, 
                      int* nrhs, int* iparm, int* msglvl, double* b, double* x, 
                      int* error);
    };


    std::vector<double> SparseMat::Solve(std::vector<double> b)
    {
        int* i;
        int* j;
        double* a;
   
        CompressedRows(i, j, a);

        void* pt[64];
        int mtype = 11;
        int iparm[64];
        int nrhs = 1;
        int maxfct, mnum, phase, error, msglvl, idum;
        double ddum;
   
        for (int k=0; k<64; ++k) {
            iparm[k] = 0;
            pt[k] = 0;
        }
   
        iparm[0] = 1;
        iparm[1] = 2;
        iparm[2] = 1;
        iparm[7] = 2;
        iparm[9] = 13;
        iparm[10] = 1;
        iparm[17] = -1;
        iparm[18] = -1;
        maxfct = 1;
        mnum = 1;
        msglvl = 1;
        error = 0;
   
        phase = 11;
        pardiso_(pt, &maxfct, &mnum, &mtype, &phase, &N, a, i, j, &idum, &nrhs,
                 iparm, &msglvl, &ddum, &ddum, &error);
   
        phase = 22;
        pardiso_(pt, &maxfct, &mnum, &mtype, &phase, &N, a, i, j, &idum, &nrhs,
                 iparm, &msglvl, &ddum, &ddum, &error);

        phase = 33;
        iparm[7] = 2;
        std::vector<double> x(N,0.);
        pardiso_(pt, &maxfct, &mnum, &mtype, &phase, &N, a, i, j, &idum, &nrhs,
                 iparm, &msglvl, &b[0], &x[0], &error);

        phase = -1;
        pardiso_(pt, &maxfct, &mnum, &mtype, &phase, &N, a, i, j, &idum, &nrhs,
                 iparm, &msglvl, &ddum, &ddum, &error);
   
        delete[] i;
        delete[] j;
        delete[] a;

        return x;
    }

    std::vector<double> operator*(const SparseMat& M, const std::vector<double> x)
    {
        std::vector<double> y(M.N,0.);
   
        for (int r=0; r<M.N; ++r) {
            for (std::map<int,double>::const_iterator i=M.ja[r].begin();
                 i != M.ja[r].end(); ++i) 
                y[r] += i->second * x[i->first];
        }
   
        return y;
    }

    extern "C" {
        void dgesv_(int* n, int* nrhs, double* A, int* lda, int* ipiv, double* b,
                    int* ldb, int* info);
    };

    void UniformMesh2D::IIM_coeffs(const double x0, const double y0, 
                                   const double nx, const double ny, const double chipp, 
                                   const double wx, const double wy, const int inside[6], double gamma[6])
    {
        int ii[6] = {-1,0,1,0,0,1};
        int jj[6] = {0,0,0,-1,1,1};
        double A[36] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
        double xi[6];
        double eta[6];
        double betam = 1.;
        double betap = 1.;
        double rho = betap/betam;
        
        for (int k=0; k<6; ++k) {
            double u=ii[k]*wx-x0;
            double v=jj[k]*wy-y0;
            eta[k] = nx*u+ny*v; 
            xi[k] = ny*u-nx*v;
        }
        
        // eqn 1: a1 + a2 = 0
        for (int j=0; j<6; ++j)
            A[j*6] = 1.;
        gamma[0] = 0.;
        
        // eqn 2: a3 +a4*rho+a8*(rho-1)*chi''+a10*(1-rho)*chi'' = 0
        for (int k=0; k<6; ++k) 
            A[1+k*6] = (inside[k] ? rho : 1.)*xi[k]
                +chipp*(1.-rho)
                *(inside[k] ? 0. : (eta[k]*eta[k]/2.-xi[k]*xi[k]/2.));
        gamma[1] = 0.;
        
        // eqn 3: a5+a6+a12*(1-rho)*chi'' = 0
        for (int k=0; k<6; ++k)
            A[2+k*6] = eta[k]+(1-rho)*chipp*(inside[k] ? 0. : eta[k]*xi[k]);
        gamma[2] = 0.;
                
        // eqn 4: a7+a8*rho = beta-
        for (int k=0; k<6; ++k)
            A[3+k*6] = xi[k]*xi[k]/2.*(inside[k] ? 1. : rho);
        gamma[3] = betam;
        
        // eqn 5: a9+a10+a8*(rho-1) = beta-
        for (int k=0; k<6; ++k)
            A[4+k*6] = eta[k]*eta[k]/2.+(inside[k] ? 0. : xi[k]*xi[k]/2.*(rho-1.));
        gamma[4] = betam;
        
        // eqn 6: a11 + a12*rho = 0;
        for (int k=0; k<6; ++k)
            A[5+k*6] = eta[k]*xi[k]*(inside[k] ? 1. : rho);
        gamma[5] = 0.;
        
        int ipiv[6];
        int n = 6;
        int nrhs = 1;
        int info;
        dgesv_(&n, &nrhs, A, &n, ipiv, gamma, &n, &info);
    }
        
                
        

    void UniformMesh2D::Poisson(const int kf, const int ku)
    {
        // First iteration:
        // kf contains boundary info on the boundary and the forcing function in
        // the interior.  Solution goes into ku.
        
        SparseMat M(maxi*maxj);
        std::vector<double> b(maxi*maxj,0.);
        
        for (int i=1; i<maxi-1; ++i) {
            M.Insert(i,i,1.);
            b[i] = data_(i,0,kf);
            M.Insert(i+(maxj-1)*maxi,i+(maxj-1)*maxi,1.);
            b[i+(maxj-1)*maxi] = data_(i,maxj-1,kf);
            for (int j=1; j<maxj-1; ++j) {
                M.Insert(i+j*maxi,i+j*maxi,-2./dx/dx-2./dy/dy);
                M.Insert(i+j*maxi,i-1+j*maxi,1/dx/dx);
                M.Insert(i+j*maxi,i+1+j*maxi,1/dx/dx);
                M.Insert(i+j*maxi,i+(j+1)*maxi,1/dy/dy);
                M.Insert(i+j*maxi,i+(j-1)*maxi,1/dy/dy);
                b[i+j*maxi] = data_(i,j,kf);
            }
        }
        
        for (int j=0; j<maxj; ++j) {
            M.Insert(j*maxi,j*maxi,1.);
            b[j*maxi] = data_(0,j,kf);
            M.Insert(maxi-1+j*maxi,maxi-1+j*maxi,1.);
            b[maxi-1+j*maxi] = data_(maxi-1,j,kf);
        }
        
        std::vector<double> x = M.Solve(b);
        
        for (int i=0; i<maxi; ++i)
            for (int j=0; j<maxj; ++j) 
                data_(i,j,ku) = x[i+j*maxi]; 
                        
    }

    double UniformMesh2D::NearestPoint(const double i, const double j, const int kphi,
                                       double& ax, double& ay, double& nx, double& ny,
                                       double& kappa)
    {
        Bicubic p;
        double dist = 10.*(dx+dy);
        double tdist, tax, tay;
        char clean;
        
        if (i>0 && j>0) {
            int score = (data_(i-1,j-1,kphi) < 0. ? 1 : 0)
                +(data_(i,j-1,kphi) < 0. ? 2 : 0)
                +(data_(i-1,j,kphi) < 0. ? 4 : 0)
                +(data_(i,j,phi) < 0. ? 8 : 0);
            if (score > 0 && score < 15) {              
                p.BuildwDeriv(data_(i-2,j-2,kphi), data_(i-1,j-2,kphi), data_(i,j-2,kphi)
                              data_(i+1,j-2,kphi), data_(i-2,j-1,kphi), data_(i-1,j-1,kphi),
                              data_(i,j-1,kphi), data_(i+1,j-1,kphi), data_(i-2,j,kphi),
                              data_(i-1,j,kphi), data_(i,j,kphi), data_(i+1,j,kphi),
                              data_(i-2,j+1,kphi), data_(i-1,j+1,kphi), data_(i,j+1,kphi),
                              data_(i+1,j+1,kphi), dx, dy);
                tdist = p.LocalDist(dx, dy, tax, tay, clean);
                if (clean && (tdist < dist)) {
                    ax = tax-dx;
                    ay = tay-dy;
                    nx = p.Dx(tax,tay);
                    ny = p.Dy(tax,tay);
                    kappa = (p.Dxx(tax,tay)*pow(p.Dy(tax,tay),2.)
                             +p.Dyy(tax,tay)*pow(p.Dx(tax,tay),2.)
                             -2.*p.Dxy(tax,tay)*p.Dx(tax,tay)*p.Dy(tax,tay))/
                        pow(pow(p.Dx(tax,tay),2.)+pow(p.Dy(tax,tay),2.),1.5);
                    dist = tdist;
                }
            }
        }
        if (i<maxi-1 && j>0) {
            int score = (data_(i,j-1,kphi) < 0. ? 1 : 0)
                +(data_(i+1,j-1,kphi) < 0. ? 2 : 0)
                +(data_(i,j,kphi) < 0. ? 4 : 0)
                +(data_(i+1,j,phi) < 0. ? 8 : 0);
            if (score > 0 && score < 15) {              
                p.BuildwDeriv(data_(i-1,j-2,kphi), data_(i,j-2,kphi), data_(i+1,j-2,kphi)
                              data_(i+2,j-2,kphi), data_(i-1,j-1,kphi), data_(i,j-1,kphi),
                              data_(i+1,j-1,kphi), data_(i+2,j-1,kphi), data_(i-1,j,kphi),
                              data_(i,j,kphi), data_(i+1,j,kphi), data_(i+2,j,kphi),
                              data_(i-1,j+1,kphi), data_(i,j+1,kphi), data_(i+1,j+1,kphi),
                              data_(i+2,j+1,kphi), dx, dy);
                tdist = p.LocalDist(dx, dy, tax, tay, clean);
                if (clean && (tdist < dist)) {
                    ax = tax;
                    ay = tay-dy;
                    nx = p.Dx(tax,tay);
                    ny = p.Dy(tax,tay);
                    kappa = (p.Dxx(tax,tay)*pow(p.Dy(tax,tay),2.)
                             +p.Dyy(tax,tay)*pow(p.Dx(tax,tay),2.)
                             -2.*p.Dxy(tax,tay)*p.Dx(tax,tay)*p.Dy(tax,tay))/
                        pow(pow(p.Dx(tax,tay),2.)+pow(p.Dy(tax,tay),2.),1.5);
                    dist = tdist;
                }
            }
        }
        if (i>0 && j<maxj-1) {
            int score = (data_(i-1,j,kphi) < 0. ? 1 : 0)
                +(data_(i,j,kphi) < 0. ? 2 : 0)
                +(data_(i-1,j+1,kphi) < 0. ? 4 : 0)
                +(data_(i,j+1,phi) < 0. ? 8 : 0);
            if (score > 0 && score < 15) {              
                p.BuildwDeriv(data_(i-2,j-1,kphi), data_(i-1,j-1,kphi), data_(i,j-1,kphi)
                              data_(i+1,j-1,kphi), data_(i-2,j,kphi), data_(i-1,j,kphi),
                              data_(i,j,kphi), data_(i+1,j,kphi), data_(i-2,j+1,kphi),
                              data_(i-1,j+1,kphi), data_(i,j+1,kphi), data_(i+1,j+1,kphi),
                              data_(i-2,j+2,kphi), data_(i-1,j+2,kphi), data_(i,j+2,kphi),
                              data_(i+1,j+2,kphi), dx, dy);
                tdist = p.LocalDist(dx, dy, tax, tay, clean);
                if (clean && (tdist < dist)) {
                    ax = tax-dx;
                    ay = tay;
                    nx = p.Dx(tax,tay);
                    ny = p.Dy(tax,tay);
                    kappa = (p.Dxx(tax,tay)*pow(p.Dy(tax,tay),2.)
                             +p.Dyy(tax,tay)*pow(p.Dx(tax,tay),2.)
                             -2.*p.Dxy(tax,tay)*p.Dx(tax,tay)*p.Dy(tax,tay))/
                        pow(pow(p.Dx(tax,tay),2.)+pow(p.Dy(tax,tay),2.),1.5);
                    dist = tdist;
                }
            }
        }
        if (i<maxi-1 && j<maxj-1) {
            int score = (data_(i,j,kphi) < 0. ? 1 : 0)
                +(data_(i+1,j,kphi) < 0. ? 2 : 0)
                +(data_(i,j+1,kphi) < 0. ? 4 : 0)
                +(data_(i+1,j+1,phi) < 0. ? 8 : 0);
            if (score > 0 && score < 15) {              
                p.BuildwDeriv(data_(i-1,j-1,kphi), data_(i,j-1,kphi), data_(i+1,j-1,kphi)
                              data_(i+2,j-1,kphi), data_(i-1,j,kphi), data_(i,j,kphi),
                              data_(i+1,j,kphi), data_(i+2,j,kphi), data_(i-1,j+1,kphi),
                              data_(i,j+1,kphi), data_(i+1,j+1,kphi), data_(i+2,j+1,kphi),
                              data_(i-1,j+2,kphi), data_(i,j+2,kphi), data_(i+1,j+2,kphi),
                              data_(i+2,j+2,kphi), dx, dy);
                tdist = p.LocalDist(dx, dy, tax, tay, clean);
                if (clean && (tdist < dist)) {
                    ax = tax;
                    ay = tay;
                    nx = p.Dx(tax,tay);
                    ny = p.Dy(tax,tay);
                    kappa = (p.Dxx(tax,tay)*pow(p.Dy(tax,tay),2.)
                             +p.Dyy(tax,tay)*pow(p.Dx(tax,tay),2.)
                             -2.*p.Dxy(tax,tay)*p.Dx(tax,tay)*p.Dy(tax,tay))/
                        pow(pow(p.Dx(tax,tay),2.)+pow(p.Dy(tax,tay),2.),1.5);
                    dist = tdist;
                }
            }
        }
        double denom = sqrt(nx*nx+ny*ny);
        nx /= denom;
        ny /= denom;

        return dist;
    }

    void UniformMesh2D::Poisson(const int kphi, const int kf, const int ku)
    {
        // First iteration:
        // kf contains boundary info on the boundary and the forcing function in
        // the interior.  Solution goes into ku.
        
        SparseMat M(maxi*maxj);
        std::vector<double> b(maxi*maxj,0.);
        
        for (int i=1; i<maxi-1; ++i) {
            M.Insert(i,i,1.);
            b[i] = data_(i,0,kf);
            M.Insert(i+(maxj-1)*maxi,i+(maxj-1)*maxi,1.);
            b[i+(maxj-1)*maxi] = data_(i,maxj-1,kf);
            for (int j=1; j<maxj-1; ++j) {
                int inside[6];
                inside[0] = data_(i-1,j,kphi) <= 0.;
                inside[1] = data_(i,j,kphi) <= 0.;
                inside[2] = data_(i+1,j,kphi) <= 0.;
                inside[3] = data_(i,j-1,kphi) <= 0.;
                inside[4] = data_(i,j+1,kphi) <= 0.;
                int total = inside[0]+2*inside[1]+4*inside[2]+8*inside[3]+16*inside[4];
                if (total > 0 && total < 31) {
                    double ax, ay, nx, ny, kappa;
                    double dist = NearestPoint(i,j,kphi,ax,ay,nx,ny,kappa);
                    ip = ax > 0. ? 1 : -1;
                    jp = ay > 0. ? 1 : -1;
                    inside[5] = data_(i+ip,j+jp,kphi);
                    double gamma[6];
                    IIM_coeffs(ax, ay, nx, ny, kappa,
                               dx, dy, inside, gamma);
                    M.Insert(i+j*maxi,i-1+j*maxi,gamma[0]);
                    M.Insert(i+j*maxi,i+j*maxi, gamma[1]);
                    M.Insert(i+j*maxi,i+1+j*maxi, gamma[2]);
                    M.Insert(i+j*maxi,i+(j-1)*maxi, gamma[3]);
                    M.Insert(i+j*maxi,i+(j+1)*maxi, gamma[4]);
                    M.insert(i+j*maxi,i+ip+(j+jp)*maxi, gamma[5]); // resume here but
                    // fix that the grid
                    // wasn't rotated.
                                
                    M.Insert(i+j*maxi,i+j*maxi,-2./dx/dx-2./dy/dy);
                    M.Insert(i+j*maxi,i-1+j*maxi,1/dx/dx);
                    M.Insert(i+j*maxi,i+1+j*maxi,1/dx/dx);
                    M.Insert(i+j*maxi,i+(j+1)*maxi,1/dy/dy);
                    M.Insert(i+j*maxi,i+(j-1)*maxi,1/dy/dy);
                    b[i+j*maxi] = data_(i,j,kf);
                }
            }
        
            for (int j=0; j<maxj; ++j) {
                M.Insert(j*maxi,j*maxi,1.);
                b[j*maxi] = data_(0,j,kf);
                M.Insert(maxi-1+j*maxi,maxi-1+j*maxi,1.);
                b[maxi-1+j*maxi] = data_(maxi-1,j,kf);
            }
        
            std::vector<double> x = M.Solve(b);
        
            for (int i=0; i<maxi; ++i)
                for (int j=0; j<maxj; ++j) 
                    data_(i,j,ku) = x[i+j*maxi]; 
                        
        }

    }
