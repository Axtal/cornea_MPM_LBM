/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2020 Sergio Galindo                                    *
 * Copyright (C) 2020 Pei Zhang                                         *
 * Copyright (C) 2020 Siqi Sun                                          *
 *                                                                      *
 * This program is free software: you can redistribute it and/or modify *
 * it under the terms of the GNU General Public License as published by *
 * the Free Software Foundation, either version 3 of the License, or    *
 * any later version.                                                   *
 *                                                                      *
 * This program is distributed in the hope that it will be useful,      *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         *
 * GNU General Public License for more details.                         *
 *                                                                      *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/
// Fish tail loaded from mesh


// MechSys
#include <mechsys/lbmmpm/Domain.h>
#include <math.h>
#include <iostream>
#include <fstream>

struct UserData
{
    double Tf; //final time
    double rhoa; //density of air
    double umax; //maximun velocity of the outlet
    double Ro;   //Radius of the outlet
};

void Setup (LBMMPM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    size_t nx = dom.LBMDOM.Ndim(0);
    size_t ny = dom.LBMDOM.Ndim(1);
    size_t nz = dom.LBMDOM.Ndim(2);
    // x- boundary
    #pragma omp parallel for schedule(static) num_threads(dom.LBMDOM.Nproc)
    for (size_t i=0; i<ny; ++i) //y-axis
    for (size_t j=0; j<nz; ++j) //z-axis
    {
        double * f = dom.LBMDOM.F[0][0][i][j];
        double ux  = -dom.LBMDOM.Cs*(-dat.rhoa+f[0]+f[3]+f[4]+f[5]+f[6]+2.0*(f[2]+f[8]+f[10]+f[12]+f[14]))/dat.rhoa;
        f[ 1] = f[ 2] +  (2.0/3.0)*dat.rhoa*ux/dom.LBMDOM.Cs;
        f[ 7] = f[ 8] + (1.0/12.0)*dat.rhoa*ux/dom.LBMDOM.Cs;
        f[ 9] = f[10] + (1.0/12.0)*dat.rhoa*ux/dom.LBMDOM.Cs;
        f[11] = f[12] + (1.0/12.0)*dat.rhoa*ux/dom.LBMDOM.Cs;
        f[13] = f[14] + (1.0/12.0)*dat.rhoa*ux/dom.LBMDOM.Cs;
        dom.LBMDOM.Vel[0][0][i][j] = OrthoSys::O;
        dom.LBMDOM.Rho[0][0][i][j] = 0.0;
        for (size_t k=0;k<dom.LBMDOM.Nneigh;k++)
        {
            dom.LBMDOM.Rho[0][0][i][j] +=  dom.LBMDOM.F[0][0][i][j][k];
            dom.LBMDOM.Vel[0][0][i][j] +=  dom.LBMDOM.F[0][0][i][j][k]*dom.LBMDOM.C[k];
        }
        dom.LBMDOM.Vel[0][0][i][j] *= dom.LBMDOM.Cs/dom.LBMDOM.Rho[0][0][i][j];
    }

    // x+ boundary
    #pragma omp parallel for schedule(static) num_threads(dom.LBMDOM.Nproc)
    for (size_t i=0; i<ny; ++i) //y-axis
    for (size_t j=0; j<nz; ++j) //z-axis
    {
        double * f = dom.LBMDOM.F[0][nx-1][i][j];
        double ux  =  dom.LBMDOM.Cs*(-dat.rhoa+f[0]+f[3]+f[4]+f[5]+f[6]+2.0*(f[1]+f[7]+f[9]+f[11]+f[13]))/dat.rhoa;
        f[ 2] = f[ 1] -  (2.0/3.0)*dat.rhoa*ux/dom.LBMDOM.Cs;
        f[ 8] = f[ 7] - (1.0/12.0)*dat.rhoa*ux/dom.LBMDOM.Cs;
        f[10] = f[ 9] - (1.0/12.0)*dat.rhoa*ux/dom.LBMDOM.Cs;
        f[12] = f[11] - (1.0/12.0)*dat.rhoa*ux/dom.LBMDOM.Cs;
        f[14] = f[13] - (1.0/12.0)*dat.rhoa*ux/dom.LBMDOM.Cs;
        dom.LBMDOM.Vel[0][nx-1][i][j] = OrthoSys::O;
        dom.LBMDOM.Rho[0][nx-1][i][j] = 0.0;
        for (size_t k=0;k<dom.LBMDOM.Nneigh;k++)
        {
            dom.LBMDOM.Rho[0][nx-1][i][j] +=  dom.LBMDOM.F[0][nx-1][i][j][k];
            dom.LBMDOM.Vel[0][nx-1][i][j] +=  dom.LBMDOM.F[0][nx-1][i][j][k]*dom.LBMDOM.C[k];
        }
        dom.LBMDOM.Vel[0][nx-1][i][j] *= dom.LBMDOM.Cs/dom.LBMDOM.Rho[0][nx-1][i][j];
    }

    // y- boundary
    #pragma omp parallel for schedule(static) num_threads(dom.LBMDOM.Nproc)
    for (size_t i=0; i<nx; ++i) //x-axis
    for (size_t j=0; j<nz; ++j) //z-axis
    {
        double * f = dom.LBMDOM.F[0][i][0][j];
        double uy  = -dom.LBMDOM.Cs*(-dat.rhoa+f[0]+f[1]+f[2]+f[5]+f[6]+2.0*(f[4]+f[8]+f[10]+f[11]+f[13]))/dat.rhoa;
        f[ 3] = f[ 4] +  (2.0/3.0)*dat.rhoa*uy/dom.LBMDOM.Cs;
        f[ 7] = f[ 8] + (1.0/12.0)*dat.rhoa*uy/dom.LBMDOM.Cs;
        f[ 9] = f[10] + (1.0/12.0)*dat.rhoa*uy/dom.LBMDOM.Cs;
        f[12] = f[11] + (1.0/12.0)*dat.rhoa*uy/dom.LBMDOM.Cs;
        f[14] = f[13] + (1.0/12.0)*dat.rhoa*uy/dom.LBMDOM.Cs;
        dom.LBMDOM.Vel[0][i][0][j] = OrthoSys::O;
        dom.LBMDOM.Rho[0][i][0][j] = 0.0;
        for (size_t k=0;k<dom.LBMDOM.Nneigh;k++)
        {
            dom.LBMDOM.Rho[0][i][0][j] +=  dom.LBMDOM.F[0][i][0][j][k];
            dom.LBMDOM.Vel[0][i][0][j] +=  dom.LBMDOM.F[0][i][0][j][k]*dom.LBMDOM.C[k];
        }
        dom.LBMDOM.Vel[0][i][0][j] *= dom.LBMDOM.Cs/dom.LBMDOM.Rho[0][i][0][j];
    }

    // y+ boundary
    #pragma omp parallel for schedule(static) num_threads(dom.LBMDOM.Nproc)
    for (size_t i=0; i<nx; ++i) //x-axis
    for (size_t j=0; j<nz; ++j) //z-axis
    {
        double * f = dom.LBMDOM.F[0][i][ny-1][j];
        double uy  =  dom.LBMDOM.Cs*(-dat.rhoa+f[0]+f[1]+f[2]+f[5]+f[6]+2.0*(f[3]+f[7]+f[9]+f[12]+f[14]))/dat.rhoa;
        f[ 4] = f[ 3] -  (2.0/3.0)*dat.rhoa*uy/dom.LBMDOM.Cs;
        f[ 8] = f[ 7] - (1.0/12.0)*dat.rhoa*uy/dom.LBMDOM.Cs;
        f[10] = f[ 9] - (1.0/12.0)*dat.rhoa*uy/dom.LBMDOM.Cs;
        f[11] = f[12] - (1.0/12.0)*dat.rhoa*uy/dom.LBMDOM.Cs;
        f[13] = f[14] - (1.0/12.0)*dat.rhoa*uy/dom.LBMDOM.Cs;
        dom.LBMDOM.Vel[0][i][ny-1][j] = OrthoSys::O;
        dom.LBMDOM.Rho[0][i][ny-1][j] = 0.0;
        for (size_t k=0;k<dom.LBMDOM.Nneigh;k++)
        {
            dom.LBMDOM.Rho[0][i][ny-1][j] +=  dom.LBMDOM.F[0][i][ny-1][j][k];
            dom.LBMDOM.Vel[0][i][ny-1][j] +=  dom.LBMDOM.F[0][i][ny-1][j][k]*dom.LBMDOM.C[k];
        }
        dom.LBMDOM.Vel[0][i][ny-1][j] *= dom.LBMDOM.Cs/dom.LBMDOM.Rho[0][i][ny-1][j];
    }

    //z+ boundary
    #pragma omp parallel for schedule(static) num_threads(dom.LBMDOM.Nproc)
    for (size_t i=0; i<nx; ++i) //x-axis
    for (size_t j=0; j<ny; ++j) //y-axis
    {
        double r = dom.LBMDOM.dx*sqrt((i-0.5*nx)*(i-0.5*nx)+(j-0.5*ny)*(j-0.5*ny)); //distance fromt he center of the face
        double rho,uz;
        double * f = dom.LBMDOM.F[0][i][j][nz-1];
        if (r<=dat.Ro)
        {
            //uz  = -dat.umax*(dat.Ro*dat.Ro-r*r)/(dat.Ro*dat.Ro)*std::min(10.0*dom.Time/dat.Tf,1.0);
            uz  = -dat.umax*(dat.Ro*dat.Ro-r*r)/(dat.Ro*dat.Ro)*sin(10.0*M_PI*dom.Time/dat.Tf);
            if (dom.Time>dat.Tf/10.0) uz = 0.0;
            rho = (f[0]+f[1]+f[2]+f[3]+f[4]+2.0*(f[5]+f[7]+f[10]+f[11]+f[14]))/(1+uz/dom.LBMDOM.Cs);
        }
        else 
        {
            rho = dat.rhoa;
            uz  =  dom.LBMDOM.Cs*(-dat.rhoa+f[0]+f[1]+f[2]+f[3]+f[4]+2.0*(f[5]+f[7]+f[10]+f[11]+f[14]))/dat.rhoa;
        }

        f[ 6] = f[ 5] -  (2.0/3.0)*rho*uz/dom.LBMDOM.Cs;
        f[ 8] = f[ 7] - (1.0/12.0)*rho*uz/dom.LBMDOM.Cs;
        f[ 9] = f[10] - (1.0/12.0)*rho*uz/dom.LBMDOM.Cs;
        f[12] = f[11] - (1.0/12.0)*rho*uz/dom.LBMDOM.Cs;
        f[13] = f[14] - (1.0/12.0)*rho*uz/dom.LBMDOM.Cs;

        dom.LBMDOM.Vel[0][i][j][nz-1] = OrthoSys::O;
        dom.LBMDOM.Rho[0][i][j][nz-1] = 0.0;
        for (size_t k=0;k<dom.LBMDOM.Nneigh;k++)
        {
            dom.LBMDOM.Rho[0][i][j][nz-1] +=  dom.LBMDOM.F[0][i][j][nz-1][k];
            dom.LBMDOM.Vel[0][i][j][nz-1] +=  dom.LBMDOM.F[0][i][j][nz-1][k]*dom.LBMDOM.C[k];
        }
        dom.LBMDOM.Vel[0][i][j][nz-1] *= dom.LBMDOM.Cs/dom.LBMDOM.Rho[0][i][j][nz-1];

    }
}

int main(int argc, char **argv) try
{ 
    LBMMPM::Domain dom;
    UserData dat;
    dom.UserData = &dat;
   
    size_t Nproc = 0.75*omp_get_max_threads();
    //size_t Nproc = 1;
    if (argc>1) Nproc = atoi(argv[1]);

    size_t ndiv = 10; //number of divisions per x lenght
    double Nu   = 0.4; //Poisson ratio
    double  E   = 6.8e5; //Young Modulus
    double  K   = E/(3.0*(1.0-2.0*Nu)); //Bulk Modulus
    double  G   = E/(2.0*(1.0+Nu)); // Shear Modulus
    double rho  = 1030.0; // Density of solid phase  
    double Cs   = sqrt(K/rho); // Speed of sound of cornea
    double dx   = 1.0e-4; //grid size in (m) for the LBM mesh
    double rhoa = 1.2041; //density of air
    double eta  = 1.48e-5; //kinematic viscosity of air
    size_t nx   = 100;     // lattice size
    size_t ny   = 100;     // lattice size
    size_t nz   = 120;     // lattice size in z equal to 1.2 cm
    double Ro   = 1.2e-3;  // Radius of the inlet in m
    double umax = 167.8;   // maximun inlet velocity
    double cair = 343.0;   // air velocity of sound
    double Tf   = 10.0*3.21e-5; // Final Time

    dom.MPMDOM.AddFromOBJMesh(-1,"eye.obj",OrthoSys::O,OrthoSys::O,OrthoSys::e0,M_PI,rho,1.0e-3,ndiv); 
    Vec3_t pos(5.0e-3,5.0e-3,0.5e-3);
    dom.MPMDOM.PosMesh(pos);
    Vec3_t Xmin,Xmax;
    dom.MPMDOM.BoundingBox(Xmin,Xmax);
    double h = dom.MPMDOM.ResizeDomainMesh(Xmin-Vec3_t(0.5e-2,0.5e-2,0.5e-2),Xmax+Vec3_t(0.5e-2,0.5e-2,0.5e-2),1.0);

    double dt   = std::min(1.0*h/Cs,dx/(sqrt(3)*cair)); //h is the minimum mpm element size and dt is the time step
    double Bc   = h;      //Mimimun level for the boundary condition

    dom.LBMmesh(D3Q15,eta,iVec3_t(nx,ny,nz),dx,dt);

    //Setting properties for the material points
    for (size_t ip=0; ip < dom.MPMDOM.Particles.Size(); ip++)
    {
        dom.MPMDOM.Particles[ip]->K = K;
        dom.MPMDOM.Particles[ip]->G = G;

        if (dom.MPMDOM.Particles[ip]->x(2)<Bc)
        {
            dom.MPMDOM.Particles[ip]->FixVeloc();
            dom.MPMDOM.Particles[ip]->Tag = -2;
        }
    }

    //Fixing positions and forces over corner points 
    for (size_t ic=0; ic < dom.MPMDOM.Corners.Size(); ic++)
    {
        if ((*dom.MPMDOM.Corners[ic]->x)(2)<Bc)
        {
            dom.MPMDOM.Corners[ic]->FixVeloc();
        }
    }

    //Setting properties for nodes
    for (size_t in=0; in < dom.MPMDOM.Nnodes; in++)
    {
        Vec3_t xn;
        dom.MPMDOM.NodePosition(in,xn);
        if (xn(2)<Bc)
        {
            dom.MPMDOM.Nodes[in].FixVeloc();
        }
    }

    //Setting initial conditions for the fluid
    for (size_t ix=0;ix<nx;ix++)
    for (size_t iy=0;iy<ny;iy++)
    for (size_t iz=0;iz<nz;iz++)
    {
        Vec3_t v(0.0,0.0,0.0);
        iVec3_t idx(ix,iy,iz);
        dom.LBMDOM.Initialize(0,idx,rhoa/*rho*/,v);
        if (iz==0)
        {
            dom.LBMDOM.IsSolid[0][ix][iy][iz] = true;
        }
    }  

    dat.Tf    = Tf;
    dat.Ro    = Ro;
    dat.rhoa  = rhoa;
    dat.umax  = umax;
    dom.Solve(Tf,Tf/100,Setup,NULL,"cornea",true,Nproc);

}
MECHSYS_CATCH
//
