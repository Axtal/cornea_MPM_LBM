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
#include <mechsys/mpm/Domain.h>
#include <math.h>
#include <iostream>
#include <fstream>

struct UserData
{
    double Tf;
    double Fb;
    Array<MPM::Corner   *> ForcePar; //Forced particles
};

void Setup (MPM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    
    double F = std::min(10.0*dat.Fb*dom.Time/dat.Tf,dat.Fb)/dat.ForcePar.Size();

    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
    for (size_t ip=0; ip<dat.ForcePar.Size(); ip++)
    {
        dat.ForcePar[ip]->h = Vec3_t(0.0,0.0,-F);
    }
}

int main(int argc, char **argv) try
{
    MPM::Domain dom;
    UserData dat;
    dom.UserData = &dat;
   
    size_t Nproc = 0.75*omp_get_max_threads();
    //size_t Nproc = 1;
    if (argc>1) Nproc = atoi(argv[1]);

    size_t ndiv = 10; //number of divisions per x lenght
    double Nu   = 0.4; //Poisson ratio
    double  E   = 2.0e5; //Young Modulus
    double  K   = E/(3.0*(1.0-2.0*Nu)); //Bulk Modulus
    double  G   = E/(2.0*(1.0+Nu)); // Shear Modulus
    double rho  = 1030.0; // Density of solid phase  
    double Cs   = sqrt(E/rho); // Speed of sound


    dom.AddFromOBJMesh(-1,"eye.obj",OrthoSys::O,OrthoSys::O,OrthoSys::e0,M_PI,rho,1.0e-3,ndiv); 
    Vec3_t pos(5.0e-3,5.0e-3,0.5e-3);
    dom.PosMesh(pos);
    Vec3_t Xmin,Xmax;
    dom.BoundingBox(Xmin,Xmax);
    double h = dom.ResizeDomainMesh(Xmin-Vec3_t(0.5e-2,0.5e-2,0.5e-2),Xmax+Vec3_t(0.5e-2,0.5e-2,0.5e-2),1.0);

    double dt   = 2.0*h/Cs; //h is the minimum mpm element size and dt is the time step
    double Bc   = h;      //Mimimun level for the boundary condition
    double Bcup = Xmax(2)-h;



    //Setting properties for the material points
    for (size_t ip=0; ip < dom.Particles.Size(); ip++)
    {
        dom.Particles[ip]->K = K;
        dom.Particles[ip]->G = G;

        if (dom.Particles[ip]->x(2)<Bc)
        {
            dom.Particles[ip]->FixVeloc();
            dom.Particles[ip]->Tag = -2;
        }
    }

    //Fixing positions and forces over corner points 
    for (size_t ic=0; ic < dom.Corners.Size(); ic++)
    {
        if ((*dom.Corners[ic]->x)(2)<Bc)
        {
            dom.Corners[ic]->FixVeloc();
        }
        if ((*dom.Corners[ic]->x)(2)>Bcup)
        {
            dat.ForcePar.Push(dom.Corners[ic]);
        }
    }

    

    //Setting properties for nodes
    for (size_t in=0; in < dom.Nnodes; in++)
    {
        Vec3_t xn;
        dom.NodePosition(in,xn);
        if (xn(2)<Bc)
        {
            dom.Nodes[in].FixVeloc();
        }
    }

    double Tf = 1.0; // Final Time
    dat.Tf    = Tf;
    dat.Fb = 1.0e-1;  // Force applied at the top
    dom.Solve(Tf,dt,Tf/100,Setup,NULL,"cornea_mpm",true,Nproc);

}
MECHSYS_CATCH
//
