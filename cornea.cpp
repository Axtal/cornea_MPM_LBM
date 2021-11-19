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
};

void Setup (MPM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
}

int main(int argc, char **argv) try
{
    MPM::Domain dom;
    UserData dat;
    dom.UserData = &dat;
   
    size_t Nproc = 0.75*omp_get_max_threads();
    //size_t Nproc = 1;
    if (argc>1) Nproc = atoi(argv[1]);
    size_t ndiv = 10; //number of divisions per z length 
    //size_t ndiv = 5; //number of divisions per z length 
    double K    = 10.0e6;
    double Nu   = 0.3;
    double  E   = (3.0*(1.0-2.0*Nu))*K;
    double  G   = E/(2.0*(1.0+Nu));
    double rho  = 3000.0;
    double Cs   = sqrt(E/rho);
    double h    = 1.0/ndiv;
    double dt   = 0.1*h/Cs;
    double Lx   = 0.6;
    double Ly   = 0.6;

    dom.AddFromOBJMesh(-1,"test.obj",OrthoSys::O,OrthoSys::O,OrthoSys::e0,0.0,rho,1.0,ndiv); //Adding the fishtile

    double Bc = 0.1;
   /* 
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

    //dom.ParticleToNode();
    double Tf = 100.0; // Final Time

    dom.WriteXDMF("cornea");
  */

    //dom.Solve(Tf,dt,Tf/100,Setup,NULL,"cornea",true,Nproc);

}
MECHSYS_CATCH

