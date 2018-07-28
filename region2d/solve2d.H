/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    shallowFoam

Description
    Transient solver for the shallow water equations

Author
    KM-Turbulenz GmbH, 2009
    Florian Mintgen, 2012

\*---------------------------------------------------------------------------*/

{

        // Transport equation for H
        fvScalarMatrix HEqn
        (
	    fvm::ddt(H)
	    + fvm::div(phi,H)
        );

        HEqn.solve(mesh.solutionDict().solver(H.name()));

        H.correctBoundaryConditions();

        // Hclip = max (H, Hdry)
        Hclip = (H-Hdry)*pos(H-Hdry) + Hdry;

	// Bottom friction from Manning-Strickler (kst = 1/n)
        alphaWSS = (g*dim_s*dim_s/dim_m)*mag(U)/pow(kst,2.0)/pow(Hclip/dim_m,1.0/3.0)/Hclip;

	// Depth-averaged parabolic eddy viscosity
        // nu_t = C_nu * u_tau * H
        nut = Cnu*sqrt(alphaWSS*mag(HU))*H;
        nut = nut*pos(nutmax-nut) + nutmax*pos(nut-nutmax);
        nut.correctBoundaryConditions();

	// Zero-Gradient for faces on wet-dry interface
	#include "zeroGradientNassTrocken.H"
	
	// Transport equation for HU
         fvVectorMatrix HUEqn
         (
             fvm::ddt(HU)                 // d(HU_i)/dt                    Temporal derivative
	     + fvm::div(phi, HU)          // + d/dx_j ( U_j * HU_i )       Convection
	     //	     + H*fvc::grad(gS)    // + H * d/dx_i ( g * S )        Source term due to bottom slope
	     //	     + H*fvc::grad(gH)    // + H * d/dx_i ( g * H )        Source term due to surface slope
             + H*gradgSpHf                // + H * d/dx_i ( g * (S + H) )  Combination of the two source terms
	     + fvm::Sp(alphaWSS,HU)       // + alphaWSS * HU               Bottom friction
	     - fvm::laplacian(nut,HU)     // - d^2/dx_j^2 (nut * HU_i)     Turbulent stresses
         );  

	HUEqn.solve();

        HU.correctBoundaryConditions();
        
///////////////////////////////////////////////////////////////////////////////

        // Calculation of U und phi for the next time step

        // if Hclip < Hdry2, then U = 0
        HU = HU*pos(Hclip - Hdry2); 
        U = HU/Hclip;
        
        phi = (fvc::interpolate(U) & mesh.Sf());

        // Froude-Number
        froude = HU / ( pow(Hclip,1.5) * pow(g,0.5) );

	// Calculation of bottom friction (with rho = 1000)
	//	tau_w = 1000 * alphaWSS * HU;

///////////////////////////////////////////////////////////////////////////////

}


// ************************************************************************* //
