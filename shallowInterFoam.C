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
    shallowInterFoam

Description

    Coupled 2D shallow water and 3D Navier-Stokes solver with free surface, based on the 2D solver shallowFoam and the 3D solver interFoam. The explicit bi-directional coupling is implemented via the boundary conditions. For a description of the solver and as reference please see:

    Mintgen, F. & Manhart, M.: A bi-directional coupling of 2D shallow water and 3D Reynolds-averaged Navier–Stokes models. March 2018. Journal of Hydraulic Research. DOI: 10.1080/00221686.2017.1419989

    A more detailed description can be found in:

    Mintgen, F.: Coupling of Shallow and Non-Shallow Flow Solvers - An Open Source Framework. 2018. PhD dissertation. Technical University of Munich.

    Copies of paper and thesis are available upon request to f.mintgen@tum.de

Disclaimer

    This is research code, if I would have the time, I would rewrite a number of things. Parts of the boundary conditions, especially on the 3D side, are redundant or obsolete and thus can get a bit messy. A cleaned-up version of the 2D solver shallowFoam can be found on gitHub.

Author
    Florian Mintgen, 2014 - 2018

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "turbulenceModel.H"
#include "regionProperties.H"

// Header files aus interFoam
#include "MULES.H"
#include "subCycle.H"
#include "interfaceProperties.H"
#include "twoPhaseMixture.H"
#include "turbulenceModel.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#include "setRootCase.H"
#include "createTime.H"

regionProperties rp(runTime);

#include "create2dMeshes.H"
#include "create3dMeshes.H"

#include "create2dFields.H"
#include "create3dFields.H"

#include "initContinuityErrs.H"

#include "readTimeControls.H"
#include "CourantNoMultiRegion.H"
#include "setInitialDeltaT.H"

Info<< "\nStarting loop over all 3D-regions\n" << endl;

 forAll(regions3d, i)
   {
     #include "set3dFields.H"
     #include "readRegion3dPISOControls.H"
     #include "correctPhi.H"
   }
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {

        #include "readTimeControls.H"
        #include "readPISOControlsSIF.H"
        #include "CourantNoMultiRegion.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

	// Pressure-velocity corrector
	for (int oCorr=0; oCorr<nOuterCorr; oCorr++)
	{

	  forAll(regions2d, i)
	    {
	      Info<< "\nSolving for 2D region " << regions2d[i].name() << endl;
              #include "set2dFields.H"
	      #include "solve2d.H"
	      Info<< "\nSolved for 2D region " << regions2d[i].name() << endl;
	    }

	  forAll(regions3d, i)
	    {
	      Info<< "\nSolving for 3D region " << regions3d[i].name() << endl;
              #include "set3dFields.H"
              #include "readRegion3dPISOControls.H"
              #include "solve3d.H"
	      Info<< "\nSolved for 3D region " << regions3d[i].name() << endl;
            }
	}

///////////////////////////////////////////////////////////////////////////////

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
