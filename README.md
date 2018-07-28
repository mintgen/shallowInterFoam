# shallowInterFoam

OpenFOAM version: foam-extend-3.1

Description

Coupled 2D shallow water and 3D Navier-Stokes solver with free surface, based on the 2D solver shallowFoam and the 3D solver interFoam. The explicit bi-directional coupling is implemented via the boundary conditions. For a description of the solver and as reference please see:

Mintgen, F. & Manhart, M.: A bi-directional coupling of 2D shallow water and 3D Reynolds-averaged Navier–Stokes models. March 2018. Journal of Hydraulic Research. DOI: 10.1080/00221686.2017.1419989

A more detailed description can be found in:

Mintgen, F.: Coupling of Shallow and Non-Shallow Flow Solvers - An Open Source Framework. 2018. PhD dissertation. Technical University of Munich.

Copies of paper and thesis are available upon request to f.mintgen@tum.de

Disclaimer

This is research code, if I would have the time, I would rewrite a number of things. Parts of the boundary conditions, especially on the 3D side, are redundant or obsolete and thus can get a bit messy. A cleaned-up version of the 2D solver shallowFoam can be found on gitHub.
