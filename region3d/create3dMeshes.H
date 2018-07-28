    PtrList<fvMesh> regions3d(rp.region3dNames().size());

    forAll(rp.region3dNames(), i)
    {
        Info<< "Create 3D mesh for region " << rp.region3dNames()[i]
            << " for time = " << runTime.timeName() << nl << endl;

        regions3d.set
        (
            i,
            new fvMesh
            (
                IOobject
                (
                    rp.region3dNames()[i],
                    runTime.timeName(),
                    runTime,
                    IOobject::MUST_READ
                )
            )
        );
    }
