    PtrList<fvMesh> regions2d(rp.region2dNames().size());

    forAll(rp.region2dNames(), i)
    {
        Info<< "Create 2D mesh for region " << rp.region2dNames()[i]
            << " for time = " << runTime.timeName() << nl << endl;

        regions2d.set
        (
            i,
            new fvMesh
            (
                IOobject
                (
                    rp.region2dNames()[i],
                    runTime.timeName(),
                    runTime,
                    IOobject::MUST_READ
                )
            )
        );
    }
