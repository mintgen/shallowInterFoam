    // Initialise fluid field pointer lists
    PtrList<volScalarField> pd3d(regions3d.size());
    PtrList<volScalarField> alpha13d(regions3d.size());
    PtrList<volVectorField> U3d(regions3d.size());
    PtrList<uniformDimensionedVectorField> g3d(regions3d.size());
    PtrList<surfaceScalarField> phi3d(regions3d.size());
    PtrList<twoPhaseMixture> twoPhaseProperties3d(regions3d.size());
    PtrList<dimensionedScalar> rho13d(regions3d.size());
    PtrList<dimensionedScalar> rho23d(regions3d.size());
    PtrList<volScalarField> rho3d(regions3d.size());
    PtrList<surfaceScalarField> rhoPhi3d(regions3d.size());
    PtrList<volScalarField> gh3d(regions3d.size());
    PtrList<surfaceScalarField> ghf3d(regions3d.size());
    PtrList<volScalarField> p3d(regions3d.size());
    List<label> pdRefCell3d(regions3d.size(),0);
    List<scalar> pdRefValue3d(regions3d.size(),0.0);
    List<label> pRefCell3d(regions3d.size(),0);
    List<scalar> pRefValue3d(regions3d.size(),0.0);
    PtrList<interfaceProperties> interface3d(regions3d.size());
    PtrList<incompressible::turbulenceModel> turbulence3d(regions3d.size());
    PtrList<volSymmTensorField> Reff3d(regions3d.size());

    // Populate 3D field pointer lists
    forAll(regions3d, i)
    {

       Info<< nl << "Reading field pd3d\n";
       pd3d.set
	 (
	  i,
	  new volScalarField
	  (
	   IOobject
	   (
	    "pd",
	    runTime.timeName(),
	    regions3d[i],
	    IOobject::MUST_READ,
	    IOobject::AUTO_WRITE
	    ),
	   regions3d[i]
	   )
	  );

       Info<< nl << "Reading field alpha13d\n";
       alpha13d.set
	 (
	  i,
	  new volScalarField
	  (
	   IOobject
	   (
	    "alpha1",
	    runTime.timeName(),
	    regions3d[i],
	    IOobject::MUST_READ,
	    IOobject::AUTO_WRITE
	    ),
	   regions3d[i]
	   )
	  );

       Info<< nl << "Reading field U3d\n";
       U3d.set
	 (
	  i,
	  new volVectorField
	  (
	   IOobject
	   (
	    "U",
	    runTime.timeName(),
	    regions3d[i],
	    IOobject::MUST_READ,
	    IOobject::AUTO_WRITE
	    ),
	   regions3d[i]
	   )
	  );
       
       Info<< nl << "Reading field g3d\n";
       g3d.set
	 (
	  i,
	  new uniformDimensionedVectorField
	  (
	   IOobject
	   (
	    "g",
	    runTime.constant(),
	    regions3d[i],
	    IOobject::MUST_READ,
	    IOobject::NO_WRITE
	    )
	   )
	  );
       
       Info<< nl << "Reading field phi3d\n";
       phi3d.set
	 (
	  i,
	  new surfaceScalarField
	  (
	   IOobject
	   (
	    "phi",
	    runTime.timeName(),
	    regions3d[i],
	    IOobject::READ_IF_PRESENT,
	    IOobject::NO_WRITE
	    ),
	   linearInterpolate(U3d[i]) & regions3d[i].Sf()
	   )
	  );
       
       Info<< nl << "Adding to twoPhaseProperties3d\n";
       twoPhaseProperties3d.set
	 (
	  i,
	  new twoPhaseMixture
	  (
	   U3d[i],
	   phi3d[i],
	   "alpha1"
	   )
	  );

	rho13d.set
	(
	    i,
	    new dimensionedScalar
	    (
	        twoPhaseProperties3d[i].rho1()
	    )
	);
       
	rho23d.set
	(
	    i,
	    new dimensionedScalar
	    (
	        twoPhaseProperties3d[i].rho2()
	    )
	);
       
       Info<< nl << "Create field rho3d\n";
       rho3d.set
	 (
	  i,
	  new volScalarField
	  (
	   IOobject
	   (
	    "rho",
	    runTime.timeName(),
	    regions3d[i],
	    IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE       
	    ),
	   alpha13d[i]*rho13d[i] + (scalar(1) - alpha13d[i])*rho23d[i],
	   alpha13d[i].boundaryField().types()
	   )
	  );
       rho3d[i].oldTime();

    // Mass flux
    // Initialisation does not matter because rhoPhi is reset after the
    // alpha1 solution before it is used in the U equation.
       Info<< nl << "Reading field phi3d\n";
       rhoPhi3d.set
	 (
	  i,
	  new surfaceScalarField
	  (
	   IOobject
	   (
	    "rhoPhi",
	    runTime.timeName(),
	    regions3d[i],
	    IOobject::NO_READ,
	    IOobject::NO_WRITE
	    ),
	   rho13d[i]*phi3d[i]
	   )
	  );
       
       Info<< nl << "Calculating field g.h\n";
       gh3d.set
	 (
	  i,
	  new volScalarField
	  (
	   "gh",
	   g3d[i] & regions3d[i].C()
	   )
	  );

       ghf3d.set
	 (
	  i,
	  new surfaceScalarField
	  (
	   "gh",
	   g3d[i] & regions3d[i].Cf()
	   )
	  );

       Info<< nl << "Create field p3d\n";
       p3d.set
	 (
	  i,
	  new volScalarField
	  (
	   IOobject
	   (
	    "p",
	    runTime.timeName(),
	    regions3d[i],
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
	    ),
	   pd3d[i] + rho3d[i]*gh3d[i]
	   )
	  );


    setRefCell
      (
       pd3d[i],
       regions3d[i].solutionDict().subDict("PISO"),
       pdRefCell3d[i],
       pdRefValue3d[i]
       );

    setRefCell
      (
       p3d[i],
       regions3d[i].solutionDict().subDict("PISO"),
       pRefCell3d[i],
       pRefValue3d[i]
       );

    /*
    //This doesn't work properly, but since we do not have a closed system, no reference pressure is needed.
    Info<< nl << "*** Reference data for pressure for region "
	<<  regions3d[i].name() << nl
	<< "pdRefCell: " << pdRefCell3d[i] << nl 
	<< "pdRefValue: " << pdRefValue3d[i] << nl
	<< "pRefCell: " << pRefCell3d[i] << nl 
	<< "pRefValue: " << pRefValue3d[i] << nl;
    */

    Info<< nl << "Adding to interfaceProperties\n";
    interface3d.set
      (
       i,
       new interfaceProperties
       (
	alpha13d[i],
	U3d[i],
	twoPhaseProperties3d[i]
	)
       );
    
    Info<< "    Adding to turbulence\n" << endl;
    turbulence3d.set
      (
       i,
       autoPtr<incompressible::turbulenceModel>
       (
	incompressible::turbulenceModel::New
	(
	 U3d[i],
	 phi3d[i],
	 twoPhaseProperties3d[i]
	 )
	).ptr()
       );
    
       Info<< nl << "Create field Reff3d\n";
       Reff3d.set
	 (
	  i,
	  new volSymmTensorField
	  (
	   IOobject
	   (
	    "Reff",
	    runTime.timeName(),
	    regions3d[i],
            IOobject::NO_READ,
            IOobject::NO_WRITE
	    ),
	    turbulence3d[i].devReff()
	   )
	  );

    }
