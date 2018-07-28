    // Initialise 2D field pointer lists
    PtrList<IOdictionary> transportProperties(regions2d.size());
    PtrList<dimensionedScalar> dim_s2d(regions2d.size());
    PtrList<dimensionedScalar> dim_m2d(regions2d.size());
    PtrList<dimensionedScalar> g2d(regions2d.size());
    PtrList<dimensionedScalar> Hdry2d(regions2d.size());
    PtrList<dimensionedScalar> Hdry22d(regions2d.size());
    PtrList<dimensionedScalar> Cnu2d(regions2d.size());
    PtrList<dimensionedScalar> nutmax2d(regions2d.size());
    PtrList<volScalarField> H2d(regions2d.size());
    PtrList<volScalarField> Hclip2d(regions2d.size());
    PtrList<volScalarField> S2d(regions2d.size());
    PtrList<volScalarField> HpS2d(regions2d.size());
    PtrList<volScalarField> kst2d(regions2d.size());
    PtrList<volVectorField> HU2d(regions2d.size());
    PtrList<volVectorField> U2d(regions2d.size());
    PtrList<surfaceScalarField> phi2d(regions2d.size());
    PtrList<volScalarField> alphaWSS2d(regions2d.size());
    PtrList<volScalarField> nut2d(regions2d.size());
    PtrList<volScalarField> nassTrocken2d(regions2d.size());
    PtrList<volVectorField> froude2d(regions2d.size());


    // Populate fluid field pointer lists
    forAll(regions2d, i)
    {
        Info<< "*** Reading 2D mesh transportPoperties for region "
            << regions2d[i].name() << nl << endl;

	transportProperties.set
	(
	    i,
	    new IOdictionary
	    (
	        IOobject
		(
		    "transportProperties",
		    runTime.constant(),
		    regions2d[i],
		    IOobject::MUST_READ,
		    IOobject::NO_WRITE
		)
	    )
	);
       
	dim_s2d.set
	(
	    i,
	    new dimensionedScalar
	    (
	        transportProperties[i].lookup("dim_s")
	    )
	);
       
	dim_m2d.set
	(
	    i,
	    new dimensionedScalar
	    (
	        transportProperties[i].lookup("dim_m")        
	    )
	);
       
	g2d.set
	(
	    i,
	    new dimensionedScalar
	    (
	        transportProperties[i].lookup("g")
	    )
	);
       
	Hdry2d.set
	(
	    i,
	    new dimensionedScalar
	    (
	        transportProperties[i].lookup("Hdry")
	    )
	);
       
	Hdry22d.set
	(
	    i,
	    new dimensionedScalar
	    (
	        transportProperties[i].lookup("Hdry2")
	    )
	);
       
	Cnu2d.set
	(
	    i,
	    new dimensionedScalar
	    (
	        transportProperties[i].lookup("Cnu")
	    )
	);
       
	nutmax2d.set
	(
	    i,
	    new dimensionedScalar
	    (
	        transportProperties[i].lookup("nutmax")
	    )
	);

	/*
       Info<< "dim_s for region " << regions2d[i].name() << ": " 
	   << dim_s2d[i] << endl;
       Info<< "dim_m for region " << regions2d[i].name() << ": " 
	   << dim_m2d[i] << endl;
       Info<< "g for region " << regions2d[i].name() << ": " 
	   << g2d[i] << endl;
       Info<< "Hdry for region " << regions2d[i].name() << ": " 
	   << Hdry2d[i] << endl;
       Info<< "Hdry2 for region " << regions2d[i].name() << ": " 
	   << Hdry22d[i] << endl;
       Info<< "Cnu for region " << regions2d[i].name() << ": " 
	   << Cnu2d[i] << endl;
       Info<< "nutmax for region " << regions2d[i].name() << ": " 
	   << nutmax2d[i] << endl;
	*/

       Info<< nl << "Adding to H2d\n";
       Info<< "Wasserstandshoehe\n" << endl;
       H2d.set
	 (
	  i,
	  new volScalarField
	  (
	   IOobject
	   (
	    "H",
	    runTime.timeName(),
	    regions2d[i],
	    IOobject::MUST_READ,
	    IOobject::AUTO_WRITE
	    ),
	   regions2d[i]
	   )
	  );

       Info<< nl << "Adding to Hclip2d\n" << endl;
       Hclip2d.set
	 (
	  i,
	  new volScalarField
	  (
	   IOobject
	   (
	    "Hclip",
	    runTime.timeName(),
	    regions2d[i],
	    IOobject::READ_IF_PRESENT,
	    IOobject::NO_WRITE
	    ),
	   (H2d[i]-Hdry2d[i])*pos(H2d[i]-Hdry2d[i]) + Hdry2d[i]        
	   )
	  );

       Info<< nl << "Adding to S2d\n";
       Info<< "Sohlverlauf\n" << endl;
       S2d.set
	 (
	  i,
	  new volScalarField
	  (
	   IOobject
	   (
	    "S",
	    runTime.timeName(),
	    regions2d[i],
	    IOobject::MUST_READ,
	    IOobject::AUTO_WRITE
	    ),
	   regions2d[i]
	   )
	  );

       Info<< nl << "Adding to HpS2d\n" << endl;
       HpS2d.set
	 (
	  i,
	  new volScalarField
	  (
	   IOobject
	   (
	    "HpS",
	    runTime.timeName(),
	    regions2d[i],
	    IOobject::READ_IF_PRESENT,
	    IOobject::NO_WRITE
	    ),
	   H2d[i]+S2d[i]        
	   )
	  );

       Info<< nl << "Adding to kst2d\n";
       Info<< "Stricklerbeiwert\n" << endl;
       kst2d.set
	 (
	  i,
	  new volScalarField
	  (
	   IOobject
	   (
	    "kst",
	    runTime.timeName(),
	    regions2d[i],
	    IOobject::MUST_READ,
	    IOobject::AUTO_WRITE
	    ),
	   regions2d[i]
	   )
	  );

       Info<< nl << "Adding to HU2d\n";
       Info<< "Spezifischer Durchfluss\n" << endl;
       HU2d.set
	 (
	  i,
	  new volVectorField
	  (
	   IOobject
	   (
	    "HU",
	    runTime.timeName(),
	    regions2d[i],
	    IOobject::MUST_READ,
	    IOobject::AUTO_WRITE
	    ),
	   regions2d[i]
	   )
	  );

       Info<< nl << "Adding to U2d\n";
       Info<< "Geschwindigkeit\n" << endl;
       U2d.set
	 (
	  i,
	  new volVectorField
	  (
	   IOobject
	   (
	    "U",
	    runTime.timeName(),
	    regions2d[i],
	    IOobject::READ_IF_PRESENT,
	    IOobject::AUTO_WRITE
	    ),
	   HU2d[i]/dim_m2d[i]
	   )
	  );

       Info<< nl << "Adding to phi2d\n" << endl;
       phi2d.set
	 (
	  i,
	  new surfaceScalarField
	  (
	   IOobject
	   (
	    "phi",
	    runTime.timeName(),
	    regions2d[i],
	    IOobject::READ_IF_PRESENT,
	    IOobject::AUTO_WRITE
	    ),
	   linearInterpolate(U2d[i]) & regions2d[i].Sf()
	   )
	  );

      Info<< nl << "Adding to alphaWSS2d\n";
       Info<< "Wandschubspannung durch U2d (oder so ähnlich...)\n" << endl;
       alphaWSS2d.set
	 (
	  i,
	  new volScalarField
	  (
	   IOobject
	   (
	    "alphaWSS",
	    runTime.timeName(),
	    regions2d[i],
	    IOobject::READ_IF_PRESENT,
	    IOobject::AUTO_WRITE
	    ),
	   (g2d[i]*dim_s2d[i]*dim_s2d[i]/dim_m2d[i])*mag(U2d[i])/pow(kst2d[i],2.0)/pow(Hclip2d[i]/dim_m2d[i],1.0/3.0)/Hclip2d[i]
	   )
	  );

      Info<< nl << "Adding to nut2d\n";
       Info<< "Turbulente Viskosität\n" << endl;
       nut2d.set
	 (
	  i,
	  new volScalarField
	  (
	   IOobject
	   (
	    "nut",
	    runTime.timeName(),
	    regions2d[i],
	    IOobject::READ_IF_PRESENT,
	    IOobject::AUTO_WRITE
	    ),
	   regions2d[i]               
	   )
	  );

       Info<< "nassTrocken\n" << endl;
       nassTrocken2d.set
	 (
	  i,
	  new volScalarField
	  (
	   IOobject
	   (
	    "nassTrocken",
	    runTime.timeName(),
	    regions2d[i],
	    IOobject::READ_IF_PRESENT,
	    IOobject::NO_WRITE
	    ),
	   pos(H2d[i]-Hdry2d[i])
	   )
	  );
       
       Info<< "Froude-Zahl\n" << endl;
       froude2d.set
	 (
	  i,
	  new volVectorField
	  (
	   IOobject
	   (
	    "froude",
	    runTime.timeName(),
	    regions2d[i],
	    IOobject::READ_IF_PRESENT,
	    IOobject::AUTO_WRITE
	    ),
               HU2d[i]/( pow(dim_m2d[i],1.5) * pow(g2d[i],0.5))
//               HU2d[i]/( pow(H2d[i],1.5) * pow(g2d[i],0.5))
               )
	  );

//    Info<< "froude for region " << regions2d[i].name() << ": " 
//	   << froude2d[i] << endl;
       
    Info<< "Felder Erzeugt\n" << endl;

    }
