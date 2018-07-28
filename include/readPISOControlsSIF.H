    dictionary piso = regions3d[0].solutionDict().subDict("PISO");

    int nOuterCorr =
        piso.lookupOrDefault<int>("nOuterCorrectors", 1);
