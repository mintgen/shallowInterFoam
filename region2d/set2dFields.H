    const fvMesh& mesh = regions2d[i];

    const dimensionedScalar& dim_s = dim_s2d[i];
    const dimensionedScalar& dim_m = dim_m2d[i];
    const dimensionedScalar& g = g2d[i];
    const dimensionedScalar& Hdry = Hdry2d[i];
    const dimensionedScalar& Hdry2 = Hdry22d[i];
    const dimensionedScalar& Cnu = Cnu2d[i];
    const dimensionedScalar& nutmax = nutmax2d[i];

    volScalarField& H = H2d[i];
    volScalarField& Hclip = Hclip2d[i];
    const volScalarField& S = S2d[i];
    volScalarField& HpS = HpS2d[i];
    const volScalarField& kst = kst2d[i];
    volVectorField& HU = HU2d[i];
    volVectorField& U = U2d[i];
    surfaceScalarField& phi = phi2d[i];
    volScalarField& alphaWSS = alphaWSS2d[i];
    volScalarField& nut = nut2d[i];
    volScalarField& nassTrocken = nassTrocken2d[i];
    volVectorField& froude = froude2d[i];
