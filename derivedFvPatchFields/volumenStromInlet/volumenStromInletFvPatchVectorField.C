/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2006-2008 OpenCFD Ltd.
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

Author
    KM-Turbulenz GmbH, 2009
    dervied from flowRateInletVelocityFvPatchVectorField.C
\*---------------------------------------------------------------------------*/

#include "volumenStromInletFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::
volumenStromInletFvPatchVectorField::
volumenStromInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    volumenStrom_(0),
    HName_("H")
{}


Foam::
volumenStromInletFvPatchVectorField::
volumenStromInletFvPatchVectorField
(
    const volumenStromInletFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    volumenStrom_(ptf.volumenStrom_),
    HName_(ptf.HName_)
{}


Foam::
volumenStromInletFvPatchVectorField::
volumenStromInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict),
    volumenStrom_(readScalar(dict.lookup("volumenStrom"))),
    HName_("H")
{

    if (dict.found("H"))
    {
        dict.lookup("H") >> HName_;
    }
}


Foam::
volumenStromInletFvPatchVectorField::
volumenStromInletFvPatchVectorField
(
    const volumenStromInletFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    volumenStrom_(ptf.volumenStrom_),
    HName_(ptf.HName_)
{}


Foam::
volumenStromInletFvPatchVectorField::
volumenStromInletFvPatchVectorField
(
    const volumenStromInletFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    volumenStrom_(ptf.volumenStrom_),
    HName_(ptf.HName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::volumenStromInletFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }


    vectorField n = patch().nf();

    const volScalarField& H = db().lookupObject<volScalarField>(HName_);
    
    const scalarField pointsZ = patch().patch().points().component(vector::Z);
    const scalar heightPatch = max(pointsZ) - min(pointsZ);
    
    scalar surfaceSumH = gSum(H.boundaryField()[patch().index()]*patch().magSf());
    scalar Umean = volumenStrom_/surfaceSumH * heightPatch;
    scalarField phip = patch().lookupPatchField<surfaceScalarField, scalar>("phi");
    scalar surfaceSumPhip = gSum(phip*H.boundaryField()[patch().index()]) / heightPatch;
    Info << "Inflow-Volumenstrom - Q: " <<  surfaceSumPhip << endl;   

    operator==(-n*Umean*H.boundaryField()[patch().index()]);

    fixedValueFvPatchField<vector>::updateCoeffs();
}


void Foam::volumenStromInletFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);

    os.writeKeyword("volumenStrom") << volumenStrom_
        << token::END_STATEMENT << nl;

    if (HName_ != "H")
    {
        os.writeKeyword("H") << HName_ << token::END_STATEMENT << nl;
    }

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       volumenStromInletFvPatchVectorField
   );
}


// ************************************************************************* //
