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

#include "abflussWasserstandOutletFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::
abflussWasserstandOutletFvPatchVectorField::
abflussWasserstandOutletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    HName_("H")
{}


Foam::
abflussWasserstandOutletFvPatchVectorField::
abflussWasserstandOutletFvPatchVectorField
(
    const abflussWasserstandOutletFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    HName_(ptf.HName_)
{}


Foam::
abflussWasserstandOutletFvPatchVectorField::
abflussWasserstandOutletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict),
    HName_("H")
{

    if (dict.found("H"))
    {
        dict.lookup("H") >> HName_;
    }
}


Foam::
abflussWasserstandOutletFvPatchVectorField::
abflussWasserstandOutletFvPatchVectorField
(
    const abflussWasserstandOutletFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    HName_(ptf.HName_)
{}


Foam::
abflussWasserstandOutletFvPatchVectorField::
abflussWasserstandOutletFvPatchVectorField
(
    const abflussWasserstandOutletFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    HName_(ptf.HName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::abflussWasserstandOutletFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

//    Info << "Inf: Abfluss-Wasserstand" << endl;

    vectorField n = patch().nf();

    scalarField Hp = patch().lookupPatchField<volScalarField, scalar>(HName_);
    scalarField kstp = patch().lookupPatchField<volScalarField, scalar>("kst");
    const volScalarField& S = db().lookupObject<volScalarField>("S");

    scalarField Sint = patch().lookupPatchField<volScalarField, scalar>("S").patchInternalField();
//    scalarField Sp = patch().lookupPatchField<volScalarField, scalar>("S");

    const volVectorField gradS = fvc::grad(S);

//    Info << "Inf: Abfluss-Wasserstand" <<  gradS << endl;
    vectorField gradSpatch = patch().lookupPatchField<volVectorField, vector>("grad(S)").patchInternalField();
//    scalarField HUp = pos(-n&gradSpatch)*sqrt(mag(n & gradSpatch))*kstp*pow(Hp,(5.0/3.0))*patch().magSf();
    scalarField phip = patch().lookupPatchField<surfaceScalarField, scalar>("phi");
    scalar surfaceSumPhip = gSum(phip*Hp);

    Info << "Abfluss - Q: " <<  surfaceSumPhip << endl;
    operator==(pos(-n&gradSpatch)*sqrt(mag(n & gradSpatch))*n*kstp*pow(Hp,(5.0/3.0)));

    fixedValueFvPatchField<vector>::updateCoeffs();
}


void Foam::abflussWasserstandOutletFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);

//     os.writeKeyword("energieHoehe") << energieHoehe_
//         << token::END_STATEMENT << nl;

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
       abflussWasserstandOutletFvPatchVectorField
   );
}


// ************************************************************************* //
