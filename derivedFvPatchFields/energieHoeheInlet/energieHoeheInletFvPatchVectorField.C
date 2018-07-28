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

#include "energieHoeheInletFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::
energieHoeheInletFvPatchVectorField::
energieHoeheInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    energieHoehe_(0),
    HName_("H")
{}


Foam::
energieHoeheInletFvPatchVectorField::
energieHoeheInletFvPatchVectorField
(
    const energieHoeheInletFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    energieHoehe_(ptf.energieHoehe_),
    HName_(ptf.HName_)
{}


Foam::
energieHoeheInletFvPatchVectorField::
energieHoeheInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict),
    energieHoehe_(readScalar(dict.lookup("energieHoehe"))),
    HName_("H")
{

    if (dict.found("H"))
    {
        dict.lookup("H") >> HName_;
    }
}


Foam::
energieHoeheInletFvPatchVectorField::
energieHoeheInletFvPatchVectorField
(
    const energieHoeheInletFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    energieHoehe_(ptf.energieHoehe_),
    HName_(ptf.HName_)
{}


Foam::
energieHoeheInletFvPatchVectorField::
energieHoeheInletFvPatchVectorField
(
    const energieHoeheInletFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    energieHoehe_(ptf.energieHoehe_),
    HName_(ptf.HName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::energieHoeheInletFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

//    Info << "Inf: " << energieHoehe_ << endl;

    vectorField n = patch().nf();

    const volScalarField& H = db().lookupObject<volScalarField>(HName_);

    scalar surfaceSumH = gSum(H.boundaryField()[patch().index()]*patch().magSf());
//    Info << "surfaceSumH: " << surfaceSumH << endl;

    scalar Hmean = surfaceSumH/gSum(patch().magSf());
//        Info << "Hmean: " << Hmean << endl;
    scalar Umean = sqrt(2*9.81*(mag(energieHoehe_-Hmean)));

//    Info << "HU " << -sign(energieHoehe_-Hmean)*n*Hmean*Umean << endl;


    operator==(-sign(energieHoehe_-Hmean)*n*Hmean*Umean);

    fixedValueFvPatchField<vector>::updateCoeffs();
}


void Foam::energieHoeheInletFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);

    os.writeKeyword("energieHoehe") << energieHoehe_
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
       energieHoeheInletFvPatchVectorField
   );
}


// ************************************************************************* //
