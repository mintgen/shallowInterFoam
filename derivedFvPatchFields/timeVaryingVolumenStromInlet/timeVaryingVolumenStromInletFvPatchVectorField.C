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

#include "timeVaryingVolumenStromInletFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::
timeVaryingVolumenStromInletFvPatchVectorField::
timeVaryingVolumenStromInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
//    volumenStrom_(0),
    timeSeries_(),
    HName_("H")
{}


Foam::
timeVaryingVolumenStromInletFvPatchVectorField::
timeVaryingVolumenStromInletFvPatchVectorField
(
    const timeVaryingVolumenStromInletFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
//    volumenStrom_(ptf.volumenStrom_),
    timeSeries_(ptf.timeSeries_),
    HName_(ptf.HName_)
{}


Foam::
timeVaryingVolumenStromInletFvPatchVectorField::
timeVaryingVolumenStromInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict),
//    volumenStrom_(readScalar(dict.lookup("volumenStrom"))),
    timeSeries_(dict),
    HName_("H")
{
   if (dict.found("value"))
   {
       fvPatchField<vector>::operator==(Field<vector>("value", dict, p.size()));
   }
   else
   {
       updateCoeffs();
   }

    if (dict.found("H"))
    {
        dict.lookup("H") >> HName_;
    }
}


Foam::
timeVaryingVolumenStromInletFvPatchVectorField::
timeVaryingVolumenStromInletFvPatchVectorField
(
    const timeVaryingVolumenStromInletFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    timeSeries_(ptf.timeSeries_),
    HName_(ptf.HName_)
{}


Foam::
timeVaryingVolumenStromInletFvPatchVectorField::
timeVaryingVolumenStromInletFvPatchVectorField
(
    const timeVaryingVolumenStromInletFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    timeSeries_(ptf.timeSeries_),
    HName_(ptf.HName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::timeVaryingVolumenStromInletFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }


    vectorField n = patch().nf();

    const volScalarField& H = db().lookupObject<volScalarField>(HName_);

    scalar volumenStrom = timeSeries_(this->db().time().timeOutputValue());

    const scalarField pointsZ = patch().patch().points().component(vector::Z);
    const scalar heightPatch = max(pointsZ) - min(pointsZ);    
    
    scalar surfaceSumH = gSum(H.boundaryField()[patch().index()]*patch().magSf());
    scalar Umean = volumenStrom/surfaceSumH * heightPatch;
    scalarField phip = patch().lookupPatchField<surfaceScalarField, scalar>("phi");
    scalar surfaceSumPhip = gSum(phip*H.boundaryField()[patch().index()]) / heightPatch;
    Info << "Inflow-Volumenstrom - Q: " <<  surfaceSumPhip << endl;   

    operator==(-n*Umean*H.boundaryField()[patch().index()]);

    fixedValueFvPatchField<vector>::updateCoeffs();
}


void Foam::timeVaryingVolumenStromInletFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);

    os.writeKeyword("timeSeries") << timeSeries_
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
       timeVaryingVolumenStromInletFvPatchVectorField
   );
}


// ************************************************************************* //
