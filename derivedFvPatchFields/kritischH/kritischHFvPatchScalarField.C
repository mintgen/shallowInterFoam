/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "kritischHFvPatchScalarField.H"
#include "freestreamFvPatchFields.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

kritischHFvPatchScalarField::kritischHFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    kritischeFliesstiefe_(0)
{}


kritischHFvPatchScalarField::kritischHFvPatchScalarField
(
    const kritischHFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    kritischeFliesstiefe_(ptf.kritischeFliesstiefe_)    
{}


kritischHFvPatchScalarField::kritischHFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict),
    kritischeFliesstiefe_(readScalar(dict.lookup("kritischeFliesstiefe")))
{}


kritischHFvPatchScalarField::kritischHFvPatchScalarField
(
    const kritischHFvPatchScalarField& wbppsf
)
:
    fixedValueFvPatchField<scalar>(wbppsf),
    kritischeFliesstiefe_(wbppsf.kritischeFliesstiefe_)
{}


kritischHFvPatchScalarField::kritischHFvPatchScalarField
(
    const kritischHFvPatchScalarField& wbppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(wbppsf, iF),
    kritischeFliesstiefe_(wbppsf.kritischeFliesstiefe_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void kritischHFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
 
//  give internal values of H next to boundary patch
    scalarField Hint = patch().lookupPatchField<volScalarField, scalar>("H").patchInternalField();
    scalarField Sp = patch().lookupPatchField<volScalarField, scalar>("S");
    
    //    operator==((Hint-kritischeFliesstiefe_)*pos(Hint-kritischeFliesstiefe_) + kritischeFliesstiefe_);
    operator==((Hint-(kritischeFliesstiefe_-Sp))*pos(Hint-(kritischeFliesstiefe_-Sp)) + (kritischeFliesstiefe_-Sp));

//    operator==(Hrand);

    fixedValueFvPatchScalarField::updateCoeffs();
}

void kritischHFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);

    os.writeKeyword("kritischeFliesstiefe") << kritischeFliesstiefe_
        << token::END_STATEMENT << nl;

    writeEntry("value", os);
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, kritischHFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
