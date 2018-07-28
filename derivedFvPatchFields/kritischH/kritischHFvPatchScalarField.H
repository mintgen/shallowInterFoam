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

Class
    Foam::kritischHFvPatchScalarField

Description
    Foam::kritischHFvPatchScalarField

SourceFiles
    kritischHFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef kritischHFvPatchScalarFields_H
#define kritischHFvPatchScalarFields_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class kritischHFvPatch Declaration
\*---------------------------------------------------------------------------*/

class kritischHFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- Inlet EnergieHoehe
        scalar kritischeFliesstiefe_;

public:

    //- Runtime type information
    TypeName("kritischH");


    // Constructors

        //- Construct from patch and internal field
        kritischHFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        kritischHFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given kritischHFvPatchScalarField onto
        //  a new patch
        kritischHFvPatchScalarField
        (
            const kritischHFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        kritischHFvPatchScalarField
        (
            const kritischHFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new kritischHFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        kritischHFvPatchScalarField
        (
            const kritischHFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new kritischHFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the flux
            scalar kritischeFliesstiefe() const
            {
                return kritischeFliesstiefe_;
            }

            //- Return reference to the flux to allow adjustment
            scalar& kritischeFliesstiefe()
            {
                return kritischeFliesstiefe_;
            }
        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

        //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
