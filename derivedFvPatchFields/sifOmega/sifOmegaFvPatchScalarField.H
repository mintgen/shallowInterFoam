/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    sifOmegaFvPatchScalarField

Description
    Mixed boundary condition for temperature, to be used by the
    conjugate heat transfer solver.
    Both sides use a mix of zeroGradient and neighbourvalue.

    Example usage:
        myInterfacePatchName
        {
            type                sifOmega;
            neighbourFieldName  Omega;
            K                   K;
            value               uniform 300;
        }

    Needs to be on underlying sifMapped(Wall)FvPatch.

    Note: runs in parallel with arbitrary decomposition. Uses sifMapped
    functionality to calculate exchange.

SourceFiles
    sifOmegaFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef sifOmegaFvPatchScalarField_H
#define sifOmegaFvPatchScalarField_H

#include "sifMappedPatchBase.H"
#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class sifOmegaFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class sifOmegaFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private data

        //- Name of field on the neighbour region
        const word neighbourFieldName_;

public:

    //- Runtime type information
    TypeName("sifOmega");


    // Constructors

        //- Construct from patch and internal field
        sifOmegaFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        sifOmegaFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  sifOmegaFvPatchScalarField onto a new patch
        sifOmegaFvPatchScalarField
        (
            const sifOmegaFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new sifOmegaFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        sifOmegaFvPatchScalarField
        (
            const sifOmegaFvPatchScalarField&,
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
                new sifOmegaFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        //- Calculate the base area of local cells
        const scalarField baseArea(const polyMesh&, const fvPatch&, const mapDistribute&) const;

        //- Calculate the flowdepth
        vectorField flowdepth3d(const polyMesh&, const fvPatch&, const mapDistribute&, const scalarField&) const;

        //- Calculate the wall shear stresses
        const vectorField wallShearStresses(const polyMesh&, const fvPatch&, const mapDistribute&, const mapDistribute&) const;

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
