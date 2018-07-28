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
    sifFlowdepthFvPatchScalarField

Description
    Mixed boundary condition for temperature, to be used by the
    conjugate heat transfer solver.
    Both sides use a mix of zeroGradient and neighbourvalue.

    Example usage:
        myInterfacePatchName
        {
            type                sifFlowdepth;
            neighbourFieldName  alpha1;
            K                   K;
            value               uniform 300;
        }

    Needs to be on underlying sifMapped(Wall)FvPatch.

    Note: runs in parallel with arbitrary decomposition. Uses sifMapped
    functionality to calculate exchange.

SourceFiles
    sifFlowdepthFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef sifFlowdepthFvPatchScalarField_H
#define sifFlowdepthFvPatchScalarField_H

#include "sifMappedPatchBase.H"
#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class sifFlowdepthFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class sifFlowdepthFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private data

        //- Name of field on the neighbour region
        const word neighbourFieldName_;

public:

    //- Runtime type information
    TypeName("sifFlowdepth");



    // Constructors

        //- Construct from patch and internal field
        sifFlowdepthFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        sifFlowdepthFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  sifFlowdepthFvPatchScalarField onto a new patch
        sifFlowdepthFvPatchScalarField
        (
            const sifFlowdepthFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new sifFlowdepthFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        sifFlowdepthFvPatchScalarField
        (
            const sifFlowdepthFvPatchScalarField&,
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
                new sifFlowdepthFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

    
        //- Calculate the base area of neighbour cells
        const scalarField baseArea(const polyMesh&, const fvPatch&, const mapDistribute&) const;

        //- Calculate the flowdepth
        vectorField flowdepth3d(const polyMesh&, const fvPatch&, const mapDistribute&, const scalarField&) const;
        
        //- Calculate the discharge in 3D region
        vectorField discharge3d(const polyMesh&, const fvPatch&, const mapDistribute&, const scalarField) const;        

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
