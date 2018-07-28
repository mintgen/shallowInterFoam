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

Class
    Foam::energieHoeheInletFvPatchVectorField

Description

    Setzt HU so dass gilt:
          energieHoehe = H + U^2/2g

    Example of the boundary condition specification:
    @verbatim
    inlet
    {
        type            energieHoeheInlet;
        energieHoehe        0.2;        // Energiehoehe [m]
        value           uniform (0 0 0); // placeholder
    }
    @endverbatim

Note
    - H am Rand muss bestimmt sein (z.B. zeroGradient)
    - Einstroemen normal zur Oberflaeche

Author
    KM-Turbulenz GmbH, 2009
    dervied from flowRateInletVelocityFvPatchVectorField.H

SourceFiles
    energieHoeheInletFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef energieHoeheInletFvPatchVectorField_H
#define energieHoeheInletFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
               Class energieHoeheInletFvPatch Declaration
\*---------------------------------------------------------------------------*/

class energieHoeheInletFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Inlet EnergieHoehe
        scalar energieHoehe_;

    //
    word HName_;


public:

   //- Runtime type information
   TypeName("energieHoeheInlet");


   // Constructors

        //- Construct from patch and internal field
        energieHoeheInletFvPatchVectorField
        (
	    const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        energieHoeheInletFvPatchVectorField
        (
	    const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  energieHoeheInletFvPatchVectorField
        //  onto a new patch
        energieHoeheInletFvPatchVectorField
        (
            const energieHoeheInletFvPatchVectorField&,
	    const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        energieHoeheInletFvPatchVectorField
        (
            const energieHoeheInletFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new energieHoeheInletFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        energieHoeheInletFvPatchVectorField
        (
            const energieHoeheInletFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new energieHoeheInletFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the flux
            scalar energieHoehe() const
            {
                return energieHoehe_;
            }

            //- Return reference to the flux to allow adjustment
            scalar& energieHoehe()
            {
                return energieHoehe_;
            }


        //- Update the coefficients associated with the patch field
//        virtual void updateCoeffs();
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
