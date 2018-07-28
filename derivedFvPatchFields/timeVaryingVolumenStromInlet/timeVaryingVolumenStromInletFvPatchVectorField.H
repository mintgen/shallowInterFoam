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
    Foam::volumenStromInletFvPatchVectorField

Description

    Setzt HU so dass gilt:
          energieHoehe = H + U^2/2g

    Example of the boundary condition specification:
    @verbatim
    inlet
    {
        type            volumenStromInlet;
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
    volumenStromInletFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef timeVaryingVolumenStromInletFvPatchVectorField_H
#define timeVaryingVolumenStromInletFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
               Class volumenStromInletFvPatch Declaration
\*---------------------------------------------------------------------------*/

class timeVaryingVolumenStromInletFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Inlet EnergieHoehe
        interpolationTable<scalar> timeSeries_;
//        scalar volumenStrom_;

    //
    word HName_;


public:

   //- Runtime type information
   TypeName("timeVaryingVolumenStromInlet");


   // Constructors

        //- Construct from patch and internal field
        timeVaryingVolumenStromInletFvPatchVectorField
        (
	    const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        timeVaryingVolumenStromInletFvPatchVectorField
        (
	    const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  volumenStromInletFvPatchVectorField
        //  onto a new patch
        timeVaryingVolumenStromInletFvPatchVectorField
        (
            const timeVaryingVolumenStromInletFvPatchVectorField&,
	    const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        timeVaryingVolumenStromInletFvPatchVectorField
        (
            const timeVaryingVolumenStromInletFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new timeVaryingVolumenStromInletFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        timeVaryingVolumenStromInletFvPatchVectorField
        (
            const timeVaryingVolumenStromInletFvPatchVectorField&,
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
                new timeVaryingVolumenStromInletFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the flux
//            scalar volumenStrom() const
              const interpolationTable<scalar>& timeSeries() const
            {
//                return volumenStrom_;
                return timeSeries_;
            }

            //- Return reference to the flux to allow adjustment
//            scalar& volumenStrom()
//            {
//                return volumenStrom_;
//            }


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
