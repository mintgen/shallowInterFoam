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

\*---------------------------------------------------------------------------*/

#include "sifDischargeFvPatchVectorField.H"
#include "sifAlpha1FvPatchScalarField.H"
#include "sifFlowdepthFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "sifMappedPatchBase.H"
#include "regionProperties.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sifDischargeFvPatchVectorField::
sifDischargeFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchVectorField(p, iF),
    neighbourFieldName_("undefined-neighbourFieldName")
{
    this->refValue() = vector(0.0,0.0,0.0);
    this->refGrad() = vector(0.0,0.0,0.0);
    this->valueFraction() = 0.0;
}


Foam::sifDischargeFvPatchVectorField::
sifDischargeFvPatchVectorField
(
    const sifDischargeFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchVectorField(ptf, p, iF, mapper),
    neighbourFieldName_(ptf.neighbourFieldName_)
{
}


Foam::sifDischargeFvPatchVectorField::
sifDischargeFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchVectorField(p, iF),
    neighbourFieldName_(dict.lookup("neighbourFieldName"))
{
    if (!isA<sifMappedPatchBase>(this->patch().patch()))
    {
        FatalErrorIn
        (
            "sifDischargeFvPatchVectorField::"
            "sifDischargeFvPatchVectorField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<vector, volMesh>& iF,\n"
            "    const dictionary& dict\n"
            ")\n"
        )   << "\n    patch type '" << p.type()
            << "' not type '" << sifMappedPatchBase::typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << dimensionedInternalField().name()
            << " in file " << dimensionedInternalField().objectPath()
            << exit(FatalError);
    }

    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));

    if (dict.found("refValue"))
    {
        // Full restart
        refValue() = vectorField("refValue", dict, p.size());
        refGrad() = vectorField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = vector(0.0,0.0,0.0);
        valueFraction() = 0.0;
    }
}


Foam::sifDischargeFvPatchVectorField::
sifDischargeFvPatchVectorField
(
    const sifDischargeFvPatchVectorField& wtcsf,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchVectorField(wtcsf, iF),
    neighbourFieldName_(wtcsf.neighbourFieldName_)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


const Foam::scalarField 
Foam::sifDischargeFvPatchVectorField::baseArea(const polyMesh& nbrMesh, const fvPatch& nbrPatch, const mapDistribute& localMap, const mapDistribute& distMap) const
{
    scalarField baseAreaField(patch().size(),0);
    
    // Get the fvMesh from polyMesh. Only fvMesh contains the necessary geometrical data (i.e. faceAreas())
    const fvMesh& nbrFvMesh = refCast<const fvMesh>(nbrMesh);
    
    forAll(baseAreaField, i)
    {
        // Find a 3d-face belonging to this 2d-face
        const label face3d = localMap.subMap()[Pstream::myProcNo()][i];
        // Get the 3d-cell belonging to the 3d-face
        const label cell3d = nbrPatch.faceCells()[face3d];
        // Get all faces of the 3d-cell
        const labelList cell3dFaces = nbrMesh.cells()[cell3d];

        scalar maxNormalZ(0);
        
        // Loop over all faces to find the one with the largest normal-component in z-direction
        forAll(cell3dFaces, j)
        {
            scalar normalZ = mag(nbrFvMesh.faceAreas()[cell3dFaces[j]][2]);
            maxNormalZ = max(maxNormalZ,normalZ);
        }
      
        baseAreaField[i] = maxNormalZ;
    }
    
    return baseAreaField;
}


Foam::vectorField 
Foam::sifDischargeFvPatchVectorField::discharge3d(const polyMesh& nbrMesh, const fvPatch& nbrPatch, const mapDistribute& distMap, const scalarField baseAreaField) const
{
    
    /*
    Concerning the [0] in the following statement, see the following thread:
     http://www.cfd-online.com/Forums/openfoam-programming-development/95076-access-patch-points-different-processor-parallel.html
    */
    const labelList& faces3d = distMap.subMap()[Pstream::myProcNo()];
    vectorField integralVelocity(patch().size(),Foam::vector(0.,0.,0.));

    const fvPatchScalarField& alphaFld = nbrPatch.lookupPatchField<volScalarField, scalar>("alpha1");
    const scalarField alphaIntFld = alphaFld.patchInternalField();
//    Info << "alphaIntFld: " << alphaIntFld << endl;    
    
    const fvPatchVectorField& velocityFld = nbrPatch.lookupPatchField<volVectorField, vector>("U");
    const vectorField velocityIntFld = velocityFld.patchInternalField();
//    Info << "velocityIntFld: " << velocityIntFld << endl;
    
    forAll(faces3d, i)
    {
        const label face2d = distMap.subMap()[Pstream::myProcNo()][i];
        const label cell3d = nbrPatch.faceCells()[i];
        integralVelocity[face2d] +=  nbrPatch.boundaryMesh().mesh().V()[cell3d] * alphaIntFld[i] * velocityIntFld[i];
        integralVelocity[face2d][2] = 0.0;
    }
    
    vectorField discharges(patch().size(),Foam::vector(0.,0.,0.));
    discharges = integralVelocity / baseAreaField;

//    Info << "discharges: " << discharges << endl;    
    return discharges;

}
 
void Foam::sifDischargeFvPatchVectorField::updateCoeffs()
{

  //    Info << "Enter sifDischarge " << endl;
    
    if (updated())
    {
        return;
    }

    // Get the coupling information from the sifMappedPatchBase
    const sifMappedPatchBase& mpp = refCast<const sifMappedPatchBase>
    (
        patch().patch()
    );
    
    const polyMesh& nbrMesh = mpp.sampleMesh();
        
    const fvPatch& nbrPatch = refCast<const fvMesh>
    (
        nbrMesh
    ).boundary()[mpp.samplePolyPatch().index()];
    
    // Get the coupling information from the sifMappedPatchBase for the neighbour patch
    const sifMappedPatchBase& mppNbr = refCast<const sifMappedPatchBase>
    (
        nbrPatch.patch()
    );
        
    // Force recalculation of mapping and schedule
    const mapDistribute& localMap = mpp.map();
    const mapDistribute& distNbrMap = mppNbr.map();

    // Calculation fo base area of neighbour cells
    scalarField baseAreaField = baseArea(nbrMesh, nbrPatch, localMap, distNbrMap);

    // Calculation of discharge in 3D region
    vectorField huNbrField = discharge3d(nbrMesh, nbrPatch, distNbrMap, baseAreaField);


    // *** Uncomment the following block and modify this->refValue()[i] below, in order to get a linear interpolation of the discharges from both sides
   
    const fvPatchVectorField& huField = patch().lookupPatchField<volVectorField, vector>("HU");
    const vectorField huIntField = huField.patchInternalField();
    
    // List of neighbouring faces of local patch faces
    const labelList& nbrFaces = mpp.map().subMap()[Pstream::myProcNo()];
    
    // Vectors from cell centres to face centres
    const vectorField localDeltas = patch().delta();
    const vectorField nbrDeltas = nbrPatch.delta();
    
    vectorField huPatchField(patch().size(),Foam::vector(0.,0.,0.));
    
    forAll(huPatchField, i)
    {
        scalar magLoc = mag(localDeltas[i]);
        vector nbrDeltaHoriz = nbrDeltas[nbrFaces[i]];
        // Vertical component set to zero, s.t. only horizontal components are used for linear interpolation on both sides.
        nbrDeltaHoriz[2] = 0.0;
        scalar magNbr = mag(nbrDeltaHoriz);

        // Linear interpolation
        huPatchField[i] = (huIntField[i] * magNbr + huNbrField[i] * magLoc) / (magLoc + magNbr);

        // HU solely from 3D
        //huPatchField[i] = huNbrField[i];
    }

    const vectorField faceNormalField = patch().Sf() / patch().magSf();

    forAll(huNbrField, i)
    {
        const label face2d = localMap.subMap()[Pstream::myProcNo()][i];
        
        bool inOrOut = pos(huPatchField[i] & faceNormalField[i]);
        // In case of outflow, set zeroGradient
        if(inOrOut == 1)
        {
            this->valueFraction()[i] = 0.0;
        }
        // In case of inflow, set discharge
        else
        {
            this->valueFraction()[i] = 1.0;
            this->refValue()[i] = huNbrField[i];
            this->refValue()[i][2] = 0.0;
        }
}

    this->refGrad() = vector(0.0,0.0,0.0);
    
    mixedFvPatchVectorField::updateCoeffs();

    //    Info << "Exit sifDischarge " << endl;


}

void Foam::sifDischargeFvPatchVectorField::write
(
    Ostream& os
) const
{
    mixedFvPatchVectorField::write(os);
    os.writeKeyword("neighbourFieldName")<< neighbourFieldName_
        << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

makePatchTypeField
(
    fvPatchVectorField,
    sifDischargeFvPatchVectorField
);

} // End namespace Foam

// ************************************************************************* //
