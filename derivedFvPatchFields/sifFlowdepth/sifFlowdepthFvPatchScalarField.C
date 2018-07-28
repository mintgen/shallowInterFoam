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

#include "sifFlowdepthFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "sifMappedPatchBase.H"
#include "regionProperties.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sifFlowdepthFvPatchScalarField::
sifFlowdepthFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    neighbourFieldName_("undefined-neighbourFieldName")
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
}


Foam::sifFlowdepthFvPatchScalarField::
sifFlowdepthFvPatchScalarField
(
    const sifFlowdepthFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    neighbourFieldName_(ptf.neighbourFieldName_)
{}


Foam::sifFlowdepthFvPatchScalarField::
sifFlowdepthFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    neighbourFieldName_(dict.lookup("neighbourFieldName"))
{
    if (!isA<sifMappedPatchBase>(this->patch().patch()))
    {
        FatalErrorIn
        (
            "sifFlowdepthFvPatchScalarField::"
            "sifFlowdepthFvPatchScalarField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<scalar, volMesh>& iF,\n"
            "    const dictionary& dict\n"
            ")\n"
        )   << "\n    patch type '" << p.type()
            << "' not type '" << sifMappedPatchBase::typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << dimensionedInternalField().name()
            << " in file " << dimensionedInternalField().objectPath()
            << exit(FatalError);
    }

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    if (dict.found("refValue"))
    {
        // Full restart
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = 0.0;
        valueFraction() = 0.0;
    }
}


Foam::sifFlowdepthFvPatchScalarField::
sifFlowdepthFvPatchScalarField
(
    const sifFlowdepthFvPatchScalarField& wtcsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(wtcsf, iF),
    neighbourFieldName_(wtcsf.neighbourFieldName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::scalarField 
Foam::sifFlowdepthFvPatchScalarField::baseArea(const polyMesh& nbrMesh, const fvPatch& nbrPatch, const mapDistribute& localMap) const
{
    scalarField baseAreaField(patch().size(),0);
    
    forAll(baseAreaField, i)
    {
        // Find a 3d-face belonging to this 2d-face
        const label face3d = localMap.subMap()[Pstream::myProcNo()][i];
        // Get the 3d-cell belonging to the 3d-face
        const label cell3d = nbrPatch.faceCells()[face3d];
        // Get all faces of the 3d-cell
        const labelList cell3dFaces = nbrMesh.cells()[cell3d];

        // Get the fvMesh from polyMesh. Only fvMesh contains the necessary geometrical data (i.e. faceAreas())
        const fvMesh& nbrFvMesh = refCast<const fvMesh>(nbrMesh);
        
        scalar maxNormalZ(0);
        
        // Loop over all faces to find the one with the largest normal-component in z-direction
        forAll(cell3dFaces, j)
        {
	    scalar normalZ = mag(nbrFvMesh.faceAreas()[cell3dFaces[j]][2]);
            maxNormalZ = max(maxNormalZ,normalZ);
        }
      
        baseAreaField[i] = maxNormalZ;
    }
    
    //    Info << "baseAreaField: " << baseAreaField << endl;    
    
    return baseAreaField;
}


Foam::vectorField 
Foam::sifFlowdepthFvPatchScalarField::flowdepth3d(const polyMesh& nbrMesh, const fvPatch& nbrPatch, const mapDistribute& distMap, const scalarField& baseAreaField) const
{

    const labelList& faces3d = distMap.subMap()[Pstream::myProcNo()];
    scalarField integralAlpha1(patch().size(),0);
    scalarField integralVolume(patch().size(),0);
    scalarField bottomCentreField(patch().size(),GREAT);
    scalarField bottomVolumeField(patch().size(),GREAT);

    const fvPatchScalarField& alphaFld = nbrPatch.lookupPatchField<volScalarField, scalar>("alpha1");
    const scalarField alphaIntFld = alphaFld.patchInternalField();
    
    forAll(faces3d, i)
    {
        const label face2d = distMap.subMap()[Pstream::myProcNo()][i];
        const label cell3d = nbrPatch.faceCells()[i];
        integralAlpha1[face2d] +=  nbrPatch.boundaryMesh().mesh().V()[cell3d] * alphaIntFld[i];
        integralVolume[face2d] +=  nbrPatch.boundaryMesh().mesh().V()[cell3d];

        if(nbrPatch.boundaryMesh().mesh().C()[cell3d][2] < bottomCentreField[face2d])
        {
            bottomCentreField[face2d] = nbrPatch.boundaryMesh().mesh().C()[cell3d][2];
            bottomVolumeField[face2d] = nbrPatch.boundaryMesh().mesh().V()[cell3d];
        }
    }
    
    vectorField depthsHeightsBottoms(patch().size(),Foam::vector(0.,0.,0.));

    forAll(depthsHeightsBottoms, i)
    {
        // -- Absolute flowdepth
        depthsHeightsBottoms[i][0] = integralAlpha1[i] / baseAreaField[i] + bottomCentreField[i] - 0.5 * bottomVolumeField[i] / baseAreaField[i];
        // Height of the cell column
        depthsHeightsBottoms[i][1] = integralVolume[i] / baseAreaField[i];
        // Absolute level of the bottom of the cell column
        depthsHeightsBottoms[i][2] = bottomCentreField[i] - 0.5 * bottomVolumeField[i] / baseAreaField[i];
    }

    return depthsHeightsBottoms;

}

Foam::vectorField 
Foam::sifFlowdepthFvPatchScalarField::discharge3d(const polyMesh& nbrMesh, const fvPatch& nbrPatch, const mapDistribute& distMap, const scalarField baseAreaField) const
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
 
void Foam::sifFlowdepthFvPatchScalarField::updateCoeffs()
{

  //    Info << "Enter sifFlowdepth " << endl;

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

/*    
    Info << "localMap.subMap(): " << localMap.subMap() << endl;
    Info << "localMap.constructMap(): " << localMap.constructMap() << endl;
    Info << "distNbrMap.subMap(): " << distNbrMap.subMap() << endl;
    Info << "distNbrMap.constructMap(): " << distNbrMap.constructMap() << endl;
*/
    
    // Calculation fo base area of neighbour cells
    scalarField baseAreaField = baseArea(nbrMesh, nbrPatch, localMap);

    // Calculation of flow depth based on alpha1 in 3D region
    vectorField nbrDepthsHeightsBottoms = flowdepth3d(nbrMesh, nbrPatch, distNbrMap, baseAreaField);
    scalarField absDepthNbrField(patch().size());
    scalarField heightNbrField(patch().size());
    scalarField bottomNbrField(patch().size());
    scalarField depthNbrField(patch().size());

    forAll(nbrDepthsHeightsBottoms, i)
    {
        absDepthNbrField[i] = nbrDepthsHeightsBottoms[i][0];
        heightNbrField[i] = nbrDepthsHeightsBottoms[i][1];
        bottomNbrField[i] = nbrDepthsHeightsBottoms[i][2];
    }

    depthNbrField = absDepthNbrField - bottomNbrField;

    // Calculation of flow depth based on alpha1 in 3D region
    vectorField dischargeNbrField = discharge3d(nbrMesh, nbrPatch, distNbrMap, baseAreaField);    

    //
    const fvPatchVectorField& dischargeField = patch().lookupPatchField<volVectorField, vector>("HU");
    const vectorField dischargeIntField = dischargeField.patchInternalField();
    const fvPatchScalarField& depthField = patch().lookupPatchField<volScalarField, scalar>("H");
    const scalarField depthIntField = depthField.patchInternalField();
    const fvPatchScalarField& bottomField = patch().lookupPatchField<volScalarField, scalar>("S");
    const scalarField bottomIntField = bottomField.patchInternalField();
    
    const scalarField absDepthIntField = depthIntField + bottomIntField;

    // List of neighbouring faces of local patch faces
    const labelList& nbrFaces = mpp.map().subMap()[Pstream::myProcNo()];
    
    // Vectors from cell centres to face centres
    const vectorField localDeltas = patch().delta();
    const vectorField nbrDeltas = nbrPatch.delta();
    
    //scalarField depthPatchField(patch().size(),0.);
    //scalarField gradPatchField(patch().size(),0.);
    
    const vectorField faceNormalField = patch().Sf() / patch().magSf();

    // The following steps and their usage later on are necessary in order to avoid a division by 0
    // Getting Hdry from the transportProperties 
    const dictionary& transportProperties = db().lookupObject<IOdictionary>
      (
       "transportProperties"
       );
    // Extracting scalar value
    dimensionedScalar HdryDim(transportProperties.lookup("Hdry"));
    const scalar Hdry = HdryDim.value();
    // depth...FieldClip = max(depth...Field,Hdry);
    const scalarField depthIntFieldClip = (depthIntField-Hdry)*pos(depthIntField-Hdry) + Hdry;
    const scalarField depthNbrFieldClip = (depthNbrField-Hdry)*pos(depthNbrField-Hdry) + Hdry;


    forAll(depthIntField, i)
    {
        scalar magLoc = mag(localDeltas[i]);
        vector nbrDeltaHoriz = nbrDeltas[nbrFaces[i]];
        // Vertical component set to zero, s.t. only horizontal components are used for linear interpolation on both sides.
        nbrDeltaHoriz[2] = 0.0;
        scalar magNbr = mag(nbrDeltaHoriz);

        // Linear interpolation
        //depthPatchField[i] = (depthIntField[i] * magNbr + depthNbrField[i] * magLoc) / (magLoc + magNbr);

	// H solely from 3D
        //depthPatchField[i] = depthNbrField[i];
	
	// Calculation of gradient
        //gradPatchField[i] = (depthNbrField[i] - depthIntField[i]) / (magLoc + magNbr);

	// Zero gradient condition
	//gradPatchField[i] = 0.0;

        // Calculation of averaged patch-normal Froude number
        vector vel2d = dischargeIntField[i] / depthIntFieldClip[i];
        vector vel3d = dischargeNbrField[i] / depthNbrFieldClip[i];
        vector averageVel = (vel2d * magNbr + vel3d * magLoc) / (magLoc + magNbr);
        scalar averageDepth = (depthIntFieldClip[i] * magNbr + depthNbrFieldClip[i] * magLoc) / (magLoc + magNbr);
        scalar froude = (averageVel & faceNormalField[i]) / (3.132092 * pow(averageDepth,0.5) );

        
        // Neumann condition for supercritical outflow
        if(froude > 1.)
        {
            this->valueFraction()[i] = 0.0;
            //this->refGrad()[i] = gradPatchField[i];
	    //            this->refGrad()[i] = (absDepthNbrField[i] - absDepthIntField[i]) / (magLoc + magNbr);
            this->refGrad()[i] = 0.0;
        }

        // Neumann condition for subcritical inflow
        else if(-1. <= froude && froude < 0.)
        {
            this->valueFraction()[i] = 0.0;
            //this->refGrad()[i] = gradPatchField[i];
            this->refGrad()[i] = (absDepthNbrField[i] - absDepthIntField[i]) / (magLoc + magNbr);
            this->refGrad()[i] = 0.0;
        }

        // Dirichlet condition for supercritical inflow
        else if(froude < -1.)
        {
            this->valueFraction()[i] = 1.0;       
	    // d h / d n = 0
            //this->refValue()[i] = depthNbrField[i];
	    // d h_abs / d n = 0
            // this->refValue()[i] = absDepthNbrField[i] - bottomField[i];
	    // d h / d n != 0
            this->refValue()[i] = (absDepthIntField[i] * magNbr + absDepthNbrField[i] * magLoc) / (magLoc + magNbr) - bottomField[i];
        }
        
        // Dirichlet condition for subcritical outflow
        else if(0. <= froude && froude <= 1)
        {
            this->valueFraction()[i] = 1.0;        
	    //this->refValue()[i] = depthNbrField[i];
            this->refValue()[i] = (absDepthIntField[i] * magNbr + absDepthNbrField[i] * magLoc) / (magLoc + magNbr) - bottomField[i];
        }
        
        else
        {
            FatalErrorIn("sifFlowdepthFvPatchScalarField.C");
        }
    }

    
    /*
    Info<< "refValue: " << this->refValue() << endl;
    Info<< "refGrad: " << this->refGrad() << endl;
    Info<< "valueFraction: " << this->valueFraction() << endl;
    */
    mixedFvPatchScalarField::updateCoeffs();

    //    Info << "Exit sifFlowdepth " << endl;

}


void Foam::sifFlowdepthFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    os.writeKeyword("neighbourFieldName")<< neighbourFieldName_
        << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

makePatchTypeField
(
    fvPatchScalarField,
    sifFlowdepthFvPatchScalarField
);

} // End namespace Foam

// ************************************************************************* //
