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

#include "sifOmegaFvPatchScalarField.H"
#include "sifFlowdepthFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "sifMappedPatchBase.H"
#include "regionProperties.H"
#include "SortableList.H"
#include "Time.H"
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sifOmegaFvPatchScalarField::
sifOmegaFvPatchScalarField
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


Foam::sifOmegaFvPatchScalarField::
sifOmegaFvPatchScalarField
(
    const sifOmegaFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    neighbourFieldName_(ptf.neighbourFieldName_)
{}


Foam::sifOmegaFvPatchScalarField::
sifOmegaFvPatchScalarField
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
            "sifOmegaFvPatchScalarField::"
            "sifOmegaFvPatchScalarField\n"
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


Foam::sifOmegaFvPatchScalarField::
sifOmegaFvPatchScalarField
(
    const sifOmegaFvPatchScalarField& wtcsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(wtcsf, iF),
    neighbourFieldName_(wtcsf.neighbourFieldName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::scalarField 
Foam::sifOmegaFvPatchScalarField::baseArea(const polyMesh& localMesh, const fvPatch& nbrPatch, const mapDistribute& distMap) const
{
    scalarField baseAreaField(nbrPatch.size(),0);

    // Get the fvMesh from polyMesh. Only fvMesh contains the necessary geometrical data (i.e. faceAreas())
    const fvMesh& locFvMesh = refCast<const fvMesh>(localMesh);
    
    forAll(baseAreaField, i)
    {
        // Find the 3d-face which is closest to this 2d-face
        const label face3d = distMap.subMap()[Pstream::myProcNo()][i];
        // Get the 3d-cell belonging to the 3d-face
        const label cell3d = patch().faceCells()[face3d];
        // Get all faces of the 3d-cell
        const labelList cell3dFaces = localMesh.cells()[cell3d];

        scalar maxNormalZ(0);
        
        // Loop over all faces to find the one with the largest normal-component in z-direction
        forAll(cell3dFaces, j)
        {
	    scalar normalZ = mag(locFvMesh.faceAreas()[cell3dFaces[j]][2]);
            maxNormalZ = max(maxNormalZ,normalZ);
        }
      
        baseAreaField[i] = maxNormalZ;
    }
    
    return baseAreaField;
}

Foam::vectorField 
Foam::sifOmegaFvPatchScalarField::flowdepth3d(const polyMesh& nbrMesh, const fvPatch& nbrPatch, const mapDistribute& localMap, const scalarField& baseAreaField) const
{
    
    /*
    Concerning the [0] in the following statement, see the following thread:
     http://www.cfd-online.com/Forums/openfoam-programming-development/95076-access-patch-points-different-processor-parallel.html
    */
    const labelList& faces3d = localMap.subMap()[Pstream::myProcNo()];
    scalarField integralAlpha1(nbrPatch.size(),0);
    scalarField integralVolume(nbrPatch.size(),0);
    scalarField bottomCentreField(nbrPatch.size(),GREAT);
    scalarField bottomVolumeField(nbrPatch.size(),GREAT);
    
    const fvPatchScalarField& alphaFld = patch().lookupPatchField<volScalarField, scalar>("alpha1");
    const scalarField alphaIntFld = alphaFld.patchInternalField();
//    Info << "alphaIntFld: " << alphaIntFld << endl;
    
    forAll(faces3d, i)
    {
        const label face2d = localMap.subMap()[Pstream::myProcNo()][i];
        const label cell3d = patch().faceCells()[i];
//        Info << "cell3d: " << cell3d << endl;
        integralAlpha1[face2d] +=  patch().boundaryMesh().mesh().V()[cell3d] * alphaIntFld[i];
        integralVolume[face2d] +=  patch().boundaryMesh().mesh().V()[cell3d];
        
        if(patch().boundaryMesh().mesh().C()[cell3d][2] < bottomCentreField[face2d])
        {
            bottomCentreField[face2d] = patch().boundaryMesh().mesh().C()[cell3d][2];
            bottomVolumeField[face2d] = patch().boundaryMesh().mesh().V()[cell3d];
        }
    }
    
//    scalarField bottoms = bottomCentreField - 0.5 * bottomVolumeField / baseAreaField;
//    Info << "bottoms: " << bottoms << endl;
    
    vectorField depthsHeightsBottoms(nbrPatch.size(),Foam::vector(0.,0.,0.));
    
    forAll(depthsHeightsBottoms, i)
    {
        // -- Absolute flowdepth
        depthsHeightsBottoms[i][0] = integralAlpha1[i] / baseAreaField[i] + bottomCentreField[i] - 0.5 * bottomVolumeField[i] / baseAreaField[i];
        // Height of the cell column
        depthsHeightsBottoms[i][1] = integralVolume[i] / baseAreaField[i];
        // Absolute level of the bottom of the cell column
        depthsHeightsBottoms[i][2] = bottomCentreField[i] - 0.5 * bottomVolumeField[i] / baseAreaField[i];
    }
    
//    Info << "depthsHeightsBottoms: " << depthsHeightsBottoms << endl;    
    return depthsHeightsBottoms;

}

const Foam::vectorField 
Foam::sifOmegaFvPatchScalarField::wallShearStresses(const polyMesh& localMesh, const fvPatch& nbrPatch, const mapDistribute& localMap, const mapDistribute& distMap) const
{
    vectorField wssField(nbrPatch.size(),vector(0,0,0));
    
    labelList innerList;
    
    // Get the fvMesh from polyMesh. Only fvMesh contains the necessary geometrical data (i.e. faceAreas())
    const fvMesh& locFvMesh = refCast<const fvMesh>(localMesh);
    
    const labelList localSize(patch().size(),0);
    
    scalarField minZ(nbrPatch.size(),GREAT);
    labelList lowestFace(nbrPatch.size());
    
    forAll(localSize, i)
    {
        // Get the 3d-cell belonging to the 3d-face
        const label cell3d = patch().faceCells()[i];
        //Info << "cell3d: " << cell3d << endl;
        // Get the 2d-face belonging to the 3d-face
        const label face2d = localMap.subMap()[Pstream::myProcNo()][i];        
        //Info << "face2d: " << face2d << endl;

        // Get all faces of the 3d-cell
        innerList = localMesh.cells()[cell3d];
        //Info << "innerList: " << innerList << endl;
        
        forAll(innerList,j)
            {
               scalar tempZ = locFvMesh.faceCentres()[innerList[j]][2];
               //Info << "tempZ: " << tempZ << endl;

            if(tempZ < minZ[face2d])
                {
                    minZ[face2d] = locFvMesh.faceCentres()[innerList[j]][2];
                    lowestFace[face2d] = innerList[j];
                }
           }

     }
    
    forAll(lowestFace,k)
    {
        //- Return patch index for a given face label
        label patchID = localMesh.boundaryMesh().whichPatch(lowestFace[k]);
        const polyPatch& bottomPatch = localMesh.boundaryMesh()[patchID];
        // Get the fvMesh from polyMesh. Only fvMesh contains the necessary geometrical data (i.e. faceAreas())
        //const wallFvPatch& bottomFvPatch = refCast<const wallFvPatch>(bottomPatch); 
        label localPatchFace = bottomPatch.whichFace(lowestFace[k]);
        
        const fvPatch& bottomFvPatch = refCast<const fvMesh>
        (
            locFvMesh
        ).boundary()[patchID];
        
        symmTensor Reff = bottomFvPatch.lookupPatchField<volSymmTensorField, tensor>("Reff")[localPatchFace];
	//        scalar nutFace = bottomFvPatch.lookupPatchField<volScalarField, scalar>("nut")[localPatchFace];
        vector normalFaceVector = -bottomFvPatch.Sf()[localPatchFace] / bottomFvPatch.magSf()[localPatchFace];
        wssField[k]= Reff & normalFaceVector;
    }
    
//    Info << "minZ: " << minZ << endl;
//    Info << "lowestFace: " << lowestFace << endl;
//    Info << "wssField: " << wssField << endl;
    
    return wssField;
}



void Foam::sifOmegaFvPatchScalarField::updateCoeffs()
{
    
    if (this->updated())
    {
        return;
    }

    // Get the coupling information from the sifMappedPatchBase
    const sifMappedPatchBase& mpp = refCast<const sifMappedPatchBase>
    (
        patch().patch()
    );
    
    const polyMesh& localMesh = this->patch().patch().boundaryMesh().mesh();
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

    // Calculation of base area of local cells
    scalarField baseAreaField = baseArea(localMesh, nbrPatch, distNbrMap);

    // Calculation of flow depth based on alpha1 in 3D region
    vectorField intDepthsHeightsBottoms = flowdepth3d(localMesh, nbrPatch, localMap, baseAreaField);
    scalarField absDepthIntField(nbrPatch.size());
    scalarField heightIntField(nbrPatch.size());
    scalarField bottomIntField(nbrPatch.size());
    scalarField depthIntField(nbrPatch.size());

    forAll(intDepthsHeightsBottoms, i)
    {
        absDepthIntField[i] = intDepthsHeightsBottoms[i][0];
        heightIntField[i] = intDepthsHeightsBottoms[i][1];
        bottomIntField[i] = intDepthsHeightsBottoms[i][2];
    }

    depthIntField = absDepthIntField - bottomIntField;

    // 
    const fvPatchVectorField& dischargeNbrField = nbrPatch.lookupPatchField<volVectorField, vector>("HU");
    const vectorField dischargeNbrIntField = dischargeNbrField.patchInternalField();
    const fvPatchScalarField& depthNbrField = nbrPatch.lookupPatchField<volScalarField, scalar>("H");
    const scalarField depthNbrIntField = depthNbrField.patchInternalField();
    const fvPatchScalarField& bottomNbrField = nbrPatch.lookupPatchField<volScalarField, scalar>("S");
    const scalarField bottomNbrIntField = bottomNbrField.patchInternalField();

    const scalarField absDepthNbrIntField = bottomNbrIntField + depthNbrIntField;
    
    // List of neighbouring faces of local patch faces
    const labelList& nbrFaces = mpp.map().subMap()[Pstream::myProcNo()];
//    Info << "nbrFaces: " << nbrFaces << endl;
    
    // Vectors from cell centres to face centres
    const vectorField localDeltas = patch().delta();
    const vectorField nbrDeltas = nbrPatch.delta();

    // Calculation of wall shear stress at the bottom
    vectorField wssField = wallShearStresses(localMesh, nbrPatch, localMap, distNbrMap);
    
    const fvPatchVectorField& velField = patch().lookupPatchField<volVectorField, vector>("U");
    const vectorField velIntField = velField.patchInternalField();
    
    const fvPatchScalarField& alpha1PatchField = patch().lookupPatchField<volScalarField, scalar>("alpha1");
    const scalarField alpha1PatchIntField = alpha1PatchField.patchInternalField();    
    
    vectorField integralVel(nbrPatch.size(),Foam::vector(0.,0.,0.));
    scalarField integralFaceArea(nbrPatch.size(),0);
    
    // Calculation of average velocity in 3D region
    forAll(velIntField, i)
    {
        const label face2d = localMap.subMap()[Pstream::myProcNo()][i];
        integralVel[face2d] += velIntField[i] * patch().boundaryMesh().mesh().V()[i] * alpha1PatchIntField[i];
        integralFaceArea[face2d] += patch().boundaryMesh().mesh().V()[i] * alpha1PatchIntField[i];
    }
    
    vectorField vel3d = integralVel / (integralFaceArea + 0.0001);

    // Calculate the friction velocity
    vectorField uStar(nbrPatch.size(),Foam::vector(0.,0.,0.));
   
    forAll(wssField,i)
    {
        // WSS is already normalized by the density, so not necessary here anymore
        // Sign comes back in here
        // WSS ist schon mit Dichte normiert, deswegen hier nicht mehr drin.
        // Vorzeichen kommt hier wieder rein.
        uStar[i][0] = pow(mag(wssField[i][0]),0.5) * sign(vel3d[i][0]);
        uStar[i][1] = pow(mag(wssField[i][1]),0.5) * sign(vel3d[i][1]);
//        uStar[i][2] = pow(mag(wssField[i][2]),0.5) * sign(velNbrIntField[i][2]);
    }
    
    //    scalarField sAvField = calcSAv(absDepthIntField, bottomIntField, localMap);
    
    scalarField alpha1NbrIntField(patch().size(),0);
    
    const vectorField faceNormalField = patch().Sf() / patch().magSf();
    //const vectorField faceCenterField = patch().Cf();

    forAll(alpha1NbrIntField, i)
    {
        const label face2d = localMap.subMap()[Pstream::myProcNo()][i];

        scalar magNbr = mag(nbrDeltas[face2d]);
        vector localDeltaHoriz = localDeltas[i];
        // Vertical component set to zero, s.t. only horizontal components are used for linear interpolation on both sides.
        localDeltaHoriz[2] = 0.0;
        scalar magLoc = mag(localDeltaHoriz);

        // Calculation of averaged velocity
        vector vel2d = dischargeNbrIntField[face2d] / (depthNbrIntField[face2d] + 0.0001);
        vel3d[face2d][2] = 0.0;
        vector averageVel = (vel2d * magNbr + vel3d[face2d] * magLoc) / (magLoc + magNbr);
        bool inOrOut = pos(averageVel & faceNormalField[i]);

        // In case of outflow, set zeroGradient
        if(inOrOut == 1)
        {
	  this->refGrad()[i] = 0.0;
          this->valueFraction()[i] = 0.0;
        }

        // In case of inflow, set profile for omega
        else
        {
	  scalar averageDepth = ((depthIntField[face2d] * magNbr + depthNbrIntField[face2d] * magLoc) / (magLoc + magNbr)) + 0.0001;
	  scalar averageBottom = ( bottomIntField[face2d] * magNbr + bottomNbrIntField[face2d] * magLoc) / (magLoc + magNbr);
	  //	  scalar zFaceRel = (faceCenterField[i][2] - averageBottom) / averageDepth;
	  scalar zFaceRel = (patch().Cf()[i][2] - averageBottom) / averageDepth;
	  this->refValue()[i] = max(6. * mag(uStar[face2d]) * pow(zFaceRel,-1.3) / averageDepth,1);
	  this->valueFraction()[i] = 1.0;
	  // The following avoids division by 0, when there is no uStar yet; i.e. when there is inflow on dry bed)
	  if(this->refValue()[i] == 0.0)
	    {
	      this->valueFraction()[i] = 0.0;
	    }
        }
   
}

    /*    
    Info<< "refValue: " << this->refValue() << endl;
    Info<< "refGrad: " << this->refGrad() << endl;
    Info<< "valueFraction: " << this->valueFraction() << endl;
    */

    mixedFvPatchScalarField::updateCoeffs();

}

    
void Foam::sifOmegaFvPatchScalarField::write
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
    sifOmegaFvPatchScalarField
);

} // End namespace Foam

// ************************************************************************* //
