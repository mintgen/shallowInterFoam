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

#include "sifAlpha1FvPatchScalarField.H"
#include "sifFlowdepthFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "sifMappedPatchBase.H"
#include "regionProperties.H"
#include "SortableList.H"
#include "Time.H"
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sifAlpha1FvPatchScalarField::
sifAlpha1FvPatchScalarField
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


Foam::sifAlpha1FvPatchScalarField::
sifAlpha1FvPatchScalarField
(
    const sifAlpha1FvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    neighbourFieldName_(ptf.neighbourFieldName_)
{}


Foam::sifAlpha1FvPatchScalarField::
sifAlpha1FvPatchScalarField
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
            "sifAlpha1FvPatchScalarField::"
            "sifAlpha1FvPatchScalarField\n"
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


Foam::sifAlpha1FvPatchScalarField::
sifAlpha1FvPatchScalarField
(
    const sifAlpha1FvPatchScalarField& wtcsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(wtcsf, iF),
    neighbourFieldName_(wtcsf.neighbourFieldName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// Function baseArea has been copied from sifFlowdepthFvPatchScalarField and just local and distant references have been interchanged. So probably this could be done in a more simple way...
const Foam::scalarField 
Foam::sifAlpha1FvPatchScalarField::baseArea(const polyMesh& localMesh, const fvPatch& nbrPatch, const mapDistribute& distMap) const
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



// For function flowdepth3d, the same holds as for function baseArea above.
Foam::vectorField 
Foam::sifAlpha1FvPatchScalarField::flowdepth3d(const polyMesh& nbrMesh, const fvPatch& nbrPatch, const mapDistribute& localMap, const scalarField& baseAreaField) const
{
    
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

Foam::vectorField 
Foam::sifAlpha1FvPatchScalarField::ghostCellCentres(const polyMesh& nbrMesh, const fvPatch& nbrPatch, const mapDistribute& localMap) const
{
    
    const vectorField deltasLoc = patch().delta();
    //    Info << "deltasLoc: " << deltasLoc << endl;
    const vectorField deltasNbr = nbrPatch.delta();
    //    Info << "deltasNbr: " << deltasNbr << endl;
    
    // Just for testing purposes:
    const fvMesh& nbrFvMesh = refCast<const fvMesh>(nbrMesh);
    // Test's end
    
    vectorField ghostCellCentreField(patch().size(),Foam::vector(0.,0.,0.));
    
    forAll(deltasLoc, i)
    {
        const label face2d = localMap.subMap()[Pstream::myProcNo()][i];
        vector centre2df = patch().Cf()[i];
	//        Info << "centre2df: " << centre2df << endl;
        vector deltaLocHor = deltasLoc[i];
        deltaLocHor[2] = 0.0;
	//        Info << "deltaLocHor: " << mag(deltaLocHor) << endl;
        ghostCellCentreField[i] = centre2df + deltasLoc[i] / mag(deltaLocHor) * mag(deltasNbr[face2d]);
    }

    //    Info << "ghostCellCentreField: " << ghostCellCentreField << endl;

    return ghostCellCentreField;
}

/*
Foam::scalarField 
Foam::sifAlpha1FvPatchScalarField::calcSAv(const scalarField& depthField, const scalarField& bottomField, const mapDistribute& localMap) const
{

    const fvPatchScalarField& alpha1Patch = patch().lookupPatchField<volScalarField, scalar>("alpha1");
    const scalarField alpha1Int = alpha1Patch.patchInternalField();

    scalarField sAvField(bottomField.size(),0);
    
    scalar argsCell = 0.0;
    scalar sCell = 0.0;
    scalarField sCellSum(bottomField.size(),0.0);
    scalarField countsCell(bottomField.size(),0);
    scalar epsS = 0.001;
    label face2d;
    label cell3d;
    
    
    scalar currentTime = this->db().time().timeOutputValue();

    string currTimeString = string(currentTime);
    


//    dataFile << "time" << " " << "face2d" << " " << "cell3d" << " " << "a1Int" << " " << "cellCZ" << " " << "depthField[face2d]" << endl;
    
    forAll(alpha1Int,i)
    {
        face2d = localMap.subMap()[Pstream::myProcNo()][i];
        cell3d = patch().faceCells()[i];
	scalar a1Int = alpha1Int[i];
	scalar cellCZ = patch().boundaryMesh().mesh().C()[cell3d].component(vector::Z);
	scalar cellCZRel = cellCZ - bottomField[face2d];

        
//        string prefix = "sAv";
//        int face2dInt = face2d;
//        string face2dString = string(face2dScalar);
//        string preFace = prefix + face2dString;
//        
//        Info << "face2dString: " << face2dString << endl;
//        Info << "preFace: " << preFace << endl;

        
//	if((0+epsS)<a1Int && a1Int<(1-epsS) && a1Int!=0.5 && cellCZ!=depthField[face2d])
	if((0+epsS)<a1Int && a1Int<(1-epsS) && a1Int!=0.5 && mag(cellCZ-depthField[face2d])>epsS)
	  {
	    argsCell = (a1Int-0.5)/0.5;
	    sCell = Foam::atanh(argsCell) / (cellCZ-depthField[face2d]);
	    sCellSum[face2d] += sCell;
	    countsCell[face2d] += 1;
            
//            OFstream dataFile
//            (
//                "preFace", std::ofstream::out | std::ofstream::app
//            );
//
//            dataFile << currentTime << " " << face2d << " " << cell3d << " " << a1Int << " " << cellCZ << " " << depthField[face2d] << " " << sCell << endl;
//            
//            Info << "bottomField: " << bottomField[face2d] << endl; 
//	    Info << "a1Int: " << a1Int << endl;
//            Info << "argsCell: " << argsCell << endl;
//            Info << "cellCZ: " << cellCZ << endl;
//            Info << "cellCZRel: " << cellCZRel << endl;
//            Info << "depthField[face2d]: " << depthField[face2d] << endl;
//	    Info << "sCell: " << sCell << "\n" << endl;

          }
        
    }
    
    sAvField = -150;

    forAll(countsCell,i)
    {
        if(countsCell[i] >= 1)
        {
            // Smooth interface -> smooth tanh-function
            sAvField[i] = sCellSum[i] / countsCell[i];
            // Limit interface extension to sAv = ...
            sAvField[i] = min(-25., sAvField[i]);
        }        
    }
    
    
    
//    Info << "sAvField alpha1: " << sAvField << endl;

    return sAvField;
   
}

*/

void Foam::sifAlpha1FvPatchScalarField::updateCoeffs()
{

  //    Info << "Entering sifAlpha1" << endl;

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
    
    // Calculation of height of local cells
    scalarField cellHeightField(patch().size(),0.);
    vectorField cellCenterField(patch().size(),Foam::vector(0.,0.,0.));
    
    forAll(cellHeightField, i)
    {
        const label face2d = localMap.subMap()[Pstream::myProcNo()][i];
        const label cell3d = patch().faceCells()[i];

        cellHeightField[i] = patch().boundaryMesh().mesh().V()[cell3d] / baseAreaField[face2d];
        cellCenterField[i] = patch().boundaryMesh().mesh().C()[cell3d];
    }     

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

    // Calculation of ghost cell centres
    vectorField ghostCellCs = ghostCellCentres(nbrMesh, nbrPatch, localMap);
    
//    Info << "ghostCellCs: " << ghostCellCs << endl;

    // 
    const fvPatchVectorField& dischargeNbrField = nbrPatch.lookupPatchField<volVectorField, vector>("HU");
    const vectorField dischargeNbrIntField = dischargeNbrField.patchInternalField();
    const fvPatchScalarField& depthNbrField = nbrPatch.lookupPatchField<volScalarField, scalar>("H");
    const scalarField depthNbrIntField = depthNbrField.patchInternalField();
    const fvPatchScalarField& bottomNbrField = nbrPatch.lookupPatchField<volScalarField, scalar>("S");
    const scalarField bottomNbrIntField = bottomNbrField.patchInternalField();

    const scalarField absDepthNbrField = bottomNbrField + depthNbrField;
    const scalarField absDepthNbrIntField = bottomNbrIntField + depthNbrIntField;
    
    // List of neighbouring faces of local patch faces
    const labelList& nbrFaces = mpp.map().subMap()[Pstream::myProcNo()];
    
    // Vectors from cell centres to face centres
    const vectorField localDeltas = patch().delta();
    const vectorField nbrDeltas = nbrPatch.delta();
    
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

    // + 0.0001 to avoid division by zero, in case of a flowdepth equal to 0
    vectorField vel3d = integralVel / (integralFaceArea + 0.0001);
    
    //    scalarField sAvField = calcSAv(absDepthIntField, bottomIntField, localMap);
    
    //scalarField alpha1NbrIntField(patch().size(),0);
    //scalarField alpha1IntFieldReconstructed(patch().size(),0);
    
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
    const scalarField depthNbrIntFieldClip = (depthNbrIntField-Hdry)*pos(depthNbrIntField-Hdry) + Hdry;

    forAll(alpha1PatchIntField, i)
    {

        const label face2d = localMap.subMap()[Pstream::myProcNo()][i];

        scalar magNbr = mag(nbrDeltas[face2d]);
        vector localDeltaHoriz = localDeltas[i];
        // Vertical component set to zero, s.t. only horizontal components are used for linear interpolation on both sides.
        localDeltaHoriz[2] = 0.0;
        scalar magLoc = mag(localDeltaHoriz);

        // Calculation of patch-normal Froude-Number: (HU & n) / ( 9.81^0.5 * H^1.5 )
        //scalar froude = (dischargeNbrField[face2d] & faceNormalField[i]) / (3.132092 * pow(depthNbrField[face2d],1.5) );

        // Calculation of averaged patch-normal Froude number
        vector vel2d = dischargeNbrIntField[face2d] / depthNbrIntFieldClip[face2d];
        vel3d[face2d][2] = 0.0;
        vector averageVel = (vel3d[face2d] * magNbr + vel2d * magLoc) / (magLoc + magNbr);
        scalar averageDepth = (depthIntFieldClip[face2d] * magNbr + depthNbrIntFieldClip[face2d] * magLoc) / (magLoc + magNbr);
        scalar froude = (averageVel & faceNormalField[i]) / (3.132092 * pow(averageDepth,0.5) );

        // Neumann condition for subcritical inflow
        if(-1. <= froude && froude <= 0.)
        {
	  // Reconstruction of alpha1-distribution in local cells, based on flow depth (This is done in order to avoid diffusion of the interface)
	  scalar centerCellLocal = cellCenterField[i][2];
	  scalar distCtoWSPLocal = centerCellLocal - absDepthIntField[face2d];

	  scalar centerCell = ghostCellCs[i][2];
	  scalar distCtoWSP = centerCell - absDepthNbrIntField[face2d];

	  scalar distCtoWSPPatch = (distCtoWSPLocal * magNbr + distCtoWSP * magLoc) / (magLoc + magNbr);

	  scalar alpha1Local = -0.5 * (distCtoWSPPatch / (cellHeightField[i]/2.) ) + 0.5;
	  
	  this->valueFraction()[i] = 1.0;
	  this->refValue()[i] = max(0.0,min(1.0,alpha1Local)); 
	  //Calculation of gradient
	  //this->refGrad()[i] = (alpha1NbrIntField[i] - alpha1PatchIntField[i]) / (magLoc + magNbr);

          this->valueFraction()[i] = 0.0;
	  this->refGrad()[i] = 0.0;  
           
        }

        // Neumann condition for supercritical outflow
        else if(froude > 1.)
        {
            this->valueFraction()[i] = 0.0;
	    this->refGrad()[i] = 0.0;            
            //Calculation of gradient
	    //this->refGrad()[i] = (alpha1NbrIntField[i] - alpha1PatchIntField[i]) / (magLoc + magNbr);            
        }

        // Dirichlet condition for supercritical inflow
        else if(froude < -1.)
        {
	  
	  // Calculation of alpha1-distribution in ghost cells        
	  scalar centerCell = ghostCellCs[i][2];
	  scalar distCtoWSP = centerCell - absDepthNbrIntField[face2d];
	  //scalar alpha1GhostCell = 0.5*Foam::tanh(sAvField[face2d]*distCtoWSP) + 0.5;
	  scalar alpha1GhostCell = -0.5 * (distCtoWSP / (cellHeightField[i]/2.) ) + 0.5;
	  //alpha1NbrIntField[i] = max(0.0,min(1.0,alpha1GhostCell));

	  this->valueFraction()[i] = 1.0;
	  this->refValue()[i] = max(0.0,min(1.0,alpha1GhostCell));
	  
	  /*
	  // Calculation of alpha1-distribution in boundary faces        
	  scalar centerFace = patch().Cf()[i][2];
	  scalar distCtoWSP = centerFace - absDepthNbrField[face2d];
	  //scalar alpha1GhostCell = 0.5*Foam::tanh(sAvField[face2d]*distCtoWSP) + 0.5;
	  scalar alpha1BoundaryFace = -0.5 * (distCtoWSP / (cellHeightField[i]/2.) ) + 0.5;
	  //alpha1NbrIntField[i] = max(0.0,min(1.0,alpha1GhostCell));

	  this->valueFraction()[i] = 1.0;
	  this->refValue()[i] = max(0.0,min(1.0,alpha1BoundaryFace));
	  */
        }

        // Neumann condition for subcritical outflow
        else if(0. < froude && froude <= 1)
        {
            this->valueFraction()[i] = 0.0;
	    this->refGrad()[i] = 0.0; 
            //Calculation of gradient
	    //this->refGrad()[i] = (alpha1NbrIntField[i] - alpha1PatchIntField[i]) / (magLoc + magNbr);

	    /*
            // Linear interpolation
	    scalar centerCell = ghostCellCs[i][2];
	    scalar distCtoWSP = centerCell - absDepthNbrIntField[face2d];
	    scalar alpha1GhostCell = -0.5 * (distCtoWSP / (cellHeightField[i]/2.) ) + 0.5;
	    alpha1GhostCell = max(0.0,min(1.0,alpha1GhostCell));
	    scalar newAlpha1Patch = (alpha1PatchIntField[i] * magNbr + alpha1GhostCell * magLoc) / (magLoc + magNbr);
	    this->refValue()[i] = max(0.0,min(1.0,newAlpha1Patch)); 
	    this->valueFraction()[i] = 1.0;
	    */

        }

        else
        {
            FatalErrorIn("sifAlpha1FvPatchScalarField.C");
        }
        
    

}

    //    Info << "Leaving sifAlpha1" << endl;

    mixedFvPatchScalarField::updateCoeffs();

}

    
void Foam::sifAlpha1FvPatchScalarField::write
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
    sifAlpha1FvPatchScalarField
);

} // End namespace Foam

// ************************************************************************* //
