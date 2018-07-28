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

#include "sifVelocityFvPatchVectorField.H"
#include "sifAlpha1FvPatchScalarField.H"
#include "sifFlowdepthFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "sifMappedPatchBase.H"
#include "regionProperties.H"
#include "RASModel.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sifVelocityFvPatchVectorField::
sifVelocityFvPatchVectorField
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


Foam::sifVelocityFvPatchVectorField::
sifVelocityFvPatchVectorField
(
    const sifVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchVectorField(ptf, p, iF, mapper),
    neighbourFieldName_(ptf.neighbourFieldName_)
{}


Foam::sifVelocityFvPatchVectorField::
sifVelocityFvPatchVectorField
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
            "sifVelocityFvPatchVectorField::"
            "sifVelocityFvPatchVectorField\n"
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


Foam::sifVelocityFvPatchVectorField::
sifVelocityFvPatchVectorField
(
    const sifVelocityFvPatchVectorField& wtcsf,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchVectorField(wtcsf, iF),
    neighbourFieldName_(wtcsf.neighbourFieldName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Function baseArea has been copied from sifFlowdepthFvPatchScalarField and just local and distant references have been interchanged. So probably this could be done in a more simple way...
const Foam::scalarField 
Foam::sifVelocityFvPatchVectorField::baseArea(const polyMesh& localMesh, const fvPatch& nbrPatch, const mapDistribute& distMap) const
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

const Foam::vectorField 
Foam::sifVelocityFvPatchVectorField::wallShearStresses(const polyMesh& localMesh, const fvPatch& nbrPatch, const mapDistribute& localMap, const mapDistribute& distMap) const
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
        // Get the 2d-face belonging to the 3d-face
        const label face2d = localMap.subMap()[Pstream::myProcNo()][i];        

        // Get all faces of the 3d-cell
        innerList = localMesh.cells()[cell3d];
        
        forAll(innerList,j)
            {
               scalar tempZ = locFvMesh.faceCentres()[innerList[j]][2];

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


	// Get the sand grain roughness Ks from the nut-patch-dictionary
	//	const fvPatchScalarField& nutField = bottomFvPatch.lookupPatchField<volScalarField, scalar>("nut");

/*        
        const scalarField& localKs = nutField.lookupPatchField
        (
            "Ks",
            reinterpret_cast<const surfaceScalarField*>(0),
            reinterpret_cast<const scalar*>(0)
        );        
*/

//        nutField.writeEntry("Ks", Info);

//        Info << "nutField:" << nutField << endl;
//        Info << "localKs:" << localKs << endl;

/*        
        Info << "patchName:" << localMesh.boundaryMesh()[patchID].name() << endl;
        Info << "localPatchFace:" << localPatchFace << endl;
        Info << "Reff:" << Reff << endl;
        Info << "nutFace:" << nutFace << endl;
        Info << "normalFaceVector:" << normalFaceVector << endl;
        Info << "wssField[k]:" << wssField[k] << endl;
*/        
    }
    
//    Info << "minZ: " << minZ << endl;
//    Info << "lowestFace: " << lowestFace << endl;
//    Info << "wssField: " << wssField << endl;

    return wssField;
}

Foam::vectorField 
Foam::sifVelocityFvPatchVectorField::ghostCellCentres(const polyMesh& nbrMesh, const fvPatch& nbrPatch, const mapDistribute& localMap) const
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

Foam::scalarField 
Foam::sifVelocityFvPatchVectorField::calcAlpha1Nbr(const vectorField& localDeltas, const vectorField& nbrDeltas, const labelList& nbrFaces) const
{
    
    scalarField alpha1Nbr(patch().size(),0.0);
    
    const fvPatchScalarField& alpha1PatchField = patch().lookupPatchField<volScalarField, scalar>("alpha1");
    const scalarField alpha1PatchIntField = alpha1PatchField.patchInternalField();    

    forAll(alpha1Nbr, i)
    {
        scalar magNbr = mag(nbrDeltas[nbrFaces[i]]);
        vector localDeltaHoriz = localDeltas[i];
        // Vertical component set to zero, s.t. only horizontal components are used for linear interpolation on both sides.
        localDeltaHoriz[2] = 0.0;
        scalar magLoc = mag(localDeltaHoriz);
        
        alpha1Nbr[i] = (alpha1PatchField[i] * (magLoc + magNbr) - alpha1PatchIntField[i] * magNbr) / magLoc;

        // Linear interpolation
//        velPatchField[i] = (velIntField[i] * magNbr + velNbrIntField[nbrFaces[i]] * magLoc) / (magLoc + magNbr);
    }
    
//    Info << "alpha1NbrVel: " << alpha1Nbr << endl;
   
    return alpha1Nbr;
}

Foam::vectorField 
Foam::sifVelocityFvPatchVectorField::ghostCellVelocities(const polyMesh& nbrMesh, const fvPatch& nbrPatch, const mapDistribute& localMap, const mapDistribute& distMap, vectorField& wssField, vectorField& ghostCellCs, scalarField& alpha1NbrField) const
{

    scalar eps = 0.0001;
    
    vectorField ghostCellVelocityField(patch().size(),Foam::vector(0.,0.,0.));
    
    const fvPatchVectorField& huNbrField = nbrPatch.lookupPatchField<volVectorField, vector>("HU");    
    const fvPatchScalarField& hNbrField = nbrPatch.lookupPatchField<volScalarField, scalar>("H");
    const fvPatchScalarField& SNbrField = nbrPatch.lookupPatchField<volScalarField, scalar>("S");    
    vectorField huNbrIntField = huNbrField.patchInternalField();
    scalarField hNbrIntField = hNbrField.patchInternalField();
    scalarField SNbrIntField = SNbrField.patchInternalField();
    
    vectorField velNbrIntField(nbrPatch.size(),Foam::vector(0.,0.,0.));
    
    velNbrIntField = huNbrIntField / (hNbrIntField + eps);
     
    // Calculate the friction velocity
    vectorField uStar(nbrPatch.size(),Foam::vector(0.,0.,0.));
   
    forAll(wssField,i)
    {
        // WSS is already normalized by the density, so not necessary here anymore
        // Sign comes back in here
        // WSS ist schon mit Dichte normiert, deswegen hier nicht mehr drin.
        // Vorzeichen kommt hier wieder rein.
        uStar[i][0] = pow(mag(wssField[i][0]),0.5) * sign(velNbrIntField[i][0]);
        uStar[i][1] = pow(mag(wssField[i][1]),0.5) * sign(velNbrIntField[i][1]);
//        uStar[i][2] = pow(mag(wssField[i][2]),0.5) * sign(velNbrIntField[i][2]);
    }
   
//    Info << "velNbrIntField: " << velNbrIntField << endl;
//    Info << "uStar: " << uStar << endl;
    
    const scalar kappa = 0.41;
    const scalar PI = 0.2;
    const scalar pi = mathematicalConstant::pi;
    const scalar B2 = 8.5;
    // Sandrauheit, muss noch softcoded werden
    //const scalar ks = 0.001;
    
    vectorField integralVel(nbrPatch.size(),Foam::vector(0.,0.,0.));
    scalarField integralFaceArea(nbrPatch.size(),0);
     
    forAll(ghostCellVelocityField, i)
    {
        const label face2d = localMap.subMap()[Pstream::myProcNo()][i];
        
        scalar relGhostCellHeight = ghostCellCs[i][2] - SNbrIntField[face2d];

	if(relGhostCellHeight < 0)
	  {
	    Info << "Check if S-field is valid!" << endl;
	  }

	ghostCellVelocityField[i] = velNbrIntField[face2d] 
	  + uStar[face2d]/kappa * (1 + Foam::log(relGhostCellHeight / (hNbrIntField[face2d] + eps)));
	//     - PI*uStar[face2d]/kappa * pow(sin( pi*relGhostCellHeight / (2.*hNbrIntField[face2d]) ),2.);
                //* pos(alpha1NbrField[i]-eps);

	// Geschwindigkeitsprofil nach Pope, Gl. 7.119, 7.148 & 7.149
	//	ghostCellVelocityField[i] = uStar[face2d] * (Foam::log(relGhostCellHeight / (ks * hNbrIntField[face2d])) / kappa + B2);

//                             Velocity wake, not in agreement with results of numerical experiments
//	           	            + 2. * PI / kappa * pow(sin( pi*relGhostCellHeight / (2.*hNbrIntField[face2d]) ),2.));

        // Geschwindigkeiten nahe der Wand mit falschem VZ eliminieren:
        ghostCellVelocityField[i][0] = max(0,(ghostCellVelocityField[i][0]*sign(velNbrIntField[face2d][0]))) * sign(velNbrIntField[face2d][0]);
        ghostCellVelocityField[i][1] = max(0,(ghostCellVelocityField[i][1]*sign(velNbrIntField[face2d][1]))) * sign(velNbrIntField[face2d][1]);
        
        // z-Komponente gleich 0 setzen
        ghostCellVelocityField[i][2] = 0.;
        
        integralVel[face2d] += ghostCellVelocityField[i] * patch().magSf()[i] * alpha1NbrField[i];
        integralFaceArea[face2d] += patch().magSf()[i] * alpha1NbrField[i];
    }
  
    vectorField averageVel(nbrPatch.size(),Foam::vector(0.,0.,0.));

	/*
    forAll(ghostCellVelocityField, i)
    {
        const label face2d = localMap.subMap()[Pstream::myProcNo()][i];
    
        if(integralFaceArea[face2d]==0)
        {
            averageVel[face2d] = Foam::vector(0.,0.,0.);
        }
        else
        {
            averageVel[face2d] = integralVel[face2d] / integralFaceArea[face2d];
        }
    }
	*/

    averageVel = integralVel / (integralFaceArea + eps);

    // correctorField is used to ensure conservativity
    vectorField correctorField(nbrPatch.size(),Foam::vector(0.,0.,0.));
    scalar maxMagCorrector = 0.;

    forAll(correctorField, face2d)
    {
        /*
        correctorField[face2d][0] = (velNbrIntField[face2d][0] - averageVel[face2d][0]) / velNbrIntField[face2d][0];
        correctorField[face2d][1] = (velNbrIntField[face2d][1] - averageVel[face2d][1]) / velNbrIntField[face2d][0];
        maxMagCorrector = max(maxMagCorrector,mag(correctorField[face2d]));
        */
        if(averageVel[face2d][0]==0)
        {
            correctorField[face2d][0] = 1.;
        }
        else
        {
            correctorField[face2d][0] = velNbrIntField[face2d][0] / averageVel[face2d][0];
        }
        
        if(averageVel[face2d][1]==0)
        {
            correctorField[face2d][1] = 1.;
        }
        else
        {
            correctorField[face2d][1] = velNbrIntField[face2d][1] / averageVel[face2d][1];
        }

        maxMagCorrector = max(maxMagCorrector, 1. - mag(correctorField[face2d]));
    }

    while(maxMagCorrector > eps)
    {
        integralVel = vector(0.,0.,0.);
        
        forAll(ghostCellVelocityField, i)
        {
            const label face2d = localMap.subMap()[Pstream::myProcNo()][i];
            ghostCellVelocityField[i][0] = ghostCellVelocityField[i][0] * correctorField[face2d][0];
            ghostCellVelocityField[i][1] = ghostCellVelocityField[i][1] * correctorField[face2d][1];
            integralVel[face2d] += ghostCellVelocityField[i] * patch().magSf()[i] * alpha1NbrField[i];
        }
        
//        Info << "ghostCellVelocityField: " << ghostCellVelocityField << endl;
        
        averageVel = integralVel / (integralFaceArea + eps);
/*        
        Info << "averageVel: " << averageVel << endl;
        Info << "integralVel: " << integralVel << endl;
        Info << "integralFaceArea: " << integralFaceArea << endl;
        Info << "velNbrIntField: " << velNbrIntField << endl; 
*/        
        maxMagCorrector = 0.;
        
        forAll(correctorField, face2d)
        {
            
            if(averageVel[face2d][0]==0)
            {
                correctorField[face2d][0] = 1.;
            }
            else
            {
                correctorField[face2d][0] = velNbrIntField[face2d][0] / averageVel[face2d][0];
            }
        
            if(averageVel[face2d][1]==0)
            {
                correctorField[face2d][1] = 1.;
            }
            else
            {
                correctorField[face2d][1] = velNbrIntField[face2d][1] / averageVel[face2d][1];
            }            

            maxMagCorrector = max(maxMagCorrector, 1. - mag(correctorField[face2d]));
        }
        
//        Info << "maxMagCorrector: " << maxMagCorrector << endl;
    }


    return ghostCellVelocityField;
}

// For function flowdepth3d, the same holds as for function baseArea above.
Foam::vectorField 
Foam::sifVelocityFvPatchVectorField::flowdepth3d(const polyMesh& nbrMesh, const fvPatch& nbrPatch, const mapDistribute& localMap, const scalarField& baseAreaField) const
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
    
    scalarField bottoms = bottomCentreField - 0.5 * bottomVolumeField / baseAreaField;
//    Info << "bottoms: " << bottoms << endl;
    
    vectorField depthsHeightsBottoms(nbrPatch.size(),Foam::vector(0.,0.,0.));
    
    forAll(depthsHeightsBottoms, i)
    {
        depthsHeightsBottoms[i][0] = integralAlpha1[i] / baseAreaField[i] + bottomCentreField[i] - 0.5 * bottomVolumeField[i] / baseAreaField[i];
        depthsHeightsBottoms[i][1] = integralVolume[i] / baseAreaField[i];
        depthsHeightsBottoms[i][2] = bottomVolumeField[i] / baseAreaField[i];
    }
    
//    Info << "depthsHeightsBottoms: " << depthsHeightsBottoms << endl;    
    return depthsHeightsBottoms;

}

void Foam::sifVelocityFvPatchVectorField::updateCoeffs()
{

  //    Info<< "Entering sifVelocity" << endl;
    
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

    // Calculation of flow depth based on alpha1 in 3D region & of height of 3D cells
    vectorField intDepthsHeightsBottoms = flowdepth3d(localMesh, nbrPatch, localMap, baseAreaField);
    scalarField depthIntField(nbrPatch.size());
    scalarField heightIntField(nbrPatch.size());
    scalarField bottomIntField(nbrPatch.size());

    forAll(intDepthsHeightsBottoms, i)
    {
        depthIntField[i] = intDepthsHeightsBottoms[i][0];
        heightIntField[i] = intDepthsHeightsBottoms[i][1];
        bottomIntField[i] = intDepthsHeightsBottoms[i][2];
    }

    
    // List of neighbouring faces of local patch faces
    const labelList& nbrFaces = mpp.map().subMap()[Pstream::myProcNo()];
    
    // Vectors from cell centres to face centres
    const vectorField localDeltas = patch().delta();
    const vectorField nbrDeltas = nbrPatch.delta();

    vectorField velPatchField(patch().size(),Foam::vector(0.,0.,0.));

    // Calculation of wall shear stress at the bottom
    vectorField wssField = wallShearStresses(localMesh, nbrPatch, localMap, distNbrMap);

    // Calculation of ghost cell centres
    vectorField ghostCellCs = ghostCellCentres(nbrMesh, nbrPatch, localMap);

    // Calculation of neighbour alpha1-Field (by extrapolation from patch)
    scalarField alpha1NbrField = calcAlpha1Nbr(localDeltas, nbrDeltas, nbrFaces);

    // Calculation of ghost cell velocities
    vectorField ghostCellVels = ghostCellVelocities(nbrMesh, nbrPatch, localMap, distNbrMap, wssField, ghostCellCs, alpha1NbrField);  

    const fvPatchVectorField& velField = patch().lookupPatchField<volVectorField, vector>("U");
    const vectorField velIntField = velField.patchInternalField();

    const fvPatchScalarField& alpha1Field = patch().lookupPatchField<volScalarField, scalar>("alpha1");

    scalarField valueFractionField(patch().size(),0.);
    
    vectorField integralVel(nbrPatch.size(),Foam::vector(0.,0.,0.));
    vectorField maxVel(nbrPatch.size(),Foam::vector(0.,0.,0.));
    scalarField integralFaceArea(nbrPatch.size(),0);

    forAll(velIntField, i)
    {

        const label face2d = localMap.subMap()[Pstream::myProcNo()][i];
        
        // Interpolation of velocities on faces
        scalar magNbr = mag(nbrDeltas[face2d]);
//        scalar magNbr = mag(nbrDeltas[nbrFaces[i]]);
        vector localDeltaHoriz = localDeltas[i];
        // Vertical component set to zero, s.t. only horizontal components are used for linear interpolation on both sides.
        localDeltaHoriz[2] = 0.0;
        scalar magLoc = mag(localDeltaHoriz);

        // Linear interpolation
        //vector newVelPatch = (velIntField[i] * magNbr + ghostCellVels[i] * magLoc) / (magLoc + magNbr);

	// Finding the largest velocity in the water phase
	if(mag(ghostCellVels[i]*alpha1Field[i]) > mag(maxVel[face2d]))
	{
	  maxVel[face2d] = ghostCellVels[i];
	}
        
        integralVel[face2d] += ghostCellVels[i] * patch().magSf()[i] * alpha1Field[i];
        integralFaceArea[face2d] += patch().magSf()[i] * alpha1Field[i];
    }

    vectorField averageVel = integralVel / (integralFaceArea + 0.00001);

    const vectorField faceNormalField = patch().Sf() / patch().magSf();

    forAll(velIntField, i)
    {
        const label face2d = localMap.subMap()[Pstream::myProcNo()][i];
        
        bool inOrOut = pos(averageVel[face2d] & faceNormalField[i]);

        // In case of outflow, set zeroGradient
        if(inOrOut == 1)
        {
            this->valueFraction()[i] = 0.0;
	    this->refGrad()[i] = vector(0.0,0.0,0.0);

	    /*
	    // Modify velocity of air phase to stabilize solution
	    if(alpha1Field[i] < 0.01)
	      {
		this->valueFraction()[i] = 1.0;
		// Zero gradient with limiter for all velocity components
		//		this->refValue()[i][0] = max(-1.0,min(1.0,velIntField[i][0]));
		//		this->refValue()[i][1] = max(-1.0,min(1.0,velIntField[i][1]));
		//		this->refValue()[i][2] = max(-1.0,min(1.0,velIntField[i][2]));
		// All velocity components set to zero
		this->refValue()[i] = vector(0.0,0.0,0.0);
	      }
	    */
        }

        // In case of inflow, set profile for velocity
        else
        {
            this->valueFraction()[i] = 1.0;
	    this->refValue()[i] = ghostCellVels[i];

	    // Set velocity of air phase to maximum vel. of water phase
	    if(alpha1Field[i] < 0.01)
	      {
		this->refValue()[i] = maxVel[face2d];
	      }

	    // Zero gradient for vertical velocity component
            this->refValue()[i][2] = velIntField[i][2];
	    // Vertical velocity component set to zero
            //this->refValue()[i][2] = 0.0;

	    /*
	    // If the flow depth in the 2D region is smaller than 0.001 m, set a zero-gradient for the velocity
	    if(depthIntField[face2d] < 0.001)
	      {
		this->valueFraction()[i] = 0.0;
		this->refGrad()[i] = vector(0.0,0.0,0.0);	
	      }
	    */
        }
}

    //    Info<< "Leaving sifVelocity" << endl;
    
    mixedFvPatchVectorField::updateCoeffs();

}

void Foam::sifVelocityFvPatchVectorField::write
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
    sifVelocityFvPatchVectorField
);

} // End namespace Foam

// ************************************************************************* //
