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

#include "sifMappedPolyPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "objectRegistry.H"
#include "Time.H"
#include "polyMesh.H"
#include "polyBoundaryMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sifMappedPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, sifMappedPolyPatch, word);
    addToRunTimeSelectionTable(polyPatch, sifMappedPolyPatch, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sifMappedPolyPatch::sifMappedPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm
)
:
    polyPatch(name, size, start, index, bm),
    sifMappedPatchBase(static_cast<const polyPatch&>(*this))
{}


Foam::sifMappedPolyPatch::sifMappedPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const word& sampleRegion,
    const sifMappedPatchBase::sampleMode mode,
    const word& samplePatch,
    const vectorField& offset,
    const polyBoundaryMesh& bm
)
:
    polyPatch(name, size, start, index, bm),
    sifMappedPatchBase
    (
        static_cast<const polyPatch&>(*this),
        sampleRegion,
        mode,
        samplePatch,
        offset
    )
{}


Foam::sifMappedPolyPatch::sifMappedPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const word& sampleRegion,
    const sifMappedPatchBase::sampleMode mode,
    const word& samplePatch,
    const vector& offset,
    const polyBoundaryMesh& bm
)
:
    polyPatch(name, size, start, index, bm),
    sifMappedPatchBase
    (
        static_cast<const polyPatch&>(*this),
        sampleRegion,
        mode,
        samplePatch,
        offset
    )
{}


Foam::sifMappedPolyPatch::sifMappedPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm
)
:
    polyPatch(name, dict, index, bm),
    sifMappedPatchBase(*this, dict)
{}


Foam::sifMappedPolyPatch::sifMappedPolyPatch
(
    const sifMappedPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    polyPatch(pp, bm),
    sifMappedPatchBase(*this, pp)
{}


Foam::sifMappedPolyPatch::sifMappedPolyPatch
(
    const sifMappedPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    polyPatch(pp, bm, index, newSize, newStart),
    sifMappedPatchBase(*this, pp)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sifMappedPolyPatch::~sifMappedPolyPatch()
{
    sifMappedPatchBase::clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Initialise the calculation of the patch geometry
void Foam::sifMappedPolyPatch::initGeometry()
{
    polyPatch::initGeometry();
}

//- Calculate the patch geometry
void Foam::sifMappedPolyPatch::calcGeometry()
{
    polyPatch::calcGeometry();
}

//- Initialise the patches for moving points
void Foam::sifMappedPolyPatch::initMovePoints(const pointField& p)
{
    polyPatch::initMovePoints(p);

    // Force recalculation of mapping with new point position
    // Note: this uses parallel communications.  HJ, 13/Mar/2012
    sifMappedPatchBase::clearOut();
    sifMappedPatchBase::map();
}

//- Correct patches after moving points
void Foam::sifMappedPolyPatch::movePoints(const pointField& p)
{
    polyPatch::movePoints(p);
}

//- Initialise the update of the patch topology
void Foam::sifMappedPolyPatch::initUpdateMesh()
{
    polyPatch::initUpdateMesh();

    // Force recalculation of mapping with new point position
    // Note: this uses parallel communications.  HJ, 13/Mar/2012
    sifMappedPatchBase::clearOut();

    // Only carry out mapping if the sampled region has been created already
    // DC, 04/Nov/2013
    if (boundaryMesh().mesh().time().found(sampleRegion()))
    {
        sifMappedPatchBase::map();
    }
}

//- Update of the patch topology
void Foam::sifMappedPolyPatch::updateMesh()
{
    polyPatch::updateMesh();
}


void Foam::sifMappedPolyPatch::write(Ostream& os) const
{
    polyPatch::write(os);
    sifMappedPatchBase::write(os);
}


// ************************************************************************* //
