/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.4.0
|    o     o     |  ENGYS Ltd. <http://engys.com/>
|       o        |
\*---------------------------------------------------------------------------
License
    This file is part of HELYXcore.
    HELYXcore is based on OpenFOAM (R) <http://www.openfoam.org/>.

    HELYXcore is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    HELYXcore is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with HELYXcore.  If not, see <http://www.gnu.org/licenses/>.

Copyright
    (c) 2024-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "nonConformal/polyPatches/nonConformalDiscreteMixing/nonConformalDiscreteMixingPolyPatch.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(nonConformalDiscreteMixingPolyPatch, 0);

    addToRunTimeSelectionTable
    (
        polyPatch,
        nonConformalDiscreteMixingPolyPatch,
        word
    );
    addToRunTimeSelectionTable
    (
        polyPatch,
        nonConformalDiscreteMixingPolyPatch,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::nonConformalDiscreteMixingPolyPatch::initCalcGeometry
(
    PstreamBuffers& pBufs
)
{
    directPolyPatch::initCalcGeometry(pBufs);

    intersectionsPtr_.clear();
    intersectionsAreValid_ = false;
}


void Foam::nonConformalDiscreteMixingPolyPatch::calcGeometry
(
    PstreamBuffers& pBufs
)
{
    calcSectors();

    setPatchTransforms();
}


void Foam::nonConformalDiscreteMixingPolyPatch::initTopoChange
(
    PstreamBuffers& pBufs
)
{
    directPolyPatch::initTopoChange(pBufs);

    intersectionsPtr_.clear();
    intersectionsAreValid_ = false;
}


void Foam::nonConformalDiscreteMixingPolyPatch::topoChange
(
    PstreamBuffers& pBufs
)
{
    directPolyPatch::topoChange(pBufs);
}


void Foam::nonConformalDiscreteMixingPolyPatch::rename(const wordList& newNames)
{
    directPolyPatch::rename(newNames);

    nbrPatch().nbrPatchName_ = newNames[index()];

    nonConformalCoupledPolyPatch::rename(newNames);
}


void Foam::nonConformalDiscreteMixingPolyPatch::reorder
(
    const labelUList& newToOldIndex
)
{
    directPolyPatch::reorder(newToOldIndex);

    if (nbrPatchID_ != -1)
    {
        nbrPatchID_ = findIndex(newToOldIndex, nbrPatchID_);
    }

    nonConformalCoupledPolyPatch::reorder(newToOldIndex);
}


void Foam::nonConformalDiscreteMixingPolyPatch::calcSectors()
{
    if (nSectors_ > 0) return;

    const scalar secAngle =
        owner() ? sectorAngles_.first() : sectorAngles_.second();

    if (mag(secAngle) > 360)
    {
        FatalErrorInFunction
            << "Invalid sector angle for patch " << name()
            << exit(FatalError);
    }

    scalar nSectors = 360.0/mag(secAngle);

    const scalar errIntSector = mag(nSectors - std::round(nSectors));

    if (errIntSector > tolIntSec_)
    {
        FatalErrorInFunction
            << "The number of sectors calculated for the non-conformal "
            << "discrete mixing patch " << name() << " is " << nSectors << ","
            << nl << "above the maximum allowed tolerance to be rounded to an "
            << "integer number of sectors. Please, check your geometry "
            << "and mesh generation."
            << exit(FatalError);
    }

    nSectors_ = std::round(nSectors);
}


void Foam::nonConformalDiscreteMixingPolyPatch::setPatchTransforms()
{
    if (patchTransforms_.size()) return;

    patchTransforms_.setSize(nSectors_);

    // Insert the initial (trivial) transform
    patchTransforms_[0] = transformer::I;

    // Then append all rotational transforms up to a full turn
    dictionary transformDict;
    transformDict.add("transformType", "rotational");
    transformDict.add("rotationAxis", rotationAxis_);
    transformDict.add("rotationCentre", rotationCentre_);

    const scalar secAngle =
        owner() ? mag(sectorAngles_.first()) : mag(sectorAngles_.second());

    scalar angle = secAngle;
    for (int i = 1; i < nSectors_; i++)
    {
        // Overwrite previous entry
        transformDict.add("rotationAngle", angle, true);

        patchTransforms_[i] = cyclicTransform(transformDict, true).transform();

        angle += secAngle;
    }
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::nonConformalDiscreteMixingPolyPatch::nonConformalDiscreteMixingPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    coupledPolyPatch(name, size, start, index, bm, patchType),
    nonConformalCoupledPolyPatch(static_cast<const polyPatch&>(*this)),
    nbrPatchName_(word::null),
    nbrPatchID_(-1),
    rotationAxis_(vector::uniform(NaN)),
    rotationCentre_(vector::uniform(NaN)),
    sectorAngles_(Pair<scalar>(NaN)),
    nSectors_(-1),
    patchTransforms_(),
    intersectionsAreValid_(false),
    intersectionsPtr_(nullptr)
{}


Foam::nonConformalDiscreteMixingPolyPatch::nonConformalDiscreteMixingPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    coupledPolyPatch(name, dict, index, bm, patchType),
    nonConformalCoupledPolyPatch(*this, dict),
    nbrPatchName_(dict.lookupOrDefault("neighbourPatch", word::null)),
    nbrPatchID_(-1),
    rotationAxis_(dict.lookup<vector>("rotationAxis")),
    rotationCentre_(dict.lookup<point>("rotationCentre")),
    sectorAngles_(dict.lookup<scalarList>("sectorAngles")),
    nSectors_(dict.lookup<label>("nSectors")),
    patchTransforms_(),
    intersectionsAreValid_(false),
    intersectionsPtr_(nullptr)
{}


Foam::nonConformalDiscreteMixingPolyPatch::nonConformalDiscreteMixingPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType,
    const word& nbrPatchName,
    const word& origPatchName,
    const dictionary& discMixingDict
)
:
    coupledPolyPatch(name, size, start, index, bm, patchType),
    nonConformalCoupledPolyPatch(*this, origPatchName),
    nbrPatchName_(nbrPatchName),
    nbrPatchID_(-1),
    rotationAxis_(discMixingDict.lookup<vector>("rotationAxis")),
    rotationCentre_(discMixingDict.lookup<point>("rotationCentre")),
    sectorAngles_
    (
        discMixingDict.lookup<scalar>("ownSectorAngle"),
        discMixingDict.lookup<scalar>("nbrSectorAngle")
    ),
    nSectors_(-1),
    patchTransforms_(),
    intersectionsAreValid_(false),
    intersectionsPtr_(nullptr)
{}


Foam::nonConformalDiscreteMixingPolyPatch::nonConformalDiscreteMixingPolyPatch
(
    const nonConformalDiscreteMixingPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(pp, bm),
    nonConformalCoupledPolyPatch(*this, pp),
    nbrPatchName_(pp.nbrPatchName_),
    nbrPatchID_(-1),
    rotationAxis_(pp.rotationAxis_),
    rotationCentre_(pp.rotationCentre_),
    sectorAngles_(pp.sectorAngles_),
    nSectors_(pp.nSectors_),
    patchTransforms_(pp.patchTransforms_),
    intersectionsAreValid_(false),
    intersectionsPtr_(nullptr)
{}


Foam::nonConformalDiscreteMixingPolyPatch::nonConformalDiscreteMixingPolyPatch
(
    const nonConformalDiscreteMixingPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart,
    const word& nbrPatchName,
    const word& origPatchName
)
:
    coupledPolyPatch(pp, bm, index, newSize, newStart),
    nonConformalCoupledPolyPatch(*this, origPatchName),
    nbrPatchName_(nbrPatchName),
    nbrPatchID_(-1),
    rotationAxis_(pp.rotationAxis_),
    rotationCentre_(pp.rotationCentre_),
    sectorAngles_(pp.sectorAngles_),
    nSectors_(pp.nSectors_),
    patchTransforms_(pp.patchTransforms_),
    intersectionsAreValid_(false),
    intersectionsPtr_(nullptr)
{
    if (nbrPatchName == name())
    {
        FatalErrorInFunction
            << "Neighbour patch name " << nbrPatchName
            << " cannot be the same as this patch: " << name()
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nonConformalDiscreteMixingPolyPatch::
~nonConformalDiscreteMixingPolyPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::nonConformalDiscreteMixingPolyPatch::nbrPatchID() const
{
    if (nbrPatchID_ == -1)
    {
        nbrPatchID_ = this->boundaryMesh().findPatchID(nbrPatchName());

        if (nbrPatchID_ == -1)
        {
            FatalErrorInFunction
                << "Illegal neighbour patch name " << nbrPatchName()
                << endl << "Valid patch names are "
                << this->boundaryMesh().names()
                << exit(FatalError);
        }

        // Check that it is a non-conformal discrete mixing patch
        const nonConformalDiscreteMixingPolyPatch& nbrPatch =
            refCast<const nonConformalDiscreteMixingPolyPatch>
            (
                this->boundaryMesh()[nbrPatchID_]
            );

        if (nbrPatch.nbrPatchName() != name())
        {
            WarningInFunction
                << "Patch " << name()
                << " specifies neighbour patch " << nbrPatchName()
                << endl << " but that in return specifies "
                << nbrPatch.nbrPatchName()
                << endl;
        }
    }

    return nbrPatchID_;
}


const Foam::PtrList<Foam::patchToPatches::intersection>&
Foam::nonConformalDiscreteMixingPolyPatch::intersections() const
{
    if (!owner())
    {
        FatalErrorInFunction
            << "The non-conformal discrete mixing intersections are only "
            << "available to the owner patch" << abort(FatalError);
    }

    if (!intersectionsAreValid_)
    {
        Info<< nl << "Creating intersections between non-conformal discrete "
            << "mixing patches for all rotational couples" << incrIndent << nl
            << indent << "Owner side: " << nSectors_
            << " sector" << (nSectors_ == 1 ? "" : "s") << nl
            << indent << "Neighbour side: " << nbrPatch().nSectors_
            << " sector" << (nbrPatch().nSectors_ == 1 ? "" : "s")
            << decrIndent << nl << endl;

        intersectionsPtr_.reset
        (
            new discreteMixingIntersection(*this, nbrPatch())
        );

        intersectionsAreValid_ = true;
    }

    return intersectionsPtr_->patchToPatchList();
}


void Foam::nonConformalDiscreteMixingPolyPatch::write(Ostream& os) const
{
    coupledPolyPatch::write(os);

    if (!nbrPatchName_.empty())
    {
        os.writeEntry("neighbourPatch", nbrPatchName_);
    }

    os.writeEntry("rotationAxis", rotationAxis_);
    os.writeEntry("rotationCentre", rotationCentre_);

    // Switch to ASCII for writing the Pair structure as a list
    auto oldFmt = os.format();
    os.format(IOstream::ASCII);

    os.writeKeyword("sectorAngles")
        << token::BEGIN_LIST
        << sectorAngles_.first() << token::SPACE << sectorAngles_.second()
        << token::END_LIST << token::END_STATEMENT << nl;

    os.format(oldFmt);

    if (nSectors_ > 0)
    {
        os.writeEntry("nSectors", nSectors_);
    }

    nonConformalCoupledPolyPatch::write(os);
}


// ************************************************************************* //
