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
    (c) 2017 OpenCFD Ltd.
    (c) 2011-2022 OpenFOAM Foundation
    (c) 2022-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "fieldValues/surfaceFieldValue/surfaceFieldValue.H"
#include "meshes/primitiveMesh/PatchTools/PatchTools.H"
#include "fvMesh/fvMeshStitchers/fvMeshStitcher/fvMeshStitcher.H"
#include "meshes/polyMesh/polyTopoChangeMap/polyTopoChangeMap.H"
#include "meshes/polyMesh/polyMeshMap/polyMeshMap.H"
#include "fvMesh/fvPatches/constraint/processorCyclic/processorCyclicFvPatch.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace fieldValues
{
    defineTypeNameAndDebug(surfaceFieldValue, 0);
    addToRunTimeSelectionTable(fieldValue, surfaceFieldValue, dictionary);
    addToRunTimeSelectionTable(functionObject, surfaceFieldValue, dictionary);
}
}
}

template<>
const char* Foam::NamedEnum
<
    Foam::functionObjects::fieldValues::surfaceFieldValue::regionTypes,
    4
>::names[] =
{
    "faceZone",
    "patch",
    "surface",
    "sampledSurface"
};

template<>
const char* Foam::NamedEnum
<
    Foam::functionObjects::fieldValues::surfaceFieldValue::operationType,
    17
>::names[] =
{
    "none",
    "sum",
    "weightedSum",
    "sumMag",
    "sumDirection",
    "sumDirectionBalance",
    "average",
    "weightedAverage",
    "areaAverage",
    "weightedAreaAverage",
    "areaIntegrate",
    "weightedAreaIntegrate",
    "min",
    "max",
    "CoV",
    "areaNormalAverage",
    "areaNormalIntegrate"
};

template<>
const char* Foam::NamedEnum
<
    Foam::functionObjects::fieldValues::surfaceFieldValue::postOperationType,
    2
>::names[] =
{
    "none",
    "sqrt"
};


const Foam::NamedEnum
<
    Foam::functionObjects::fieldValues::surfaceFieldValue::regionTypes,
    4
> Foam::functionObjects::fieldValues::surfaceFieldValue::regionTypeNames_;

const Foam::NamedEnum
<
    Foam::functionObjects::fieldValues::surfaceFieldValue::operationType,
    17
> Foam::functionObjects::fieldValues::surfaceFieldValue::operationTypeNames_;

const Foam::NamedEnum
<
    Foam::functionObjects::fieldValues::surfaceFieldValue::postOperationType,
    2
>
Foam::functionObjects::fieldValues::surfaceFieldValue::postOperationTypeNames_;


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

const Foam::objectRegistry&
Foam::functionObjects::fieldValues::surfaceFieldValue::obr() const
{
    if (regionType_ == stSurface)
    {
        return mesh_.lookupObject<objectRegistry>(regionName_);
    }
    else
    {
        return mesh_;
    }
}


void Foam::functionObjects::fieldValues::surfaceFieldValue::setFaceZoneFaces()
{
    const label zoneId = mesh_.faceZones().findZoneID(regionName_);

    if (zoneId < 0)
    {
        FatalErrorInFunction
            << type() << " " << name() << ": "
            << regionTypeNames_[regionType_] << "(" << regionName_ << "):" << nl
            << "    Unknown face zone name: " << regionName_
            << ". Valid face zones are: " << mesh_.faceZones().names()
            << nl << exit(FatalError);
    }

    // Ensure addressing is built on all processes
    mesh_.polyBFacePatches();
    mesh_.polyBFacePatchFaces();

    const faceZone& zone = mesh_.faceZones()[zoneId];

    DynamicList<label> faceIds(zone.size());
    DynamicList<label> facePatchIds(zone.size());
    DynamicList<bool> faceFlips(zone.size());

    forAll(zone, zoneFacei)
    {
        const label facei = zone[zoneFacei];
        const label faceFlip = zone.flipMap()[zoneFacei] ? true : false;

        if (mesh_.isInternalFace(facei))
        {
            faceIds.append(facei);
            facePatchIds.append(-1);
            faceFlips.append(faceFlip);
        }
        else
        {
            const label bFacei = facei - mesh_.nInternalFaces();

            const labelUList patches = mesh_.polyBFacePatches()[bFacei];
            const labelUList patchFaces = mesh_.polyBFacePatchFaces()[bFacei];

            forAll(patches, i)
            {
                // Don't include processor patch faces twice
                const fvPatch& fvp = mesh_.boundary()[patches[i]];
                if
                (
                    isType<processorFvPatch>(fvp)
                 && refCast<const processorFvPatch>(fvp).neighbour()
                )
                {
                    continue;
                }

                faceIds.append(patchFaces[i]);
                facePatchIds.append(patches[i]);
                faceFlips.append(faceFlip);
            }
        }
    }

    faceId_.transfer(faceIds);
    facePatchId_.transfer(facePatchIds);
    faceFlip_.transfer(faceFlips);

    nFaces_ = returnReduce(faceId_.size(), sumOp<label>());
}


void Foam::functionObjects::fieldValues::surfaceFieldValue::setPatchFaces()
{
    const label patchId = mesh_.boundaryMesh().findPatchID(regionName_);

    if (patchId < 0)
    {
        FatalErrorInFunction
            << type() << " " << name() << ": "
            << regionTypeNames_[regionType_] << "(" << regionName_ << "):" << nl
            << "    Unknown patch name: " << regionName_
            << ". Valid patch names are: "
            << mesh_.boundaryMesh().names() << nl
            << exit(FatalError);
    }

    const fvPatch& fvp = mesh_.boundary()[patchId];

    faceId_ = identity(fvp.size());
    facePatchId_ = labelList(fvp.size(), patchId);
    faceFlip_ = boolList(fvp.size(), false);

    // If we have selected a cyclic, also include any associated processor
    // cyclic faces.
    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& fvp = mesh_.boundary()[patchi];

        if
        (
            isA<processorCyclicFvPatch>(fvp)
         && refCast<const processorCyclicFvPatch>(fvp).referPatchID() == patchId
        )
        {
            faceId_.append(identity(fvp.size()));
            facePatchId_.append(labelList(fvp.size(), patchi));
            faceFlip_.append(boolList(fvp.size(), false));
        }
    }

    nFaces_ = returnReduce(faceId_.size(), sumOp<label>());
}


void Foam::functionObjects::fieldValues::surfaceFieldValue::combineMeshGeometry
(
    faceList& faces,
    pointField& points
) const
{
    List<faceList> allFaces(Pstream::nProcs());
    List<pointField> allPoints(Pstream::nProcs());

    labelList globalFacesIs(faceId_);
    forAll(globalFacesIs, i)
    {
        if (facePatchId_[i] != -1)
        {
            label patchi = facePatchId_[i];
            globalFacesIs[i] += mesh_.boundaryMesh()[patchi].start();
        }
    }

    // Add local faces and points to the all* lists
    indirectPrimitivePatch pp
    (
        IndirectList<face>(mesh_.faces(), globalFacesIs),
        mesh_.points()
    );
    allFaces[Pstream::myProcNo()] = pp.localFaces();
    allPoints[Pstream::myProcNo()] = pp.localPoints();

    Pstream::gatherList(allFaces);
    Pstream::gatherList(allPoints);

    // Renumber and flatten
    label nFaces = 0;
    label nPoints = 0;
    forAll(allFaces, proci)
    {
        nFaces += allFaces[proci].size();
        nPoints += allPoints[proci].size();
    }

    faces.setSize(nFaces);
    points.setSize(nPoints);

    nFaces = 0;
    nPoints = 0;

    // My own data first
    {
        const faceList& fcs = allFaces[Pstream::myProcNo()];
        forAll(fcs, i)
        {
            const face& f = fcs[i];
            face& newF = faces[nFaces++];
            newF.setSize(f.size());
            forAll(f, fp)
            {
                newF[fp] = f[fp] + nPoints;
            }
        }

        const pointField& pts = allPoints[Pstream::myProcNo()];
        forAll(pts, i)
        {
            points[nPoints++] = pts[i];
        }
    }

    // Other proc data follows
    forAll(allFaces, proci)
    {
        if (proci != Pstream::myProcNo())
        {
            const faceList& fcs = allFaces[proci];
            forAll(fcs, i)
            {
                const face& f = fcs[i];
                face& newF = faces[nFaces++];
                newF.setSize(f.size());
                forAll(f, fp)
                {
                    newF[fp] = f[fp] + nPoints;
                }
            }

            const pointField& pts = allPoints[proci];
            forAll(pts, i)
            {
                points[nPoints++] = pts[i];
            }
        }
    }

    // Merge
    labelList oldToNew;
    pointField newPoints;
    bool hasMerged = mergePoints
    (
        points,
        SMALL,
        false,
        oldToNew,
        newPoints
    );

    if (hasMerged)
    {
        if (debug)
        {
            Pout<< "Merged from " << points.size()
                << " down to " << newPoints.size() << " points" << endl;
        }

        points.transfer(newPoints);
        forAll(faces, i)
        {
            inplaceRenumber(oldToNew, faces[i]);
        }
    }
}


void Foam::functionObjects::fieldValues::surfaceFieldValue::
combineSurfaceGeometry
(
    faceList& faces,
    pointField& points
) const
{
    if (regionType_ == stSurface)
    {
        const surfMesh& s = dynamicCast<const surfMesh>(obr());

        if (Pstream::parRun())
        {
            // Dimension as fraction of surface
            const scalar mergeDim = 1e-10*boundBox(s.points(), true).mag();

            labelList pointsMap;

            PatchTools::gatherAndMerge
            (
                mergeDim,
                primitivePatch
                (
                    SubList<face>(s.faces(), s.faces().size()),
                    s.points()
                ),
                points,
                faces,
                pointsMap
            );
        }
        else
        {
            faces = s.faces();
            points = s.points();
        }
    }
    else if (surfacePtr_.valid())
    {
        const sampledSurface& s = surfacePtr_();

        if (Pstream::parRun())
        {
            // Dimension as fraction of mesh bounding box
            const scalar mergeDim = 1e-10*mesh_.bounds().mag();

            labelList pointsMap;

            PatchTools::gatherAndMerge
            (
                mergeDim,
                primitivePatch
                (
                    SubList<face>(s.faces(), s.faces().size()),
                    s.points()
                ),
                points,
                faces,
                pointsMap
            );
        }
        else
        {
            faces = s.faces();
            points = s.points();
        }
    }
}


Foam::scalar
Foam::functionObjects::fieldValues::surfaceFieldValue::totalArea() const
{
    scalar totalArea;

    if (regionType_ == stSurface)
    {
        const surfMesh& s = dynamicCast<const surfMesh>(obr());

        totalArea = gSum(s.magSf());
    }
    else if (surfacePtr_.valid())
    {
        totalArea = gSum(surfacePtr_().magSf());
    }
    else
    {
        totalArea = gSum(filterField(mesh_.magSf()));
    }

    return totalArea;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::fieldValues::surfaceFieldValue::needsSf() const
{
    // Many operations use the Sf field
    switch (operation_)
    {
        case opNone:
        case opSum:
        case opSumMag:
        case opAverage:
        case opMin:
        case opMax:
        {
            return false;
        }
        default:
        {
            return true;
        }
    }
}


bool Foam::functionObjects::fieldValues::surfaceFieldValue::needsWeight() const
{
    // Only a few operations use weight field
    switch (operation_)
    {
        case opWeightedSum:
        case opWeightedAverage:
        case opWeightedAreaAverage:
        case opWeightedAreaIntegrate:
        {
            return true;
        }
        default:
        {
            return false;
        }
    }
}


void Foam::functionObjects::fieldValues::surfaceFieldValue::initialise
(
    const dictionary& dict
)
{
    regionName_ = dict.lookup<word>("name");

    switch (regionType_)
    {
        case stFaceZone:
        {
            setFaceZoneFaces();
            surfacePtr_.clear();
            break;
        }
        case stPatch:
        {
            setPatchFaces();
            surfacePtr_.clear();
            break;
        }
        case stSurface:
        {
            const surfMesh& s = dynamicCast<const surfMesh>(obr());
            nFaces_ = returnReduce(s.size(), sumOp<label>());

            faceId_.clear();
            facePatchId_.clear();
            faceFlip_.clear();
            surfacePtr_.clear();
            break;
        }
        case stSampledSurface:
        {
            faceId_.clear();
            facePatchId_.clear();
            faceFlip_.clear();

            surfacePtr_ = sampledSurface::New
            (
                name(),
                mesh_,
                dict.subDict("sampledSurfaceDict")
            );
            surfacePtr_().update();
            nFaces_ =
                returnReduce(surfacePtr_().faces().size(), sumOp<label>());

            break;
        }
        default:
        {
            FatalErrorInFunction
                << type() << " " << name() << ": "
                << int(regionType_) << "(" << regionName_ << "):"
                << nl << "    Unknown region type. Valid region types are:"
                << regionTypeNames_ << nl
                << exit(FatalError);
        }
    }

    if (nFaces_ == 0 && (!mesh_.stitcher().stitches() || !mesh_.conformal()))
    {
        FatalErrorInFunction
            << type() << " " << name() << ": "
            << regionTypeNames_[regionType_] << "(" << regionName_ << "):" << nl
            << "    Region has no faces" << exit(FatalError);
    }

    if (surfacePtr_.valid())
    {
        surfacePtr_().update();
    }

    totalArea_ = totalArea();

    Info<< type() << " " << name() << ":" << nl
        << "    operation     = ";

    if (postOperation_ != postOpNone)
    {
        Info<< postOperationTypeNames_[postOperation_] << '('
            << operationTypeNames_[operation_] << ')'  << nl;
    }
    else
    {
        Info<< operationTypeNames_[operation_] << nl;
    }
    Info<< "    total faces   = " << nFaces_ << nl
        << "    total area    = " << totalArea_ << nl;

    weightFieldNames_.clear();
    if (needsWeight())
    {
        bool missing = true;
        if (dict.readIfPresent("weightFields", weightFieldNames_))
        {
            missing = false;
        }
        else
        {
            weightFieldNames_.resize(1);

            if (dict.readIfPresent("weightField", weightFieldNames_.first()))
            {
                missing = false;
                if ("none" == weightFieldNames_.first())
                {
                    // "none" == no weighting
                    weightFieldNames_.clear();
                }
            }
        }

        if (missing)
        {
            // Suggest possible alternative unweighted operation?
            FatalIOErrorInFunction(dict)
                << "The '" << operationTypeNames_[operation_]
                << "' operation is missing a weightField." << nl
                << "Either provide the weightField, "
                << "use weightField 'none' to suppress weighting," << nl
                << "or use a different operation."
                << exit(FatalIOError);
        }

        Info<< "    weight field  = ";
        if (weightFieldNames_.empty())
        {
            Info<< "none" << nl;
        }
        else
        {
            Info<< flatOutput(weightFieldNames_) << nl;
        }
    }

    // Backwards compatibility for v1612+ and older
    List<word> orientedFields;
    if (dict.readIfPresent("orientedFields", orientedFields))
    {
        WarningInFunction
            << "The 'orientedFields' option is deprecated.  These fields can "
            << "and have been added to the standard 'fields' list."
            << endl;

        fields_.append(orientedFields);
    }


    surfaceWriterPtr_.clear();
    if (writeFields_)
    {
        const word surfaceFormat(dict.lookup("surfaceFormat"));

        if (surfaceFormat != "none")
        {
            surfaceWriterPtr_.reset
            (
                surfaceWriter::New
                (
                    surfaceFormat,
                    dict.subOrEmptyDict("formatOptions").
                        subOrEmptyDict(surfaceFormat)
                ).ptr()
            );

            Info<< "    surfaceFormat = " << surfaceFormat << nl;
        }
    }

    Info<< nl << endl;
}


void Foam::functionObjects::fieldValues::surfaceFieldValue::writeFileHeader
(
    Ostream& os
) const
{
    if (operation_ != opNone)
    {
        writeCommented(os, "Region type : ");
        os << regionTypeNames_[regionType_] << " " << regionName_ << endl;

        writeHeaderValue(os, "Faces", nFaces_);
        writeHeaderValue(os, "Area", totalArea_);
        writeHeaderValue(os, "Scale factor", scaleFactor_);

        if (weightFieldNames_.size())
        {
            writeHeaderValue
            (
                os,
                "Weight field",
                flatOutput(weightFieldNames_, FlatOutput::BareComma{})
            );
        }

        writeCommented(os, "Time");
        if (writeArea_)
        {
            os  << tab << "Area";
        }

        forAll(fields_, i)
        {
            os  << tab << operationTypeNames_[operation_]
                << "(" << fields_[i] << ")";
        }

        os  << endl;
    }
}


template<>
Foam::scalar
Foam::functionObjects::fieldValues::surfaceFieldValue::processValues
(
    const Field<scalar>& values,
    const vectorField& Sf,
    const scalarField& weightField
) const
{
    switch (operation_)
    {
        case opSumDirection:
        {
            vector n(dict_.lookup("direction"));
            return gSum(pos0(values*(Sf & n))*mag(values));
        }
        case opSumDirectionBalance:
        {
            vector n(dict_.lookup("direction"));
            const scalarField nv(values*(Sf & n));

            return gSum(pos0(nv)*mag(values) - neg(nv)*mag(values));
        }
        default:
        {
            // Fall through to other operations
            return processSameTypeValues(values, Sf, weightField);
        }
    }
}


template<>
Foam::vector
Foam::functionObjects::fieldValues::surfaceFieldValue::processValues
(
    const Field<vector>& values,
    const vectorField& Sf,
    const scalarField& weightField
) const
{
    switch (operation_)
    {
        case opSumDirection:
        {
            vector n(dict_.lookup("direction"));
            n /= mag(n) + ROOTVSMALL;

            const scalarField nv(n & values);
            return gSum(pos0(nv)*n*(nv));
        }
        case opSumDirectionBalance:
        {
            vector n(dict_.lookup("direction"));
            n /= mag(n) + ROOTVSMALL;

            const scalarField nv(n & values);
            return gSum(pos0(nv)*n*(nv));
        }
        case opAreaNormalAverage:
        {
            const scalar val = gSum(values & Sf)/gSum(mag(Sf));
            return vector(val, 0, 0);
        }
        case opAreaNormalIntegrate:
        {
            const scalar val = gSum(values & Sf);
            return vector(val, 0, 0);
        }
        default:
        {
            // Fall through to other operations
            return processSameTypeValues(values, Sf, weightField);
        }
    }
}


template<>
Foam::tmp<Foam::scalarField>
Foam::functionObjects::fieldValues::surfaceFieldValue::weightingFactor
(
    const Field<scalar>& weightField
)
{
    // pass through
    return weightField;
}


template<>
Foam::tmp<Foam::scalarField>
Foam::functionObjects::fieldValues::surfaceFieldValue::weightingFactor
(
    const Field<scalar>& weightField,
    const vectorField& Sf
)
{
    // scalar * Area
    if (returnReduce(weightField.empty(), andOp<bool>()))
    {
        return mag(Sf);
    }
    else
    {
        return weightField * mag(Sf);
    }
}


template<>
Foam::tmp<Foam::scalarField>
Foam::functionObjects::fieldValues::surfaceFieldValue::weightingFactor
(
    const Field<vector>& weightField,
    const vectorField& Sf
)
{
    // vector (dot) Area
    if (returnReduce(weightField.empty(), andOp<bool>()))
    {
        return mag(Sf);
    }
    else
    {
        return weightField & Sf;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldValues::surfaceFieldValue::surfaceFieldValue
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldValue(name, runTime, dict, typeName),
    regionType_(regionTypeNames_.read(dict.lookup("regionType"))),
    operation_(operationTypeNames_.read(dict.lookup("operation"))),
    postOperation_
    (
        postOperationTypeNames_
        [dict.lookupOrDefault<word>("postOperation", "none")]
    ),
    weightFieldNames_(),
    writeArea_(dict.lookupOrDefault("writeArea", false)),
    nFaces_(0),
    faceId_(),
    facePatchId_(),
    faceFlip_()
{
    read(dict);
    writeFileHeader(file());
}


Foam::functionObjects::fieldValues::surfaceFieldValue::surfaceFieldValue
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    fieldValue(name, obr, dict, typeName),
    regionType_(regionTypeNames_.read(dict.lookup("regionType"))),
    operation_(operationTypeNames_.read(dict.lookup("operation"))),
    postOperation_
    (
        postOperationTypeNames_
        [dict.lookupOrDefault<word>("postOperation", "none")]
    ),
    weightFieldNames_(),
    writeArea_(dict.lookupOrDefault("writeArea", false)),
    nFaces_(0),
    faceId_(),
    facePatchId_(),
    faceFlip_()
{
    read(dict);
    writeFileHeader(file());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::fieldValues::surfaceFieldValue::~surfaceFieldValue()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fieldValues::surfaceFieldValue::read
(
    const dictionary& dict
)
{
    fieldValue::read(dict);
    initialise(dict);

    return true;
}


bool Foam::functionObjects::fieldValues::surfaceFieldValue::write()
{
    if (surfacePtr_.valid())
    {
        surfacePtr_().update();
    }

    if (operation_ != opNone)
    {
        fieldValue::write();

        if (Pstream::master())
        {
            writeTime(file());
        }
    }

    if (writeArea_)
    {
        totalArea_ = totalArea();
        Log << "    total area = " << totalArea_ << endl;

        if (operation_ != opNone && Pstream::master())
        {
            file() << tab << totalArea_;
        }
    }

    // Many operations use the Sf field
    vectorField Sf;
    if (needsSf())
    {
        if (regionType_ == stSurface)
        {
            const surfMesh& s = dynamicCast<const surfMesh>(obr());
            Sf = s.Sf();
        }
        else if (surfacePtr_.valid())
        {
            Sf = surfacePtr_().Sf();
        }
        else
        {
            Sf = filterField(mesh_.Sf());
        }
    }

    // Faces and points for surface format (if specified)
    faceList faces;
    pointField points;

    if (surfaceWriterPtr_.valid())
    {
        if (regionType_ == stSurface || surfacePtr_.valid())
        {
            combineSurfaceGeometry(faces, points);
        }
        else
        {
            combineMeshGeometry(faces, points);
        }
    }

    meshedSurfRef surfToWrite(points, faces);

    // Default is a zero-size scalar weight field (ie, weight = 1)
    scalarField scalarWeights;
    vectorField vectorWeights;

    for (const word& weightName : weightFieldNames_)
    {
        if (validField<scalar>(weightName))
        {
            tmp<scalarField> tfld = getFieldValues<scalar>(weightName, true);

            if (scalarWeights.empty())
            {
                scalarWeights = tfld;
            }
            else
            {
                scalarWeights *= tfld;
            }
        }
        else if (validField<vector>(weightName))
        {
            tmp<vectorField> tfld = getFieldValues<vector>(weightName, true);

            if (vectorWeights.empty())
            {
                vectorWeights = tfld;
            }
            else
            {
                FatalErrorInFunction
                    << "weightField " << weightName
                    << " - only one vector weight field allowed. " << nl
                    << "weights: " << flatOutput(weightFieldNames_) << nl
                    << abort(FatalError);
            }
        }
        else if (weightName != "none")
        {
            // Silently ignore "none", flag everything else as an error

            // TBD: treat missing "rho" like incompressible with rho = 1
            // and/or provided rhoRef value

            FatalErrorInFunction
                << "weightField " << weightName
                << " not found or an unsupported type" << nl
                << abort(FatalError);
        }
    }

    // Process the fields
    if (returnReduce(!vectorWeights.empty(), orOp<bool>()))
    {
        if (scalarWeights.size())
        {
            vectorWeights *= scalarWeights;
        }

        writeAll(Sf, vectorWeights, surfToWrite);
    }
    else
    {
        writeAll(Sf, scalarWeights, surfToWrite);
    }

    if (operation_ != opNone && Pstream::master())
    {
        file() << endl;
    }

    Log << endl;

    return true;
}


void Foam::functionObjects::fieldValues::surfaceFieldValue::movePoints
(
    const polyMesh& mesh
)
{
    bool isNonConformalPatch = false;
    if (regionType_ == stPatch)
    {
        const label patchId = mesh_.boundaryMesh().findPatchID(regionName_);
        const fvPatch& fvp = mesh_.boundary()[patchId];

        if (isA<nonConformalFvPatch>(fvp))
        {
            isNonConformalPatch = true;
        }
    }

    // It may be necessary to reset if the mesh moves
    if ((&mesh == &mesh_) && isNonConformalPatch)
    {
        // If non-conformal faces are selected, reset the function object.
        initialise(dict_);
    }
    else
    {
        // Otherwise, just recalculate patch or region area.
        totalArea_ = totalArea();
    }
}


void Foam::functionObjects::fieldValues::surfaceFieldValue::topoChange
(
    const polyTopoChangeMap& map
)
{
    if (&map.mesh() == &mesh_)
    {
        initialise(dict_);
    }
}


void Foam::functionObjects::fieldValues::surfaceFieldValue::mapMesh
(
    const polyMeshMap& map
)
{
    if (&map.mesh() == &mesh_)
    {
        initialise(dict_);
    }
}


// ************************************************************************* //
