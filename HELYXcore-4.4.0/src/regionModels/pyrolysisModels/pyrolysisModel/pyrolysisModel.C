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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "pyrolysisModels/pyrolysisModel/pyrolysisModel.H"
#include "pyrolysisModels/thermo/pyrolysisThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
    defineTypeNameAndDebug(pyrolysisModel, 0);
    defineRunTimeSelectionTable(pyrolysisModel, mesh);
    defineRunTimeSelectionTable(pyrolysisModel, dictionary);
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::regionModels::pyrolysisModel::initialise(const word& modelType)
{
    if (basicThermo::dictName == basicThermo::matDictName)
    {
        thermoPtr_ =
            &refCast<multiphaseThermo>
            (
                basicThermo::lookupOrCreate(regionMesh().thisDb())
            );
        solidThermo_ = &const_cast<solidThermo&>(thermoPtr_->sThermo());
    }
    else
    {
        if (modelType == pyrolysisThermo::typeName)
        {
            solidThermo_ = solidThermo::New(regionMesh()).ptr();
            (*solidThermo_).db().store
            (
                dynamic_cast<basicThermo::implementation*>(solidThermo_)
            );
        }
        else
        {
            solidThermo_ = solidMulticomponentThermo::New(regionMesh()).ptr();
            (*solidThermo_).db().store
            (
                dynamic_cast<basicThermo::implementation*>(solidThermo_)
            );
        }
    }

    if (modelType != pyrolysisThermo::typeName)
    {
        solidChemistry_.reset
        (
            basicSolidChemistryModel::New
            (
                dynamic_cast<solidMulticomponentThermo&>(*solidThermo_)
            ).ptr()
        );
    }

    radiation_.reset
    (
        radiationModel::New
        (
            pyrolysisModel::sThermo().T(),
            pyrolysisModel::sThermo().TRef()
        ).ptr()
    );

    nMagSfPtr_.reset
    (
        new surfaceScalarField
        (
            IOobject("nMagSf", time().timeName(), regionMesh()),
            regionMesh(),
            dimensionedScalar(dimArea, 0)
        )
    );

    // Calculate boundaryFaceFaces and boundaryFaceCells
    DynamicList<label> faceIDs;
    DynamicList<label> cellIDs;

    label localPyrolysisFacei = 0;

    const polyBoundaryMesh& rbm = regionMesh().boundaryMesh();

    forAll(intCoupledPatchIDs_, i)
    {
        const label patchi = intCoupledPatchIDs_[i];
        const polyPatch& ppCoupled = rbm[patchi];
        localPyrolysisFacei += ppCoupled.size();
    }

    boundaryFaceOppositeFace_.setSize(localPyrolysisFacei);
    boundaryFaceFaces_.setSize(localPyrolysisFacei);
    boundaryFaceCells_.setSize(localPyrolysisFacei);

    localPyrolysisFacei = 0;

    forAll(intCoupledPatchIDs_, i)
    {
        const label patchi = intCoupledPatchIDs_[i];
        const polyPatch& ppCoupled = rbm[patchi];
        forAll(ppCoupled, localFacei)
        {
            label facei = ppCoupled.start() + localFacei;
            label celli = -1;
            label nCells = 0;
            do
            {
                label ownCelli = regionMesh().faceOwner()[facei];
                if (ownCelli != celli)
                {
                    celli = ownCelli;
                }
                else
                {
                    celli = regionMesh().faceNeighbour()[facei];
                }
                nCells++;
                cellIDs.append(celli);
                const cell& cFaces = regionMesh().cells()[celli];
                faceIDs.append(facei);
                label face0 =
                    cFaces.opposingFaceLabel(facei, regionMesh().faces());
                facei = face0;
            } while (regionMesh().isInternalFace(facei));

            boundaryFaceOppositeFace_[localPyrolysisFacei] = facei;
            boundaryFaceFaces_[localPyrolysisFacei].transfer(faceIDs);
            boundaryFaceCells_[localPyrolysisFacei].transfer(cellIDs);

            localPyrolysisFacei++;
            nLayers_ = nCells;
        }
    }
    faceIDs.clear();
    cellIDs.clear();

    surfaceScalarField& nMagSf = nMagSfPtr_();

    surfaceScalarField::Boundary& nMagSfBf = nMagSf.boundaryFieldRef();

    localPyrolysisFacei = 0;

    forAll(intCoupledPatchIDs_, i)
    {
        const label patchi = intCoupledPatchIDs_[i];
        const polyPatch& ppCoupled = rbm[patchi];
        const vectorField& pNormals = ppCoupled.faceNormals();

        nMagSfBf[patchi] = regionMesh().Sf().boundaryField()[patchi] & pNormals;

        forAll(pNormals, localFacei)
        {
            const vector n = pNormals[localFacei];
            const labelList& faces = boundaryFaceFaces_[localPyrolysisFacei++];
            forAll(faces, facei)
            {
                // facei = 0 is on boundary
                if (facei > 0)
                {
                    const label faceID = faces[facei];
                    nMagSf[faceID] = regionMesh().Sf()[faceID] & n;
                }
            }
        }
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::labelField> Foam::regionModels::pyrolysisModel::moveMesh
(
    const scalarList& deltaV,
    const scalar minDelta
)
{
    tmp<labelField> tcellMoveMap(new labelField(regionMesh().nCells(), 0));
    labelField& cellMoveMap = tcellMoveMap.ref();

    if (!moveMesh_)
    {
        return cellMoveMap;
    }

    pointField oldPoints = regionMesh().points();
    pointField newPoints = oldPoints;

    const polyBoundaryMesh& bm = regionMesh().boundaryMesh();

    label totalFacei = 0;
    forAll(intCoupledPatchIDs_, localPatchi)
    {
        label patchi = intCoupledPatchIDs_[localPatchi];
        const polyPatch& pp = bm[patchi];

        forAll(pp, patchFacei)
        {
            const labelList& faces = boundaryFaceFaces_[totalFacei];
            const labelList& cells = boundaryFaceCells_[totalFacei];
            const label oFace = boundaryFaceOppositeFace_[totalFacei];

            const vector n = pp.faceNormals()[patchFacei];
            const vector sf = pp.faceAreas()[patchFacei];

            List<point> oldCf(faces.size() + 1, vector::zero);
            List<bool> frozen(faces.size(), false);

            forAll(faces, i)
            {
                oldCf[i] = regionMesh().faceCentres()[faces[i]];
            }

            oldCf[faces.size()] = regionMesh().faceCentres()[oFace];

            forAll(faces, i)
            {
                const label celli = cells[i];

                if (mag(oldCf[i + 1] - oldCf[i]) < minDelta)
                {
                    frozen[i] = true;
                    cellMoveMap[celli] = 1;
                }
            }

            vectorField newDelta(cells.size() + 1, vector::zero);

            label j = 0;
            forAllReverse (cells, i)
            {
                const label celli = cells[i];
                newDelta[j+1] = (deltaV[celli]/mag(sf))*n + newDelta[j];
                j++;
            }

            forAll(faces, i)
            {
                const label facei = faces[i];
                const face f = regionMesh().faces()[facei];

                if (!frozen[i])
                {
                    forAll(f, pti)
                    {
                        const label pointi = f[pti];

                        newPoints[pointi] =
                            oldPoints[pointi]
                          + newDelta[newDelta.size() - 1 - i];
                    }
                }
            }

            totalFacei++;
        }
    }

    // Move points
    regionMesh().movePoints(newPoints);

    return tcellMoveMap;
}


bool Foam::regionModels::pyrolysisModel::read()
{
    moveMesh_.readIfPresent("moveMesh", coeffs_);
    if (regionModel::read())
    {
        maxDiff_ = time().controlDict().lookup<scalar>("maxDi");
        moveMesh_.readIfPresent("moveMesh", coeffs_);
        return true;
    }
    else
    {
        return false;
    }
}


bool Foam::regionModels::pyrolysisModel::read(const dictionary& dict)
{
    moveMesh_.readIfPresent("moveMesh", coeffs_);
    if (regionModel::read(dict))
    {
        maxDiff_ = time().controlDict().lookup<scalar>("maxDi");
        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionModels::pyrolysisModel::pyrolysisModel
(
    const word& modelType,
    const fvMesh& mesh,
    const word& regionType
)
:
    regionModel(mesh, regionType, modelType, false),
    pimple_(regionMesh()),
    thermoPtr_(nullptr),
    nLayers_(0),
    nMagSfPtr_(nullptr),
    moveMesh_(true),
    maxDiff_(10)
{
    if (active())
    {
        pyrolysisModel::read();
        initialise(modelType);
    }
}


Foam::regionModels::pyrolysisModel::pyrolysisModel
(
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& regionType
)
:
    regionModel(mesh, regionType, modelType, dict, true),
    pimple_(regionMesh()),
    thermoPtr_(nullptr),
    nLayers_(0),
    nMagSfPtr_(nullptr),
    moveMesh_(false),
    maxDiff_(10)
{
    if (active())
    {
        pyrolysisModel::read(dict);
        initialise(modelType);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionModels::pyrolysisModel::~pyrolysisModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::regionModels::pyrolysisModel::addMassSources
(
    const label patchi,
    const label facei
)
{
    return 0.0;
}


Foam::scalar Foam::regionModels::pyrolysisModel::solidRegionDiffNo() const
{
    return -GREAT;
}


Foam::scalar Foam::regionModels::pyrolysisModel::maxDiff() const
{
    return maxDiff_;
}


// ************************************************************************* //
