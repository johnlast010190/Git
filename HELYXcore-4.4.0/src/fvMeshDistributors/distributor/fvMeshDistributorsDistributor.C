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
    (c) 2021 OpenFOAM Foundation
    (c) 2023 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "distributor/fvMeshDistributorsDistributor.H"
#include "decompositionMethod/decompositionMethod.H"
#include "fvMeshDistribute/fvMeshDistribute.H"
#include "meshes/polyMesh/polyDistributionMap/polyDistributionMap.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshDistributors
{
    defineTypeNameAndDebug(distributor, 0);
    addToRunTimeSelectionTable
    (
        fvMeshDistributor,
        distributor,
        fvMesh
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fvMeshDistributors::distributor::readDict()
{
    const dictionary& distributorDict(dict());

    redistributionInterval_ =
        distributorDict.lookupOrDefault("redistributionInterval", 10);

    maxImbalance_ =
        distributorDict.lookupOrDefault<scalar>("maxImbalance", 0.1);
}


void Foam::fvMeshDistributors::distributor::redecompose
(
    const scalarField& weights
)
{
    fvMesh& mesh = this->mesh();

    IOdictionary decompositionDict
    (
        IOobject
        (
            "decomposeParDict",
            mesh.time().system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    labelList finalDecomp;

    // Create decompositionMethod and new decomposition
    {
        autoPtr<decompositionMethod> distributor_
        (
            decompositionMethod::New
            (
                decompositionDict
            )
        );

        if (!distributor_().parallelAware())
        {
            WarningInFunction
                << "You have selected decomposition method "
                << distributor_().typeName
                << " which does" << endl
                << "not synchronise the decomposition across"
                << " processor patches." << endl
                << "    You might want to select a decomposition method which"
                << " is aware of this. Continuing."
                << endl;
        }

        finalDecomp = distributor_().decompose(mesh, weights);
    }

    // Mesh distribution engine
    fvMeshDistribute distributor(mesh, globalMeshData::matchTol_*mesh.bounds().mag());

    // Do actual sending/receiving of mesh
    autoPtr<polyDistributionMap> map
    (
        distributor.distribute(finalDecomp)
    );

    mesh.distribute(map);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshDistributors::distributor::distributor(fvMesh& mesh)
:
    fvMeshDistributor(mesh),
    redistributionInterval_(1),
    maxImbalance_(0.1),
    timeIndex_(-1)
{
    readDict();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshDistributors::distributor::~distributor()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fvMeshDistributors::distributor::update()
{
    const fvMesh& mesh = this->mesh();

    bool redistributed = false;

    if
    (
        Pstream::nProcs() > 1
     && mesh.time().timeIndex() > 1
     && timeIndex_ != mesh.time().timeIndex()
     && mesh.time().timeIndex() % redistributionInterval_ == 0
    )
    {
        timeIndex_ = mesh.time().timeIndex();

        const scalar idealNCells =
            mesh.globalData().nTotalCells()/Pstream::nProcs();

        const scalar imbalance = returnReduce
        (
            mag(1 - mesh.nCells()/idealNCells),
            maxOp<scalar>()
        );

        if (imbalance > maxImbalance_)
        {
            if (debug)
            {
                Info<< "Redistributing mesh with imbalance "
                    << imbalance << endl;
            }

            redecompose(scalarField());

            redistributed = true;
        }
    }

    return redistributed;
}


void Foam::fvMeshDistributors::distributor::topoChange(const polyTopoChangeMap&)
{}


void Foam::fvMeshDistributors::distributor::mapMesh(const polyMeshMap&)
{}


void Foam::fvMeshDistributors::distributor::distribute
(
    const polyDistributionMap&
)
{}


bool Foam::fvMeshDistributors::distributor::write(const bool write) const
{
    return true;
}


// ************************************************************************* //
