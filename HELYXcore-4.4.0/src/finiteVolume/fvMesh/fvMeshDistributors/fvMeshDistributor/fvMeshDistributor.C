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
    (c) 2021-2023 OpenFOAM Foundation
    (c) 2023-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "fvMesh/fvMeshDistributors/fvMeshDistributor/fvMeshDistributor.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fvMeshDistributor, 0);
    defineRunTimeSelectionTable(fvMeshDistributor, fvMesh);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshDistributor::fvMeshDistributor(fvMesh& mesh)
:
    mesh_(mesh),
    dynamicMeshDict_
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                mesh.time().constant(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                false
            )
        )
    )
{}


Foam::fvMeshDistributor::velocityMotionCorrection::velocityMotionCorrection
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    velocityFields_(dict.lookupOrDefault("velocityFields", wordList()))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshDistributor::~fvMeshDistributor()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fvMeshDistributor::velocityMotionCorrection::update() const
{
    forAll(velocityFields_, i)
    {
        if (mesh_.foundObject<volVectorField>(velocityFields_[i]))
        {
            mesh_.lookupObjectRef<volVectorField>
            (
                velocityFields_[i]
            ).correctBoundaryConditions();
        }
    }
}


// ************************************************************************* //
