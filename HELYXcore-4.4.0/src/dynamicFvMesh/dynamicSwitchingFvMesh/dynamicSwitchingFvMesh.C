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
    (c) 2019 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "dynamicSwitchingFvMesh/dynamicSwitchingFvMesh.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/surfaceFields/surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dynamicSwitchingFvMesh, 0);
    addToRunTimeSelectionTable
    (
        dynamicFvMesh,
        dynamicSwitchingFvMesh,
        IOobject
    );
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::wordList Foam::dynamicSwitchingFvMesh::readRegistries
(
    const IOobject& io
)
{
    wordList reg(0);

    const dictionary meshDict
    (
        dynamicMeshDict().optionalSubDict(typeName + "Coeffs")
    );

    forAllConstIter(dictionary, meshDict, iter)
    {
        if (iter().isDict())
        {
            reg.append(iter().dict().dictName());
        }
    }

    return reg;
}


Foam::boolListList Foam::dynamicSwitchingFvMesh::readPatches
(
    const IOobject& io,
    const wordList& registries
)
{
    boolListList regPatchList
    (
        registries.size(),
        List<bool>(boundaryMesh().size(), true)
    );

    const dictionary meshDict
    (
        dynamicMeshDict().optionalSubDict(typeName + "Coeffs")
    );

    forAll(registries, rI)
    {
        dictionary reg = meshDict.subDict(registries[rI]);

        if (reg.found("offPatches"))
        {
            wordList offPatches = reg.lookup("offPatches");

            forAll(boundaryMesh(), pI)
            {
                forAll(offPatches, opI)
                {
                    if (offPatches[opI] == boundaryMesh()[pI].name())
                    {
                        regPatchList[rI][pI] = false;
                    }
                }
            }
        }
    }

    return regPatchList;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicSwitchingFvMesh::dynamicSwitchingFvMesh(const IOobject& io)
:
    dynamicFvMesh(io),
    currentRegistryName_(this->name()),
    registries_(readRegistries(io)),
    regPatchState_(readPatches(io, registries_))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dynamicSwitchingFvMesh::~dynamicSwitchingFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dynamicSwitchingFvMesh::update()
{
    return true;
}


bool Foam::dynamicSwitchingFvMesh::update(const word& registryName)
{
    if (currentRegistryName_ == registryName)
    {
        return true;
    }

    label regI = -1;
    forAll(registries_, rI)
    {
        if (registries_[rI] == registryName)
        {
            regI = rI;
            break;
        }
    }

    if (regI == -1)
    {
        FatalErrorInFunction
            << "Request for mesh configuration for registry "
            << registryName << " failed. "
            << "List of available registries is: " << registries_
            << exit(FatalError);
    }

    // Calculate centres and areas; currently there is no function for only
    // calculating the face areas (for either the whole mesh or boundaryMesh)
    vectorField fCentres(nFaces());
    vectorField faceAreas(nFaces());
    scalarField magfaceAreas(nFaces());
    makeFaceCentresAndAreas(points(), fCentres, faceAreas, magfaceAreas);

    surfaceVectorField& Sf = const_cast<surfaceVectorField&>(this->Sf());
    surfaceScalarField& magSf = const_cast<surfaceScalarField&>(this->magSf());

    forAll(Sf.boundaryField(), pI)
    {
        label startI = boundaryMesh()[pI].start();

        vectorField& Sfp = Sf.boundaryFieldRef()[pI];
        scalarField& magSfp = magSf.boundaryFieldRef()[pI];

        forAll(Sfp, fI)
        {
            Sfp[fI] = faceAreas[startI + fI];
            magSfp[fI] = mag(Sfp[fI]);
        }
    }

    const fvBoundaryMesh& boundary = this->boundary();
    forAll(boundary, pI)
    {
        if (!regPatchState_[regI][pI])
        {
            vectorField& Sfp = Sf.boundaryFieldRef()[pI];
            scalarField& magSfp = magSf.boundaryFieldRef()[pI];

            forAll(Sfp, fI)
            {
                Sfp[fI] = vector::zero;
                magSfp[fI] *= ROOTVSMALL;
            }
        }
    }

    currentRegistryName_ = registryName;

    return true;
}


// ************************************************************************* //
