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
    (c) 2014-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "eulerianWallLubricationModel.H"
#include "../../../eulerianPhasePair/eulerianPhasePair/eulerianPhasePair.H"
#include "finiteVolume/fvc/fvcFlux.H"
#include "interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolate.H"
#include "fvMesh/fvPatches/derived/wall/wallFvPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(eulerianWallLubricationModel, 0);
    defineRunTimeSelectionTable(eulerianWallLubricationModel, dictionary);
}

const Foam::dimensionSet Foam::eulerianWallLubricationModel::dimF(1, -2, -2, 0, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField> Foam::eulerianWallLubricationModel::zeroGradWalls
(
    tmp<volVectorField> tFi
) const
{
    volVectorField& Fi = tFi.ref();
    const fvPatchList& patches =  Fi.mesh().boundary();

    volVectorField::Boundary& FiBf = Fi.boundaryFieldRef();

    forAll(patches, patchi)
    {
        if (isA<wallFvPatch>(patches[patchi]))
        {
            fvPatchVectorField& Fiw = FiBf[patchi];
            Fiw = Fiw.patchInternalField();
        }
    }

    return tFi;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eulerianWallLubricationModel::eulerianWallLubricationModel
(
    const dictionary& dict,
    const eulerianPhasePair& pair
)
:
    eulerianWallDependentModel(pair.phase1().mesh()),
    pair_(pair)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::eulerianWallLubricationModel::~eulerianWallLubricationModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField> Foam::eulerianWallLubricationModel::F() const
{
    return pair_.dispersed().volFrac()*Fi();
}


Foam::tmp<Foam::surfaceScalarField> Foam::eulerianWallLubricationModel::Ff() const
{
    return fvc::interpolate(pair_.dispersed().volFrac())*fvc::flux(Fi());
}


// ************************************************************************* //
