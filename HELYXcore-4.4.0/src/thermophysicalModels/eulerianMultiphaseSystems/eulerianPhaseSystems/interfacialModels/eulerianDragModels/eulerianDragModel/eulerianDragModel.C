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
    (c) 2011-2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "eulerianDragModel.H"
#include "../../../eulerianPhasePair/eulerianPhasePair/eulerianPhasePair.H"
#include "../../eulerianSwarmCorrections/eulerianSwarmCorrection/eulerianSwarmCorrection.H"
#include "interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolate.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(eulerianDragModel, 0);
    defineRunTimeSelectionTable(eulerianDragModel, dictionary);
}

const Foam::dimensionSet Foam::eulerianDragModel::dimK(1, -3, -1, 0, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eulerianDragModel::eulerianDragModel
(
    const eulerianPhasePair& pair,
    const bool registerObject
)
:
    regIOobject
    (
        IOobject
        (
            IOobject::groupName(typeName, pair.name()),
            pair.phase1().mesh().time().timeName(),
            pair.phase1().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            registerObject
        )
    ),
    pair_(pair)
{}


Foam::eulerianDragModel::eulerianDragModel
(
    const dictionary& dict,
    const eulerianPhasePair& pair,
    const bool registerObject
)
:
    regIOobject
    (
        IOobject
        (
            IOobject::groupName(typeName, pair.name()),
            pair.phase1().mesh().time().timeName(),
            pair.phase1().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            registerObject
        )
    ),
    pair_(pair),
    swarmCorrection_
    (
        eulerianSwarmCorrection::New(eulerianSwarmCorrection::typeName, dict, pair)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::eulerianDragModel::~eulerianDragModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::eulerianDragModel::Ki() const
{
    return
        0.75
       *CdRe()
       *swarmCorrection_->Cs()
       *pair_.continuous().rho()
       *pair_.continuous().nu()
       /sqr(pair_.dispersed().d());
}


Foam::tmp<Foam::volScalarField> Foam::eulerianDragModel::K() const
{
    return
        max
        (
            pair_.dispersed().volFrac(),
            pair_.dispersed().residualAlpha()
        )*Ki();
}


Foam::tmp<Foam::surfaceScalarField> Foam::eulerianDragModel::Kf() const
{
    return
        max
        (
            fvc::interpolate(pair_.dispersed().volFrac()),
            pair_.dispersed().residualAlpha()
        )*fvc::interpolate(Ki());
}


bool Foam::eulerianDragModel::writeData(Ostream& os) const
{
    return os.good();
}


// ************************************************************************* //
