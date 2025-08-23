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
    (c) 2014-2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "eulerianVirtualMassModel.H"
#include "../../../eulerianPhasePair/eulerianPhasePair/eulerianPhasePair.H"
#include "interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolate.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(eulerianVirtualMassModel, 0);
    defineRunTimeSelectionTable(eulerianVirtualMassModel, dictionary);
}

const Foam::dimensionSet Foam::eulerianVirtualMassModel::dimK(dimDensity);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eulerianVirtualMassModel::eulerianVirtualMassModel
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
    pair_(pair)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::eulerianVirtualMassModel::~eulerianVirtualMassModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::eulerianVirtualMassModel::Ki() const
{
    return Cvm()*pair_.continuous().rho();
}


Foam::tmp<Foam::volScalarField> Foam::eulerianVirtualMassModel::K() const
{
    return pair_.dispersed().volFrac()*Ki();
}


Foam::tmp<Foam::surfaceScalarField> Foam::eulerianVirtualMassModel::Kf() const
{
    return
        fvc::interpolate(pair_.dispersed().volFrac())*fvc::interpolate(Ki());
}


bool Foam::eulerianVirtualMassModel::writeData(Ostream& os) const
{
    return os.good();
}


// ************************************************************************* //
