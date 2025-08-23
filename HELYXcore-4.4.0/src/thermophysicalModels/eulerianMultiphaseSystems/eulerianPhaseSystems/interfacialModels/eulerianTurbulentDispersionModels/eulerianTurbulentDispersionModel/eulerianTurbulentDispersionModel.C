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

#include "eulerianTurbulentDispersionModel.H"
#include "../../../eulerianPhasePair/eulerianPhasePair/eulerianPhasePair.H"
#include "finiteVolume/fvc/fvcGrad.H"
#include "interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolate.H"
#include "finiteVolume/fvc/fvcSnGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(eulerianTurbulentDispersionModel, 0);
    defineRunTimeSelectionTable(eulerianTurbulentDispersionModel, dictionary);
}

const Foam::dimensionSet Foam::eulerianTurbulentDispersionModel::dimD(1, -1, -2, 0, 0);
const Foam::dimensionSet Foam::eulerianTurbulentDispersionModel::dimF(1, -2, -2, 0, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eulerianTurbulentDispersionModel::eulerianTurbulentDispersionModel
(
    const dictionary& dict,
    const eulerianPhasePair& pair
)
:
    pair_(pair)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::eulerianTurbulentDispersionModel::~eulerianTurbulentDispersionModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField>
Foam::eulerianTurbulentDispersionModel::F() const
{
    return D()*fvc::grad(pair_.dispersed().volFrac());
}


Foam::tmp<Foam::surfaceScalarField>
Foam::eulerianTurbulentDispersionModel::Ff() const
{
    return fvc::interpolate(D())*fvc::snGrad(pair_.dispersed().volFrac());
}


// ************************************************************************* //
