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

#include "eulerianLiftModel.H"
#include "../../../eulerianPhasePair/eulerianPhasePair/eulerianPhasePair.H"
#include "finiteVolume/fvc/fvcCurl.H"
#include "finiteVolume/fvc/fvcFlux.H"
#include "interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolate.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(eulerianLiftModel, 0);
    defineRunTimeSelectionTable(eulerianLiftModel, dictionary);
}

const Foam::dimensionSet Foam::eulerianLiftModel::dimF(1, -2, -2, 0, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eulerianLiftModel::eulerianLiftModel
(
    const dictionary& dict,
    const eulerianPhasePair& pair
)
:
    pair_(pair)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::eulerianLiftModel::~eulerianLiftModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField> Foam::eulerianLiftModel::Fi() const
{
    return
        Cl()
       *pair_.continuous().rho()
       *(
            pair_.Ur() ^ fvc::curl(pair_.continuous().U())
        );
}


Foam::tmp<Foam::volVectorField> Foam::eulerianLiftModel::F() const
{
    return pair_.dispersed().volFrac()*Fi();
}


Foam::tmp<Foam::surfaceScalarField> Foam::eulerianLiftModel::Ff() const
{
    return fvc::interpolate(pair_.dispersed().volFrac())*fvc::flux(Fi());
}


// ************************************************************************* //
