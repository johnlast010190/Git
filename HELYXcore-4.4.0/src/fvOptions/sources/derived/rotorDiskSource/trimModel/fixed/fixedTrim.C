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
    (c) 2012-2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "fixedTrim.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "global/unitConversion/unitConversion.H"
#include "global/constants/mathematical/mathematicalConstants.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fixedTrim, 0);

    addToRunTimeSelectionTable(trimModel, fixedTrim, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedTrim::fixedTrim
(
    const fv::rotorDiskSource& rotor,
    const dictionary& dict
)
:
    trimModel(rotor, dict, typeName),
    thetag_(rotor.cells().size(), 0.0)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::fixedTrim::~fixedTrim()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fixedTrim::read(const dictionary& dict)
{
    trimModel::read(dict);

    scalar theta0 = degToRad(coeffs_.lookup<scalar>("theta0"));
    scalar theta1c = degToRad(coeffs_.lookup<scalar>("theta1c"));
    scalar theta1s = degToRad(coeffs_.lookup<scalar>("theta1s"));

    const List<point>& x = rotor_.x();
    forAll(thetag_, i)
    {
        scalar psi = x[i].y();
        thetag_[i] = theta0 + theta1c*cos(psi) + theta1s*sin(psi);
    }
}


Foam::tmp<Foam::scalarField> Foam::fixedTrim::thetag() const
{
    return tmp<scalarField>(thetag_);
}


void Foam::fixedTrim::correct
(
    const vectorField& U,
    vectorField& force
)
{}


void Foam::fixedTrim::correct
(
    const volScalarField rho,
    const vectorField& U,
    vectorField& force)
{}


// ************************************************************************* //
