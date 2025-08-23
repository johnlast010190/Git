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
    (c) 2011 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "tetherPotential/derived/pitchForkRing/pitchForkRing.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace tetherPotentials
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pitchForkRing, 0);

addToRunTimeSelectionTable
(
    tetherPotential,
    pitchForkRing,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pitchForkRing::pitchForkRing
(
    const word& name,
    const dictionary& tetherPotentialProperties
)
:
    tetherPotential(name, tetherPotentialProperties),
    pitchForkRingCoeffs_
    (
        tetherPotentialProperties.subDict(typeName + "Coeffs")
    ),
    mu_(pitchForkRingCoeffs_.lookup<scalar>("mu")),
    alpha_(pitchForkRingCoeffs_.lookup<scalar>("alpha")),
    rOrbit_(pitchForkRingCoeffs_.lookup<scalar>("rOrbit"))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

scalar pitchForkRing::energy(const vector r) const
{
    scalar p = sqrt(r.x()*r.x() + r.y()*r.y());

    scalar pMinusRSqr = sqr(p - rOrbit_);

    return
       -0.5 * mu_ * pMinusRSqr
      + 0.25 * pMinusRSqr * pMinusRSqr
      + 0.5 * alpha_ * r.z() * r.z();
}


vector pitchForkRing::force(const vector r) const
{
    scalar p = sqrt(r.x()*r.x() + r.y()*r.y());

    scalar pMinusR = (p - rOrbit_);

    return vector
    (
        (mu_ - sqr(pMinusR)) * pMinusR * r.x()/(p + VSMALL),
        (mu_ - sqr(pMinusR)) * pMinusR * r.y()/(p + VSMALL),
      - alpha_ * r.z()
    );
}


bool pitchForkRing::read(const dictionary& tetherPotentialProperties)
{
    tetherPotential::read(tetherPotentialProperties);

    pitchForkRingCoeffs_ =
        tetherPotentialProperties.subDict(typeName + "Coeffs");

    mu_ = pitchForkRingCoeffs_.lookup<scalar>("mu");
    alpha_ = pitchForkRingCoeffs_.lookup<scalar>("alpha");
    rOrbit_ = pitchForkRingCoeffs_.lookup<scalar>("rOrbit");

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace tetherPotentials
} // End namespace Foam

// ************************************************************************* //
