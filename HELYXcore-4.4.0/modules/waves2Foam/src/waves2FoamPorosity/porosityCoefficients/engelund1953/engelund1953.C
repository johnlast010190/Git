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
    (c) held by original author

\*---------------------------------------------------------------------------*/

#include "engelund1953.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(engelund1953, 0);
addToRunTimeSelectionTable(porosityCoefficient, engelund1953, porosityCoefficient);

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

engelund1953::engelund1953
(
    const dictionary & poroProp
)
:
    porosityCoefficient( poroProp )
{
    dimensionedScalar d50(poroProperties_.lookup("d50"));

    dimensionedScalar alpha(poroProperties_.lookup("alpha"));

    dimensionedScalar beta(poroProperties_.lookup("beta"));

    dimensionedScalar KC
    (
    	poroProperties_.lookupOrDefault<dimensionedScalar>
        (
    	    "KC",
    	    dimensionedScalar("KC", dimless, 10000)
    	)
    );

//    scalar poro(poroProperties_.lookup<scalar>("porosity"));
    scalar poro(readResistancePorosity(poroProperties_));

    // Compute linear resistance coefficient
    dimensionedVector d( alpha * Foam::pow3( 1 - poro ) / Foam::sqr( poro ) / Foam::sqr(d50) * vector::one);

    linearCoefficient_.value() = d.value();
    linearCoefficient_.dimensions().reset( d.dimensions() );

    // Compute quadratic resistance coefficient
    dimensionedVector f( 2.0 * beta * (1 - poro) / Foam::pow3( poro ) / d50 * vector::one);

    quadraticCoefficient_.value() = f.value();
    quadraticCoefficient_.dimensions().reset( f.dimensions() );

    Info<< "Coefficients (engelund1953): " << linearCoefficient_.value() << "\t" << quadraticCoefficient_.value() << "\n" << endl;
}

engelund1953::~engelund1953()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
