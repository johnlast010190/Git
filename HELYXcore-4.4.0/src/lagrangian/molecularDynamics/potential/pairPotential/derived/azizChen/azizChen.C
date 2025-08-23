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

#include "pairPotential/derived/azizChen/azizChen.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace pairPotentials
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(azizChen, 0);

addToRunTimeSelectionTable
(
    pairPotential,
    azizChen,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

azizChen::azizChen
(
    const word& name,
    const dictionary& azizChen
)
:
    pairPotential(name, azizChen),
    azizChenCoeffs_(azizChen.subDict(typeName + "Coeffs")),
    epsilon_(azizChenCoeffs_.lookup<scalar>("epsilon")),
    rm_(azizChenCoeffs_.lookup<scalar>("rm")),
    A_(azizChenCoeffs_.lookup<scalar>("A")),
    alpha_(azizChenCoeffs_.lookup<scalar>("alpha")),
    C6_(azizChenCoeffs_.lookup<scalar>("C6")),
    C8_(azizChenCoeffs_.lookup<scalar>("C8")),
    C10_(azizChenCoeffs_.lookup<scalar>("C10")),
    D_(azizChenCoeffs_.lookup<scalar>("D")),
    gamma_(azizChenCoeffs_.lookup<scalar>("gamma"))
{
    setLookupTables();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

scalar azizChen::unscaledEnergy(const scalar r) const
{
    scalar x = r/rm_;

    scalar F = 1.0;

    if (x < D_)
    {
        F = exp(-pow(((D_ / x) - 1.0),2));
    }

    return
        epsilon_
       *(
            A_ * Foam::pow(x, gamma_)*exp(-alpha_*x)
          - (
                (C6_/ Foam::pow(x, 6))
              + (C8_/ Foam::pow(x, 8))
              + (C10_/ Foam::pow(x, 10))
            )
           *F
    );
}


bool azizChen::read(const dictionary& azizChen)
{
    pairPotential::read(azizChen);

    azizChenCoeffs_ = azizChen.subDict(typeName + "Coeffs");

    epsilon_ = azizChenCoeffs_.lookup<scalar>("epsilon");
    rm_ = azizChenCoeffs_.lookup<scalar>("rm");
    A_ = azizChenCoeffs_.lookup<scalar>("A");
    alpha_ = azizChenCoeffs_.lookup<scalar>("alpha");
    C6_ = azizChenCoeffs_.lookup<scalar>("C6");
    C8_ = azizChenCoeffs_.lookup<scalar>("C8");
    C10_ = azizChenCoeffs_.lookup<scalar>("C10");
    D_ = azizChenCoeffs_.lookup<scalar>("D");
    gamma_ = azizChenCoeffs_.lookup<scalar>("gamma");

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pairPotentials
} // End namespace Foam

// ************************************************************************* //
