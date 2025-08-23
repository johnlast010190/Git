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

#include "randomPhase.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(randomPhase, 0);
addToRunTimeSelectionTable(phases, randomPhase, phases);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


randomPhase::randomPhase
(
    const Time& rT,
    dictionary& dict
)
:
    phases(rT, dict)
{
    Info<< "\nConstructing: " << this->type() << endl;
    srand(time(nullptr));

    if (dict.found("seedForRandomPhase"))
    {
        label seed = dict.lookup<label>("seedForRandomPhase");
        srand(seed);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


scalar randomPhase::phase(const scalar& freq, const vector& k)
{
    return
    (
        2.0*M_PI*static_cast<scalar>(rand())
       /static_cast<scalar>(RAND_MAX)
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
