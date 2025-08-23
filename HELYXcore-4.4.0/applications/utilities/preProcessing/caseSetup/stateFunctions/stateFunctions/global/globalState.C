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
    (c) 2016, 2020 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "globalState.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{

namespace stateFunctions
{

defineTypeNameAndDebug(globalState, 0);
addToRunTimeSelectionTable(stateFunction, globalState, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

globalState::globalState
(
    const Time& time,
    const dictionary& input,
    const dictionary& defaults,
    const bool distributed,
    const bool collated,
    const stateIndex& index
)
:
    stateFunction(input, defaults, distributed, collated, index),
    time_(&time)
{}

globalState::globalState
(
    word region,
    const dictionary& input,
    const dictionary& defaults,
    const stateFunction& master,
    const stateIndex& index,
    word meshName
)
:
    stateFunction(region, input, defaults, master, index, meshName),
    time_(nullptr)
{
    FatalErrorInFunction << "Invalid contructor" << endl;
}

// * * * * * * * * * * * * * * * * Functions     * * * * * * * * * * * * * * //

const Time& globalState::runTime() const
{
    if (time_ == nullptr)
    {
        FatalErrorInFunction << "Pointer to runTime is nullptr."
            << endl;
    }

    return *time_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace stateFunction
} // End namespace Foam

// ************************************************************************* //
