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
    (c) YEAR Engys Ltd

\*---------------------------------------------------------------------------*/

#include "FUNCTIONOBJECT.H"
#include "Time.H"
#include "fvMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(FUNCTIONOBJECT, 0);
    addToRunTimeSelectionTable(functionObject, FUNCTIONOBJECT, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::FUNCTIONOBJECT::FUNCTIONOBJECT
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    wordData_(dict.lookupOrDefault<word>("wordData", "defaultWord")),
    scalarData_(dict.lookup<scalar>("scalarData")),
    labelData_(dict.lookup<label>("labelData"))
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::FUNCTIONOBJECT::~FUNCTIONOBJECT()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::FUNCTIONOBJECT::read(const dictionary& dict)
{
    dict.readIfPresent("wordData", wordData_);
    scalarData_ = dict.lookup<scalar>("scalarData");
    labelData_ = dict.lookup<label>("labelData");

    return true;
}


bool Foam::functionObjects::FUNCTIONOBJECT::execute()
{
    return true;
}


bool Foam::functionObjects::FUNCTIONOBJECT::end()
{
    return true;
}


bool Foam::functionObjects::FUNCTIONOBJECT::write()
{
    return true;
}


// ************************************************************************* //
