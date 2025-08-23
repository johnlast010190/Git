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
    (c) 2019-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "patchSwitch/switches/tableSwitch/tableSwitch.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshGIBChangers
{
    defineTypeNameAndDebug(tableSwitch, 0);
    addToRunTimeSelectionTable(GIBSwitch, tableSwitch, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::fvMeshGIBChangers::tableSwitch::enableCondition()
{
    bool tableValue = boolTable_->value(mesh_.time().value());

    bool condition = true;

    if (condition != tableValue)
    {
        condition = false;
    }

    Info<< this->type() << " for GIB faceZone "
         << mesh_.faceZones()[zoneID_].name()
         << ": open, disableCondition - boolean table: " << tableValue << endl;

    return condition;
}


bool Foam::fvMeshGIBChangers::tableSwitch::disableCondition()
{
    bool tableValue = boolTable_->value(mesh_.time().value());

    bool condition = false;

    if (condition == tableValue)
    {
        condition = true;
    }

    Info<< this->type() << " for GIB faceZone "
         << mesh_.faceZones()[zoneID_].name()
         << ": closed, enableCondition - boolean table: " << tableValue << endl;

    return condition;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshGIBChangers::tableSwitch::tableSwitch(const fvMesh& mesh, const dictionary& dict)
:
    GIBSwitch(mesh, dict),
    boolTable_(Function1<bool>::New("boolTable", dict))
{}


// ************************************************************************* //
