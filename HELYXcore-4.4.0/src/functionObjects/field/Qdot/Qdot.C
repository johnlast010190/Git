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
    (c) 2019 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "Qdot/Qdot.H"
#include "combustionModel/combustionModel.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(Qdot, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        Qdot,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::Qdot::Qdot
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    phaseName_(word::null)
{
    Qdot::read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::Qdot::~Qdot()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::Qdot::read
(
    const dictionary& dict
)
{
    fvMeshFunctionObject::read(dict);

    phaseName_ = dict.lookupOrDefault<word>("phase", word::null);

    return true;
}


bool Foam::functionObjects::Qdot::execute()
{
    return true;
}


bool Foam::functionObjects::Qdot::write()
{
    const word modelName
    (
        IOobject::groupName
        (
            combustionModel::combustionPropertiesName,
            phaseName_
        )
    );
    return mesh_.lookupObject<combustionModel>(modelName).Qdot()().write();
}


// ************************************************************************* //