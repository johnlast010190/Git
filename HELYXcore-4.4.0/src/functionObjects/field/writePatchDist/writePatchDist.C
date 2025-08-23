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
    (c) 2021 Engys Ltd.
    (c) 2016 OpenCFD Ltd.
    (c) 2016-2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "writePatchDist/writePatchDist.H"
#include "fvMesh/wallDist/wallDist/wallDist.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(writePatchDist, 0);
    addToRunTimeSelectionTable(functionObject, writePatchDist, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::writePatchDist::writePatchDist
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    mesh_(refCast<const fvMesh>(obr_)),
    patchSet_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::writePatchDist::~writePatchDist()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::writePatchDist::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    patchSet_ =
        mesh_.boundaryMesh().patchSet
        (
            wordReList(dict.lookupOrDefault("patches", wordReList()))
        );

    Info<< type() << " " << name() << ":" << nl;


    if (patchSet_.empty())
    {
        WarningInFunction
            << "Requested patch distance but no patches specified "<< endl;
    }

    return true;
}


bool Foam::functionObjects::writePatchDist::execute()
{
    return true;
}


bool Foam::functionObjects::writePatchDist::write()
{
    wallDist d(mesh_, "meshWave", patchSet_, "patch");

    volScalarField y(d.y());

    y.write();

    return true;
}


// ************************************************************************* //
