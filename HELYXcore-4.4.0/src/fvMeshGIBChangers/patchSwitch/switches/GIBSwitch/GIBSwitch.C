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

#include "patchSwitch/switches/GIBSwitch/GIBSwitch.H"
#include "db/IOstreams/Fstreams/IFstream.H"
#include "meshes/polyMesh/polyPatches/indirectPolyPatch/indirectPolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshGIBChangers
{
    defineTypeNameAndDebug(GIBSwitch, 0);
    defineRunTimeSelectionTable(GIBSwitch, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fvMeshGIBChangers::GIBSwitch::initialise()
{
    word zoneName(dict_.lookup("gibFaceZone"));
    zoneID_ = mesh_.faceZones().findZoneID(zoneName);
    if (zoneID_ == -1)
    {
        FatalIOErrorInFunction
        (
            dict_
        )   << "Cannot find cellZone named " << zoneName
            << ". Valid zones are " << mesh_.faceZones().names()
            << exit(FatalIOError);
    }

    findPatchIDs();

    if (!filePath().empty())
    {
        IFstream is(filePath());
        dictionary storedDataDict(is);
        bool isActive = true;
        storedDataDict.readIfPresent
        (
            "active",
            isActive
        );
        if (!isActive)
        {
            disableBoundary();
        }
    }
    else
    {
        // Probably in the 1st outer loop
        if (!startActive_)
        {
            disableBoundary();
            curCondInterval_ = 0;
        }
    }
}


void Foam::fvMeshGIBChangers::GIBSwitch::findPatchIDs()
{
    forAll(mesh_.boundary(), pI)
    {
        const polyPatch& poly = mesh_.boundary()[pI].patch();
        if (isA<indirectPolyPatch>(poly))
        {
            const indirectPolyPatch& gibPolyPatch =
                refCast<const indirectPolyPatch>(poly);
            label zoneId = gibPolyPatch.zoneId();
            const word& inPolyType = gibPolyPatch.indirectPolyPatchType();

            if (zoneID_ == zoneId)
            {
                if (inPolyType == "master")
                {
                    masterID_ = pI;
                }
                else if (inPolyType == "slave")
                {
                    slaveID_ = pI;
                }
            }
        }
    }
}


Foam::IOobject Foam::fvMeshGIBChangers::GIBSwitch::createIOobject
(
    const fvMesh& mesh,
    const dictionary& dict
) const
{
    word zoneName(dict.lookup("gibFaceZone"));

    label zoneID = mesh.faceZones().findZoneID(zoneName);
    if (zoneID == -1)
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "Cannot find cellZone named " << zoneName
            << ". Valid zones are " << mesh.faceZones().names()
            << exit(FatalIOError);
    }

    IOobject io
    (
        zoneName+"_Coeffs",
        mesh.time().timeName(),
        "uniform",
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    if (io.typeHeaderOk<IOdictionary>(true))
    {
        io.readOpt() = IOobject::MUST_READ_IF_MODIFIED;
        return io;
    }
    else
    {
        io.readOpt() = IOobject::NO_READ;
        return io;
    }
}


void Foam::fvMeshGIBChangers::GIBSwitch::enableBoundary()
{
     fvPatch& mP = const_cast<fvPatch&>(mesh_.boundary()[masterID_]);
     bool& mPact = mP.active();
     mPact = true;

     fvPatch& sP = const_cast<fvPatch&>(mesh_.boundary()[slaveID_]);
     bool& sPact = sP.active();
     sPact = true;
}


void Foam::fvMeshGIBChangers::GIBSwitch::disableBoundary()
{
     fvPatch& mP = const_cast<fvPatch&>(mesh_.boundary()[masterID_]);
     bool& mPact = mP.active();
     mPact = false;

     fvPatch& sP = const_cast<fvPatch&>(mesh_.boundary()[slaveID_]);
     bool& sPact = sP.active();
     sPact = false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshGIBChangers::GIBSwitch::GIBSwitch(const fvMesh& mesh, const dictionary& dict)
:
    IOdictionary(createIOobject(mesh, dict)),
    mesh_(mesh),
    dict_(dict),
    zoneID_(-1),
    masterID_(-1),
    slaveID_(-1),
    startActive_(dict.lookupOrDefault<Switch>("startActive", true)),
    startTime_(dict.lookupOrDefault<scalar>("startTime", -1)),
    condInterval_(dict.lookupOrDefault<label>("conditionInterval", 1)),
    curCondInterval_(1)
{
    if (!filePath().empty())
    {
        IFstream is(filePath());
        dictionary storedWeightsDict(is);
        storedWeightsDict.readIfPresent
        (
            "curConditionInterval",
            curCondInterval_
        );
    }

    initialise();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fvMeshGIBChangers::GIBSwitch::update()
{
    bool switchB = false;

    if (mesh_.time().value() >= startTime_)
    {
        bool currectStatus = mesh_.boundary()[masterID_].isActive();
        if (!currectStatus)
        {
            if (enableCondition())
            {
                if (curCondInterval_ >= condInterval_)
                {
                    enableBoundary();
                    curCondInterval_ = 1;
                    switchB = true;

                    Info<< "GIBboundaries for faceZone "
                         << mesh_.faceZones()[zoneID_].name()
                         << " enabled" << endl;
                }
                else
                {
                    curCondInterval_++;
                }
            }
        }
        else
        {
            if (disableCondition())
            {
                if (curCondInterval_ >= condInterval_)
                {
                    disableBoundary();
                    curCondInterval_ = 1;
                    switchB = true;

                    Info<< "GIBboundaries for faceZone "
                         << mesh_.faceZones()[zoneID_].name()
                         << " disabled" << endl;
                }
                else
                {
                    curCondInterval_++;
                }
            }
        }
    }
    else
    {
        Info<< "GIBboundaries for faceZone "
             << mesh_.faceZones()[zoneID_].name()
             << " starts updating at " << startTime_ << endl;
    }

    this->updateDict();

    return switchB;
}


void Foam::fvMeshGIBChangers::GIBSwitch::updateDict()
{
    this->add("active", mesh_.boundary()[masterID_].isActive(), true);
    this->add("curConditionInterval", curCondInterval_, true);
}


// ************************************************************************* //
