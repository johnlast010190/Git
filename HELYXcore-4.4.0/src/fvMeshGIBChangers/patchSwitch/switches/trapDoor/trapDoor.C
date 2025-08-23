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

#include "patchSwitch/switches/trapDoor/trapDoor.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fields/volFields/volFields.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshGIBChangers
{
    defineTypeNameAndDebug(trapDoor, 0);
    addToRunTimeSelectionTable(GIBSwitch, trapDoor, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::fvMeshGIBChangers::trapDoor::enableCondition()
{
    // Be careful here! The GIBs are disabled.
    // If you compute values on the GIB, they will depend on the BC.

    // Check if in first timestep the U field was found in memory.
    // This can happen because the user might want to use a field that is
    // constructed and stored after the mesh (like using a FO).
    if (mesh_.time().restartTimeIndex() == 1)
    {
        if (!mesh_.foundObject<volVectorField>(UName_))
        {
            Info<< "trapDoor "
                 << mesh_.faceZones()[zoneID_].name()
                 << ": velocity field " << UName_
                 << " was not found in memory. "
                 << "Leaving the GIB disabled for this timestep." << endl;

            return false;
        }
    }

    const volVectorField& U = mesh_.lookupObject<volVectorField>(UName_);

    const surfaceScalarField& magSf = mesh_.magSf();
    const surfaceVectorField& Sf = mesh_.Sf();

    const scalarField& magSfm = magSf.boundaryField()[masterID_];
    const vectorField& Sfm = Sf.boundaryField()[masterID_];

    const labelList& fcs = mesh_.boundary()[masterID_].faceCells();
    vectorField Umc(U.boundaryField()[masterID_].size(), vector::zero);

    forAll(Umc, cI)
    {
        Umc[cI] = U[fcs[cI]];
    }

    scalar meanUm = 0.0;

    if (mesh_.foundObject<volScalarField>(fractionName_))
    {
        // Assume only the velocities of the fuel
        const volScalarField& alpha =
            mesh_.lookupObject<volScalarField>(fractionName_);
        const scalarField& alpham = alpha.boundaryField()[masterID_];

        meanUm = gSum((Umc*alpham) & Sfm)/gSum(magSfm);
    }
    else
    {
        // Assume single phase (alpha does not exist or user did not specify
        // the phase name correctly).
        meanUm = gSum(Umc & Sfm)/gSum(magSfm);
    }

    if (invertDirection_)
    {
        meanUm = -meanUm;
    }

    Info<< this->type() << " for GIB faceZone "
         << mesh_.faceZones()[zoneID_].name()
         << ": open, disableCondition - mean normal fluid velocity: "
         << meanUm << endl;

    return (meanUm >- fVel_);
}


bool Foam::fvMeshGIBChangers::trapDoor::disableCondition()
{
    // Check if in first timestep the U field was found in memory.
    // This can happen because the user might want to use a field that is
    // constructed and stored after the mesh (like using a FO).
    if (mesh_.time().restartTimeIndex() == 1)
    {
        if (!mesh_.foundObject<volScalarField>(pName_))
        {
            Info<< "trapDoor "
                 << mesh_.faceZones()[zoneID_].name()
                 << ": pressure field " << pName_
                 << " was not found in memory. "
                 << "Leaving the GIB enabled for this timestep." << endl;
            return false;
        }
    }

    const volScalarField& p = mesh_.lookupObject<volScalarField>(pName_);
    const surfaceScalarField& magSf = mesh_.magSf();

    const scalarField& magSfm = magSf.boundaryField()[masterID_];
    const scalarField& magSfs = magSf.boundaryField()[slaveID_];

    const scalarField& pm = p.boundaryField()[masterID_];
    const scalarField& ps = p.boundaryField()[slaveID_];

    scalar meanPm = gSum(pm*magSfm)/gSum(magSfm);

    scalar meanPs = gSum(ps*magSfs)/gSum(magSfs);

    scalar meanDP = meanPs - meanPm;

    if (invertDirection_)
    {
        meanDP = -meanDP;
    }

    Info<< this->type() << " for GIB faceZone "
         << mesh_.faceZones()[zoneID_].name()
         << ": closed, enableCondition - pressureDrop: "
         << meanDP << endl;

    return (meanDP > dp_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshGIBChangers::trapDoor::trapDoor(const fvMesh& mesh, const dictionary& dict)
:
    GIBSwitch(mesh, dict),
    invertDirection_(dict.lookupOrDefault<Switch>("invertDirection", false)),
    fVel_(dict.lookup<scalar>("closeDoorFuelVelocity")),
    dp_(dict.lookup<scalar>("openDoorPressureDrop")),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    pName_(dict.lookupOrDefault<word>("p", "p")),
    fractionName_(dict.lookupOrDefault<word>("alpha", "alpha"))
{}


// ************************************************************************* //
