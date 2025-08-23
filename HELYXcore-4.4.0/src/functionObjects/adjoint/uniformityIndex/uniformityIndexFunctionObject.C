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

#include "uniformityIndex/uniformityIndexFunctionObject.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(uniformityIndexFunctionObject, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        uniformityIndexFunctionObject,
        dictionary
    );
}
}


// * * * * * * * * * * * * * *Private Member Functions * * * * * * * * * * * //

scalar
Foam::functionObjects::uniformityIndexFunctionObject::
objectivePatchArea() const
{
    scalar outletArea = 0;

    forAll(mesh_.boundary(), patchI)
    {
        const fvPatch& p = mesh_.boundary()[patchI];

        if
        (
            (
                p.type() == "outlet" ||
                p.patch().physicalType() == "outlet" ||
                p.type() == "cyclicAMI"
            ) &&
            (
                objectivePatch_[patchI]
            )
        )
        {
            outletArea += sum(p.magSf());
        }
    }

    reduce(outletArea, sumOp<scalar>());

    return outletArea;
}


scalar
Foam::functionObjects::uniformityIndexFunctionObject::
targetVelocity() const
{
    scalar Vd = 0;

    forAll(mesh_.boundary(), patchI)
    {
        const fvPatch& p = mesh_.boundary()[patchI];

        if
        (
            (
                p.type() == "outlet" ||
                p.patch().physicalType() == "outlet" ||
                p.type() == "cyclicAMI"
            ) &&
            (
                objectivePatch_[patchI]
            )
        )
        {
            Vd += sum(phi().boundaryField()[patchI]);
        }
    }

    reduce(Vd, sumOp<scalar>());

    Vd /= objPatchA_;

    //to shore things up if the initial outlet velocity is negative or zero
    Vd = sign(Vd)*(max(mag(Vd), SMALL));

    return Vd;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::uniformityIndexFunctionObject::
uniformityIndexFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& objectiveDict,
    bool useAdjointFileFormat
)
:
    objectiveFunctionObject(runTime, const_cast<dictionary&>(objectiveDict)),
    objPatchA_(objectivePatchArea())
{
    createFiles(useAdjointFileFormat);
}


Foam::functionObjects::uniformityIndexFunctionObject::
uniformityIndexFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& objectiveDict
)
:
    uniformityIndexFunctionObject(name, runTime, objectiveDict, false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::uniformityIndexFunctionObject::
~uniformityIndexFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool
Foam::functionObjects::uniformityIndexFunctionObject::read
(
    const dictionary& dict
)
{
    objectiveFunctionObject::read(dict);

    return true;
}


bool
Foam::functionObjects::uniformityIndexFunctionObject::execute()
{
    scalar Vd = targetVelocity();
    objectiveValue_ = 0;

    forAll(mesh_.boundary(), patchI)
    {
        const fvPatch& p = mesh_.boundary()[patchI];

        if
        (
            (
                p.type() == "outlet" ||
                p.patch().physicalType() == "outlet" ||
                p.type() == "cyclicAMI"
            ) &&
            (
                objectivePatch_[patchI]
            )
        )
        {
            scalarField Udev
            (
                mag(phi().boundaryField()[patchI]
                /p.magSf()/Vd - 1.0)
            );

            objectiveValue_ += sum
            (
                p.magSf()*Udev
            );
        }
    }

    reduce(objectiveValue_, sumOp<scalar>());
    objectiveValue_ = 1 - 0.5*objectiveValue_/objPatchA_;

    Info<< type() << " " << name() << " execute:" << nl
        << "Outlet uniformity index = " << objectiveValue_
        << nl << endl;

    return true;
}


bool
Foam::functionObjects::uniformityIndexFunctionObject::write()
{
    writeOutputToFile();

    return true;
}


// ************************************************************************* //
