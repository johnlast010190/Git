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
    (c) 2019 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "rearLift/rearLiftObjectiveFunctionObject.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fvMesh/fvPatches/derived/wall/wallFvPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(rearLiftObjectiveFunctionObject, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        rearLiftObjectiveFunctionObject,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::rearLiftObjectiveFunctionObject::
rearLiftObjectiveFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& objectiveDict,
    bool useAdjointFileFormat
)
:
    objectiveFunctionObject(runTime, const_cast<dictionary&>(objectiveDict)),
    forceDirection_(objectiveDict.lookup("forceDirection")),
    momentDirection_(objectiveDict.lookup("momentDirection")),
    rotationCentre_(objectiveDict.lookup("rotationCentre")),
    axisLength_(objectiveDict.lookup<scalar>("axisLength"))
{
    createFiles(useAdjointFileFormat);
    forceDirection_ /= mag(forceDirection_);
    momentDirection_ /= mag(momentDirection_);
}


Foam::functionObjects::rearLiftObjectiveFunctionObject::
rearLiftObjectiveFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& objectiveDict
)
:
    rearLiftObjectiveFunctionObject(name, runTime, objectiveDict, false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::rearLiftObjectiveFunctionObject::
~rearLiftObjectiveFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool
Foam::functionObjects::rearLiftObjectiveFunctionObject::read
(
    const dictionary& dict
)
{
    objectiveFunctionObject::read(dict);

    return true;
}


bool
Foam::functionObjects::rearLiftObjectiveFunctionObject::execute()
{
    tmp<volScalarField> pp = objectiveFunctionObject::P();
    const volScalarField::Boundary& Pw
        = pp().boundaryField();

    tmp<volSymmTensorField> tau = devRhoReff();

    const volSymmTensorField::Boundary& tauw
        = tau().boundaryField();

    const surfaceVectorField::Boundary& Sfb =
        mesh_.Sf().boundaryField();

    vector totalForce  = vector::zero;
    vector totalMoment = vector::zero;

    forAll(objectivePatch_, pI)
    {
        if (objectivePatch_[pI] && isA<wallFvPatch>(mesh_.boundary()[pI]))
        {
            const vectorField& Cf = mesh_.boundary()[pI].Cf();
            //forces
            vectorField pf( Sfb[pI]*Pw[pI] );
            vectorField vf( Sfb[pI] & tauw[pI] );

            totalForce += sum(pf + vf);
            //moments
            vectorField pM( ((Cf - rotationCentre_)^Sfb[pI])*Pw[pI] );
            vectorField vM( (Cf - rotationCentre_)^(tauw[pI] & Sfb[pI]) );
            totalMoment += sum(pM + vM);
        }
    }

    reduce(totalForce , sumOp<vector>());
    reduce(totalMoment, sumOp<vector>());
    objectiveValue_ = 0.5*(totalForce & forceDirection_)
                    - (totalMoment & momentDirection_)/axisLength_;

    Info<< type() << " " << name() << " execute:" << nl
        << "Total rear lift = " << objectiveValue_ << " [N]" << nl << endl;

    return true;
}


bool
Foam::functionObjects::rearLiftObjectiveFunctionObject::write()
{
    writeOutputToFile();

    return true;
}


// ************************************************************************* //
