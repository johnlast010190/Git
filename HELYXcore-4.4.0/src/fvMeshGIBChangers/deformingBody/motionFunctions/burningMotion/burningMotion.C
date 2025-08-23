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
    (c) 2018-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "deformingBody/motionFunctions/burningMotion/burningMotion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshGIBChangers
{
    defineTypeNameAndDebug(burningMotion, 0);

    addToRunTimeSelectionTable
    (
        motionFunction,
        burningMotion,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::fvMeshGIBChangers::burningMotion::burnRate(const label& pI) const
{
    tmp<scalarField> Unt
    (
        new scalarField(mesh_.boundary()[pI].size(), 0)
    );
    scalarField& Un = Unt.ref();

    if (mode_ == "uniform")
    {
        forAll(Un, fI)
        {
            Un[fI] = uniformBurnRate_().value(mesh_.time().value());
        }
    }
    else if (mode_ == "userField")
    {
        if (mesh_.foundObject<volScalarField>(fieldName_))
        {
            const volScalarField& uf =
                mesh_.lookupObject<volScalarField>(fieldName_);
            const volScalarField::Boundary& bf = uf.boundaryField();

            Un = bf[pI].patchInternalField();
        }
    }
    else if (mode_ == "userExpression")
    {
        #ifndef HELYX_SWAK4FOAM_DISABLED
        const volScalarField userExpression(expression_()());
        const volScalarField::Boundary& bf = userExpression.boundaryField();

        Un = bf[pI].patchInternalField();
        #endif
    }
    else
    {
        FatalErrorInFunction
            << "Mode has to be specified inside dynamicFvMeshDict" << nl
            << "Supported modes are: uniform / userField / userExpression"
            << abort(FatalError);
    }

    return Unt;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshGIBChangers::burningMotion::burningMotion
(
    const fvMesh& mesh,
    const dictionary& DBMFCoeffs,
    const Time& runTime
)
:
    motionFunction(mesh, DBMFCoeffs, runTime),
    mode_("none"),
    rhos_(0),
    fieldName_("interfaceRegressionSpeed"),
    scaleFactor_(1)
{
    burningMotion::read(DBMFCoeffs);
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::fvMeshGIBChangers::burningMotion::~burningMotion() {}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::fvMeshGIBChangers::burningMotion::read
(
    const dictionary& DBMFCoeffs
)
{
    mode_ = DBMFCoeffs.lookupOrDefault<word>("mode", "userField");
    rhos_ = DBMFCoeffs.lookup<scalar>("solidDensity");

    if (mode_ == "uniform")
    {
        uniformBurnRate_.reset
        (
            Function1<scalar>::New("uniformBurnRate", DBMFCoeffs)
        );
    }
    else if (mode_ == "userField")
    {
        fieldName_ =
            DBMFCoeffs.lookupOrDefault<word>
            (
                "fieldName", "interfaceRegressionSpeed"
            );
    }
    else if (mode_ == "userExpression")
    {
        // Read expression from dict
        #ifndef HELYX_SWAK4FOAM_DISABLED
        Info<< "Reading expression" << nl << endl;

        expression_.reset
        (
            new expressionSource<scalar>
            (
                DBMFCoeffs,
                mesh_
            )
        );
        #else
        FatalErrorInFunction
            << "HELYX was compiled without swak4Foam, which is required to "
            << "support user expressions."
            << exit(FatalError);
        #endif
    }

    scaleFactor_ = DBMFCoeffs.lookupOrDefault<scalar>("scaleFactor", 1);

    return true;
}


Foam::tmp<Foam::vectorField>
Foam::fvMeshGIBChangers::burningMotion::boundaryVelocity
(
    const label& pMaster
) const
{
    tmp<vectorField> tVelocityBoundary
    (
        new vectorField(mesh_.boundary()[pMaster].size(), vector::zero)
    );
    vectorField& velocityBoundary = tVelocityBoundary.ref();

    tmp<vectorField> tInterSpeed = interfaceVelocity(pMaster);
    const vectorField& interSpeed = tInterSpeed();

    const volScalarField& rho = mesh_.lookupObject<volScalarField>("rho");

    velocityBoundary = (1-rhos_/rho.boundaryField()[pMaster])*interSpeed;

    velocityBoundary /= scaleFactor_;

    return tVelocityBoundary;
}


Foam::tmp<Foam::vectorField>
Foam::fvMeshGIBChangers::burningMotion::interfaceVelocity
(
    const label& pMaster
) const
{
    tmp<scalarField> Unt = burnRate(pMaster);

    tmp<vectorField> tSpeedInter(new vectorField(Unt->size(), vector::zero));

    vectorField& speedInter = tSpeedInter.ref();

    speedInter = Unt()*mesh_.boundary()[pMaster].nf()();

    speedInter *= scaleFactor_;

    return tSpeedInter;
}


// ************************************************************************* //
