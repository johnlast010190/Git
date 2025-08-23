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
    (c) 2016 OpenCFD Ltd.
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "acousticDampingSource.H"
#include "fvMesh/fvMesh.H"
#include "fvMatrices/fvMatrices.H"
#include "finiteVolume/fvm/fvmSup.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "global/constants/mathematical/mathematicalConstants.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(acousticDampingSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        acousticDampingSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fv::acousticDampingSource::setBlendingFactor()
{
    blendFactor_.primitiveFieldRef() = 1;

    const vectorField& Cf = mesh_.C();

    const scalar pi = constant::mathematical::pi;

    forAll(cells_, i)
    {
        label celli = cells_[i];
        scalar d = mag(Cf[celli] - x0_);

        if (d < r1_)
        {
            blendFactor_[celli] = 0.0;
        }
        else if ((d >= r1_) && (d <= r2_))
        {
            blendFactor_[celli] = (1.0 - cos(pi*mag(d - r1_)/(r2_ - r1_)))/2.0;
        }
    }

    blendFactor_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::acousticDampingSource::acousticDampingSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    cellSetOption(name, modelType, dict, obr),
    frequency_("frequency", dimless/dimTime, 0),
    blendFactor_
    (
        volScalarField
        (
            IOobject
            (
                name_ + ":blend",
                mesh_.time().timeName(),
                obr,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimless, 1),
            zeroGradientFvPatchField<vector>::typeName
        )
    ),
    URefName_("unknown-URef"),
    x0_(Zero),
    r1_(0),
    r2_(0),
    w_(20)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::acousticDampingSource::sourceFields(wordList& fieldNames)
{
    if (coeffs_.found("UNames"))
    {
        fieldNames = coeffs_.lookup<wordList>("UNames");
    }
    else if (coeffs_.found("UName"))
    {
        fieldNames = wordList(1, coeffs_.lookup<word>("UName"));
    }
    else
    {
        fieldNames = wordList(1, "U");
    }
}


void Foam::fv::acousticDampingSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    const volVectorField& U = eqn.psi();
    const volScalarField coeff(name_ + ":coeff", w_*frequency_*blendFactor_);
    const volVectorField& URef(obr_.lookupObject<volVectorField>(URefName_));

    fvMatrix<vector> dampingEqn
    (
        fvm::Sp(coeff, U) - coeff*URef
    );
    eqn -= dampingEqn;
}


void Foam::fv::acousticDampingSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    const volVectorField& U = eqn.psi();
    const volScalarField coeff(name_ + ":coeff", w_*frequency_*blendFactor_);
    const volVectorField& URef(obr_.lookupObject<volVectorField>(URefName_));

    fvMatrix<vector> dampingEqn
    (
        fvm::Sp(rho*coeff, U) - rho*coeff*URef
    );
    eqn -= dampingEqn;
}


void Foam::fv::acousticDampingSource::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    const volVectorField& U = eqn.psi();
    const volScalarField coeff(name_ + ":coeff", w_*frequency_*blendFactor_);
    const volVectorField& URef(obr_.lookupObject<volVectorField>(URefName_));

    fvMatrix<vector> dampingEqn
    (
        fvm::Sp(alpha*rho*coeff, U) - alpha*rho*coeff*URef
    );
    eqn -= dampingEqn;
}


bool Foam::fv::acousticDampingSource::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        frequency_.value() = coeffs_.lookup<scalar>("frequency");
        URefName_ = coeffs_.lookup<word>("URef");
        x0_ = coeffs_.lookup<point>("centre");
        r1_ = coeffs_.lookup<scalar>("radius1");
        r2_ = coeffs_.lookup<scalar>("radius2");

        if (coeffs_.readIfPresent("w", w_))
        {
            Info<< name_ << ": Setting stencil width to " << w_ << endl;
        }

        setBlendingFactor();

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
