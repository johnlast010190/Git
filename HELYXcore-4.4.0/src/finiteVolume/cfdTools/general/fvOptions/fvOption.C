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
    (c) 2011-2019 OpenFOAM Foundation
    (c) 2010-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/fvOptions/fvOption.H"
#include "fields/volFields/volFields.H"
#include "fvMatrices/fvMatrix/fvMatrix.H"
#include "fvSolutionRegistry/fvSolutionRegistry.H"
#include "db/dynamicLibrary/dlLibraryTable/dlLibraryTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
    template<>
    const char* Foam::NamedEnum< fv::option::execHookType, 6 >::names[] =
    {
        "constrain",
        "correct",
        "operator",
        "outerCorrect",
        "solve",
        "none"
    };

    namespace fv
    {
        defineTypeNameAndDebug(option, 0);
        defineRunTimeSelectionTable(option, dictionary);

        const NamedEnum<fv::option::execHookType, 6>
            fv::option::execHookTypeNames_;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::option::option
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    name_(name),
    modelType_(modelType),
    mesh_(fvSolutionRegistry::getMesh(obr)),
    obr_(obr),
    dict_(dict),
    coeffs_(dict.optionalSubDict(modelType + "Coeffs")),
    active_(dict_.lookupOrDefault<Switch>("active", true)),
    MRF_(dict_.lookupOrDefault<Switch>("MRF", false)),
    execHook_(execHookTypeNames_[dict_.lookupOrDefault<word>("hookOp", "none")]),
    fieldNames_(),
    applied_(),
    boundarySourcePatchIDs_(),
    boundaryApplied_()
{
    Info<< incrIndent << indent << "Source: " << name_ << endl << decrIndent;
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fv::option> Foam::fv::option::New
(
    const word& name,
    const dictionary& coeffs0,
    const objectRegistry& obr
)
{
    // Backward compatibility
    dictionary coeffs = coeffs0;
    if (word(coeffs.lookup("type")) == "solverOption")
    {
        WarningInFunction
            << "Using deprecated fvOption type 'solverOption'. "
            << "Please use solver type directly as an fvOption." << nl << endl;
        coeffs.merge(coeffs0.subDict("solverOptionCoeffs"));
    }

    word modelType(coeffs.lookup("type"));

    Info<< indent
        << "Selecting finite volume options model type " << modelType << endl;

    libs.open(coeffs, "libs");

    const auto ctor =
        ctorTableLookup
        (
            "model type", dictionaryConstructorTable_(), modelType
        );
    return autoPtr<option>(ctor(name, modelType, coeffs, obr));
}


Foam::fv::option::~option()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::option::isActive()
{
    return active_;
}

bool Foam::fv::option::isMRF() const
{
    return MRF_;
}

Foam::label Foam::fv::option::applyToField
(
    const word& fieldName, const word& regionName
) const
{
    forAll(fieldNames_, i)
    {
        if (fieldName == fieldNames_[i] && regionName == regionNames_[i])
        {
            return i;
        }
    }
    return -1;
}

Foam::label Foam::fv::option::applyToBoundaryFieldAndPatch
(
    const word& fieldName,
    const label patchID
) const
{
    if (!boundarySourcePatchIDs_.found(fieldName))
    {
        return -1;
    }
    else
    {
        return findIndex(boundarySourcePatchIDs_[fieldName], patchID);
    }
}


void Foam::fv::option::checkApplied() const
{
    forAll(applied_, i)
    {
        if (!applied_[i])
        {
            WarningInFunction
                << "Source " << name_ << " defined for field "
                << fieldNames_[i] << " but never used" << endl;
        }
    }
}


void Foam::fv::option::checkBoundaryApplied() const
{
    forAllConstIters(boundaryApplied_, iter)
    {
        forAll(iter(), i)
        {
            if (!iter()[i])
            {
                label patchi = boundarySourcePatchIDs_[iter.key()][i];
                WarningInFunction
                    << "Boundary source " << name_ << " defined for field "
                    << iter.key() << ", " << "patch "
                    << mesh_.boundary()[patchi].name()
                    << " but never used" << endl;
            }
        }
    }
}


void Foam::fv::option::setSourceNames()
{
    fieldNames_.clear();
    regionNames_.clear();
    sourceFields(fieldNames_, regionNames_);
    applied_ = boolList(fieldNames_.size(), false);

    boundarySourcePatchIDs_.clear();
    boundarySourceFieldsAndPatches(boundarySourcePatchIDs_);
    boundaryApplied_.clear();
    forAllConstIters(boundarySourcePatchIDs_, iter)
    {
        boundaryApplied_.set
        (
            iter.key(), boolList(iter().size(), false)
        );
    }
}


void Foam::fv::option::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{}


void Foam::fv::option::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{}


void Foam::fv::option::addSup
(
    fvMatrix<sphericalTensor>& eqn,
    const label fieldi
)
{}


void Foam::fv::option::addSup
(
    fvMatrix<symmTensor>& eqn,
    const label fieldi
)
{}


void Foam::fv::option::addSup
(
    fvMatrix<tensor>& eqn,
    const label fieldi
)
{}


void Foam::fv::option::addSup
(
    fvBlockMatrix<vector>& eqn,
    const label fieldi
)
{}


void Foam::fv::option::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{}


void Foam::fv::option::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{}


void Foam::fv::option::addSup
(
    const volScalarField& rho,
    fvMatrix<sphericalTensor>& eqn,
    const label fieldi
)
{}


void Foam::fv::option::addSup
(
    const volScalarField& rho,
    fvMatrix<symmTensor>& eqn,
    const label fieldi
)
{}


void Foam::fv::option::addSup
(
    const volScalarField& rho,
    fvMatrix<tensor>& eqn,
    const label fieldi
)
{}


void Foam::fv::option::addSup
(
    const volScalarField& rho,
    fvBlockMatrix<vector>& eqn,
    const label fieldi
)
{}


void Foam::fv::option::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    addSup(alpha*rho, eqn, fieldi);
}


void Foam::fv::option::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    addSup(alpha*rho, eqn, fieldi);
}


void Foam::fv::option::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<sphericalTensor>& eqn,
    const label fieldi
)
{
    addSup(alpha*rho, eqn, fieldi);
}


void Foam::fv::option::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<symmTensor>& eqn,
    const label fieldi
)
{
    addSup(alpha*rho, eqn, fieldi);
}


void Foam::fv::option::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<tensor>& eqn,
    const label fieldi
)
{
    addSup(alpha*rho, eqn, fieldi);
}


void Foam::fv::option::constrain(fvMatrix<scalar>& eqn, const label fieldi)
{}


void Foam::fv::option::constrain(fvMatrix<vector>& eqn, const label fieldi)
{}


void Foam::fv::option::constrain
(
    fvMatrix<sphericalTensor>& eqn,
    const label fieldi
)
{}


void Foam::fv::option::constrain
(
    fvMatrix<symmTensor>& eqn,
    const label fieldi
)
{}


void Foam::fv::option::constrain(fvMatrix<tensor>& eqn, const label fieldi)
{}


void Foam::fv::option::correct(volScalarField& field)
{}


void Foam::fv::option::correct(volVectorField& field)
{}


void Foam::fv::option::correct(volSphericalTensorField& field)
{}


void Foam::fv::option::correct(volSymmTensorField& field)
{}


void Foam::fv::option::correct(volTensorField& field)
{}


void Foam::fv::option::correct()
{}


// ************************************************************************* //
