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
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "acousticSource.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "global/constants/mathematical/mathematicalConstants.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(acousticSource, 0);
        addToRunTimeSelectionTable
        (
            option,
            acousticSource,
            dictionary
        );
    }

    template<>
    const char* Foam::NamedEnum
    <
        Foam::fv::acousticSource::volumeModeType,
        2
    >::names[] =
    {
        "absolute",
        "specific"
    };

    template<>
    const char* Foam::NamedEnum
    <
        Foam::fv::acousticSource::acousticSourceModeType,
        3
    >::names[] =
    {
        "sinusoidal",
        "sinusoidalSuperposition",
        "table"
    };
}

const Foam::NamedEnum<Foam::fv::acousticSource::volumeModeType, 2>
Foam::fv::acousticSource::volumeModeTypeNames_;

const Foam::NamedEnum<Foam::fv::acousticSource::acousticSourceModeType, 3>
Foam::fv::acousticSource::acousticSourceModeTypeNames_;

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void Foam::fv::acousticSource::updateSource(const fvMatrix<scalar>& eqn)
{
    const scalar t = obr_.time().value();

    if (t <= duration_)
    {
        if (acousticSourceMode_ == sinusoidal)
        {
            A_ = AF1_->value(t)*sin(mathematical::twoPi*frequency_*t);
        }
        else if (acousticSourceMode_ == sinusoidalSuperposition)
        {
            A_ = 0.0;
            scalar df = (fMax_ - fMin_)/nComps_;

            for (label i=0; i < nComps_; i++)
            {
                A_ += AF1_->value(t)*sin(mathematical::twoPi*(fMin_ + i*df)*t);
            }
        }
        else if (acousticSourceMode_ == table)
        {
            A_ = AF1_->value(t);
        }
        else
        {
            FatalErrorInFunction
                << "Unknown acousticSourceMode type "
                << acousticSourceModeTypeNames_[acousticSourceMode_]
                << ". Valid acousticSourceMode types are:" << nl
                << acousticSourceModeTypeNames_
                << exit(FatalError);
        }
    }
    else
    {
        A_ = 0.0;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::acousticSource::acousticSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    cellSetOption(name, modelType, dict, obr),
    volumeMode_(vmAbsolute),
    VDash_(1.0),
    acousticSourceMode_(sinusoidal),
    AF1_(nullptr),
    A_(0.0),
    amplitude_(0.0),
    frequency_(0.0),
    duration_(GREAT),
    fMin_(0.0),
    fMax_(0.0),
    nComps_(1)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::acousticSource::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldI
)
{
    updateSource(eqn);

    volScalarField::Internal Su
    (
        IOobject
        (
            "Su",
            mesh_.time().timeName(),
            obr_
        ),
        mesh_,
        dimensionedScalar(eqn.dimensions()/dimVolume, 0)
    );

    UIndirectList<scalar>(Su, cells_) = A_/VDash_;

    eqn += Su;
}


void Foam::fv::acousticSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldI
)
{
    addSup(eqn, fieldI);
}


void Foam::fv::acousticSource::sourceFields
(
    wordList& fieldNames
)
{
    fieldNames.setSize(1);
    fieldNames[0] = word(coeffs_.lookup("fieldName"));
}


bool Foam::fv::acousticSource::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        volumeMode_ = volumeModeTypeNames_.read(coeffs_.lookup("volumeMode"));
        acousticSourceMode_
            = acousticSourceModeTypeNames_.read(coeffs_.lookup("acousticSourceMode"));

        // Set volume normalisation
        if (volumeMode_ == vmAbsolute)
        {
            VDash_ = V_;
        }

        AF1_ = Function1<scalar>::New("amplitude", coeffs_);
        duration_ = coeffs_.lookupOrDefault<scalar>("duration", GREAT);

        if (acousticSourceMode_ == sinusoidal)
        {
            frequency_ = coeffs_.lookup<scalar>("frequency");
        }
        else if (acousticSourceMode_ == sinusoidalSuperposition)
        {
            fMin_ = coeffs_.lookupOrDefault<scalar>("fMin", 1.0);
            fMax_ = coeffs_.lookupOrDefault<scalar>("fMax", 1.0);
            nComps_ = coeffs_.lookupOrDefault<scalar>("nComps", 1);
        }

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
