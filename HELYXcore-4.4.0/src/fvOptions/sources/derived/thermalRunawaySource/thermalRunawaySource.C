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
    (c) 2025 Engys Ltd

\*---------------------------------------------------------------------------*/
#include "thermalRunawaySource.H"
#include "fvMesh/fvMesh.H"
#include "finiteVolume/fvm/fvmSup.H"
#include "fvMatrices/fvMatrix/fvMatrix.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "global/unitConversion/unitConversion.H"
#include "fvMatrices/fvMatrix/fvMatrix.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {

        defineTypeNameAndDebug(thermalRunawaySource, 0);
        addToRunTimeSelectionTable
        (
            option,
            thermalRunawaySource,
            dictionary
        );
    }

    template<>
    const char* Foam::NamedEnum
    <
        Foam::fv::thermalRunawaySource::volumeModeType,
        2
    >::names[] =
    {
        "absolute",
        "specific"
    };
}
const Foam::NamedEnum<Foam::fv::thermalRunawaySource::volumeModeType, 2>
Foam::fv::thermalRunawaySource::volumeModeTypeNames_;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::thermalRunawaySource::thermalRunawaySource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    cellSetOption(name, modelType, dict, obr),
    stateFunctionObject(name, obr.time()),
    volumeMode_
    (
        volumeModeTypeNames_
        [
            coeffs_.lookupOrDefault<word>("volumeMode", "absolute")
        ]
    ),
    VDash_(1.0),
    state_(coeffs_.lookupOrDefault<bool>("state", false)),
    TR_time_(0),
    Tcritical_(coeffs_.lookup<scalar>("Tcritical")),
    HRR_(nullptr)
{
    read(dict);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::thermalRunawaySource::~thermalRunawaySource()
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::fv::thermalRunawaySource::sourceFields
(
    wordList& fieldNames
)
{
    fieldNames = coeffs_.lookup<wordList>("fieldNames");
}


void Foam::fv::thermalRunawaySource::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    const Time& time(obr_.time());
    const scalarField& V(mesh_.V());
    const volScalarField& T(obr_.lookupObject<volScalarField>("T"));

    scalarField& heSource = eqn.source();

    scalar Tmean(0);
    scalar totalVolume(0);

    forAll(cells_, i)
    {
        const label celli = cells_[i];
        Tmean += T[celli]*V[celli];
        totalVolume += V[celli];
    }
    label nCells = cells_.size();
    reduce(Tmean, sumOp<scalar>());
    reduce(totalVolume, sumOp<scalar>());
    reduce(nCells, sumOp<label>());

    if (nCells>0)
    {
        Tmean /= totalVolume;
        if (debug)
            Info<< "Tmean is " << Tmean
                << ", nCells " << nCells << endl;

        if (Tmean>Tcritical_)
        {
            if (state_ == false)
            {
                TR_time_ = time.value();
                state_ = true;
                Info<< "Starting thermal runaway for "
                    << option::name_ << endl;
            }

            scalar tEff = time.value()-TR_time_;

            heSource -= HRR_->value(tEff)*V/VDash_;
            if (debug)
                Info<< "Heat release from table "
                    << HRR_->value(tEff) << endl;
        }
    }

    if (state_)
        setResult("state", state_);
    if (TR_time_>0)
        setResult("TR_time", TR_time_);
}


void Foam::fv::thermalRunawaySource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    const Time& time(obr_.time());
    const scalarField& V(mesh_.V());
    const volScalarField& T(obr_.lookupObject<volScalarField>("T"));

    scalarField& heSource = eqn.source();

    scalar Tmean(0);
    scalar totalVolume(0);

    forAll(cells_, i)
    {
        const label celli = cells_[i];
        Tmean += T[celli]*V[celli];
        totalVolume += V[celli];
    }
    label nCells = cells_.size();
    reduce(Tmean, sumOp<scalar>());
    reduce(totalVolume, sumOp<scalar>());
    reduce(nCells, sumOp<label>());

    if (nCells>0)
    {
        Tmean /= totalVolume;
        if (debug)
            Info<< "Tmean is " << Tmean
                << ", nCells " << nCells << endl;

        if (Tmean>Tcritical_)
        {
            if (state_ == false)
            {
                TR_time_ = time.value();
                state_ = true;
                Info<< "Starting thermal runaway for "
                    << option::name_ << endl;
            }

            scalar tEff = time.value()-TR_time_;

            heSource -= HRR_->value(tEff)*V/VDash_;
            if (debug)
                Info<< "Heat release from table "
                    << HRR_->value(tEff) << endl;
        }
    }

    if (state_)
        setResult("state", state_);
    if (TR_time_>0)
        setResult("TR_time", TR_time_);
}


bool Foam::fv::thermalRunawaySource::write()
{
    return true;
}

bool Foam::fv::thermalRunawaySource::execute()
{
    return true;
}


bool Foam::fv::thermalRunawaySource::read(const dictionary& dict)
{
    stateFunctionObject::read(dict);
    state_ = getResult<bool>("state");
    TR_time_ = getResult<scalar>("TR_time");

    if (cellSetOption::read(dict))
    {
        volumeMode_ = volumeModeTypeNames_
            [
                coeffs_.lookupOrDefault<word>("volumeMode", "absolute")
            ];

        // Set volume normalisation
        if (volumeMode_ == vmAbsolute)
        {
            VDash_ = V_;
        }

        Tcritical_ = coeffs_.lookup<scalar>("Tcritical");
        HRR_ = Function1<scalar>::New("HRR", coeffs_);

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //


