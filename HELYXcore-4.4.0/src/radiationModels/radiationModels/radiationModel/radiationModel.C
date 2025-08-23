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
    (c) 2011-2022 OpenFOAM Foundation
    (c) 2016 OpenCFD Ltd.
    (c) 2016-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "radiationModels/radiationModel/radiationModel.H"
#include "absorptionEmissionModels/absorptionEmissionModel/absorptionEmissionModel.H"
#include "scatterModels/scatterModel/scatterModel.H"
#include "sootModels/sootModel/sootModel.H"
#include "finiteVolume/fvm/fvmSup.H"
#include "basicThermo/basicThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(radiationModel, 0);
    defineRunTimeSelectionTable(radiationModel, T);
    defineRunTimeSelectionTable(radiationModel, dictionary);
}

const Foam::word Foam::radiationModel::dictName = "radiationProperties";
const Foam::word Foam::radiationModel::externalRadHeatFieldName_ = "QrExt";


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::IOobject Foam::radiationModel::createIOobject
(
    const objectRegistry& obr
) const
{
    IOobject io
    (
        dictName,
        obr.time().constant(),
        obr,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
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


void Foam::radiationModel::initialise()
{
    solverFreq_ = max(1, lookupOrDefault<label>("solverFreq", 1));

    participating_ = lookupOrDefault<Switch>("participating", true);

    absorptionEmission_.reset
    (
        radiationModels::absorptionEmissionModel::New(*this, mesh_).ptr()
    );

    scatter_.reset(radiationModels::scatterModel::New(*this, mesh_).ptr());

    soot_.reset(radiationModels::sootModel::New(*this, mesh_).ptr());

    transmissivity_.reset(radiationModels::transmissivityModel::New(*this, mesh_).ptr());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiationModel::radiationModel
(
    const volScalarField& T,
    const dimensionedScalar& TRef
)
:
    IOdictionary
    (
        IOobject
        (
            dictName,
            T.time().constant(),
            T.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(T.mesh()),
    time_(T.time()),
    T_(T),
    TRef_(TRef),
    participating_(true),
    coeffs_(dictionary::null),
    solverFreq_(0),
    firstIter_(true),
    absorptionEmission_(nullptr),
    scatter_(nullptr),
    soot_(nullptr),
    transmissivity_(nullptr)
{}


Foam::radiationModel::radiationModel
(
    const word& type,
    const volScalarField& T,
    const dimensionedScalar& TRef
)
:
    IOdictionary(createIOobject(T.db())),
    mesh_(T.mesh()),
    time_(T.time()),
    T_(T),
    TRef_(TRef),
    participating_(lookupOrDefault("participating", true)),
    coeffs_(subOrEmptyDict(type + "Coeffs")),
    solverFreq_(1),
    firstIter_(true),
    absorptionEmission_(nullptr),
    scatter_(nullptr),
    soot_(nullptr),
    transmissivity_(nullptr)
{
    initialise();
}


Foam::radiationModel::radiationModel
(
    const word& type,
    const dictionary& dict,
    const volScalarField& T,
    const dimensionedScalar& TRef
)
:
    IOdictionary
    (
        IOobject
        (
            dictName,
            T.time().constant(),
            T.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dict
    ),
    mesh_(T.mesh()),
    time_(T.time()),
    T_(T),
    TRef_(TRef),
    participating_(lookupOrDefault("participating", true)),
    coeffs_(subOrEmptyDict(type + "Coeffs")),
    solverFreq_(1),
    firstIter_(true),
    absorptionEmission_(nullptr),
    scatter_(nullptr),
    soot_(nullptr),
    transmissivity_(nullptr)
{
    initialise();
}


// * * * * * * * * * * * * * * * * Destructor    * * * * * * * * * * * * * * //

Foam::radiationModel::~radiationModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool Foam::radiationModel::correct()
{
    bool didRecalc = false;
    if (firstIter_ || (time_.timeIndex() % solverFreq_ == 0))
    {
        calculate();
        firstIter_ = false;
        didRecalc = true;
    }

    if (!soot_.empty())
    {
        soot_->correct();
    }

    return didRecalc;
}


bool Foam::radiationModel::read()
{
    if (regIOobject::read())
    {
        coeffs_ = subOrEmptyDict(type() + "Coeffs");

        participating_ = lookupOrDefault<Switch>("participating", true);
        solverFreq_ = lookupOrDefault<label>("solverFreq", 1);
        solverFreq_ = max(1, solverFreq_);

        return true;
    }
    else
    {
        return false;
    }
}


Foam::tmp<Foam::fvScalarMatrix> Foam::radiationModel::Sh
(
    const basicThermo& thermo,
    const volScalarField& he
) const
{
    const volScalarField& Cpv = thermo.Cpv();
    const volScalarField T3(pow3(T_ + TRef_));

    return
    (
        Ru()
      - fvm::Sp(4.0*Rp()*T3/Cpv, he)
      - Rp()*T3*(T_ + TRef_ - 4.0*he/Cpv)
    );
}


Foam::tmp<Foam::fvScalarMatrix> Foam::radiationModel::ST
(
    const dimensionedScalar& rhoCp,
    volScalarField& T
) const
{
    const volScalarField C(Rp()*pow3(T_ + TRef_)/rhoCp);
    return
    (
        Ru()/rhoCp
      - fvm::Sp(4.0*C, T)
      + 3.0*C*(T_+TRef_)
      - 4.0*C*TRef_
    );
}


const Foam::radiationModels::absorptionEmissionModel&
Foam::radiationModel::absorptionEmission() const
{
    if (!absorptionEmission_.valid())
    {
        FatalErrorInFunction
            << "Requested radiation absorptionEmission model, but model is "
            << "not activate" << abort(FatalError);
    }

    return absorptionEmission_();
}


const Foam::radiationModels::sootModel&
Foam::radiationModel::soot() const
{
    if (!soot_.valid())
    {
        FatalErrorInFunction
            << "Requested radiation sootModel model, but model is "
            << "not activate" << abort(FatalError);
    }

    return soot_();
}


const Foam::radiationModels::transmissivityModel&
Foam::radiationModel::transmissivity() const
{
    if (!transmissivity_.valid())
    {
        FatalErrorInFunction
            << "Requested radiation sootModel model, but model is "
            << "not activated" << abort(FatalError);
    }

    return transmissivity_();
}


// ************************************************************************* //
