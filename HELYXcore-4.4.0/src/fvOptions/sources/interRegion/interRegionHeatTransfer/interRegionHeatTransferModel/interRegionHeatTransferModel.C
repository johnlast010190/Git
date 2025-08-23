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
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "interRegionHeatTransferModel.H"
#include "basicThermo/basicThermo.H"
#include "finiteVolume/fvm/fvmSup.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchFields.H"
#include "finiteVolume/fvc/fvcVolumeIntegrate.H"
#include "cfdTools/general/fvOptions/fvOptionList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(interRegionHeatTransferModel, 0);
}
}


// * * * * * * * * * * * *  Protected member functions * * * * * * * * * * * //

void Foam::fv::interRegionHeatTransferModel::setNbrModel()
{
    if (!firstIter_)
    {
        return;
    }

    const fvMesh& nbrMesh = mesh_.time().lookupObject<fvMesh>(nbrRegionName_);

    const optionList& fvOptions = nbrMesh.lookupObject<optionList>("fvOptions");

    bool nbrModelFound = false;

    forAll(fvOptions, i)
    {
        if (fvOptions[i].name() == nbrModelName_)
        {
            nbrModel_ = &const_cast<interRegionHeatTransferModel&>
            (
                refCast<const interRegionHeatTransferModel>(fvOptions[i])
            );
            nbrModelFound = true;
            break;
        }
    }

    if (!nbrModelFound)
    {
        FatalErrorInFunction
            << "Neighbour model not found" << nbrModelName_
            << " in region " << nbrMesh.name() << nl
            << exit(FatalError);
    }

    firstIter_ = false;

    // Set nbr model's nbr model to avoid construction order problems
    nbrModel_->setNbrModel();
}


void Foam::fv::interRegionHeatTransferModel::correct()
{
    if (master_)
    {
        if (mesh_.time().timeIndex() != timeIndex_)
        {
            const volScalarField& T = obr_.lookupObject<volScalarField>(TName_);
            tmp<volScalarField> tTmapped = getTmapped(T);
            volScalarField& Tmapped = tTmapped.ref();
            const scalar deltaT = fvc::domainIntegrate(T - Tmapped).value()/meshInterp().V();
            calculateHtc(deltaT);
        }
    }
    else
    {
        nbrModel().correct();
        interpolate(nbrModel().htc(), htc_);
    }
}


Foam::tmp<Foam::volScalarField>
Foam::fv::interRegionHeatTransferModel::getTmapped
(
    const volScalarField& T
) const
{
    tmp<volScalarField> tTmapped(volScalarField::New(type() + "Tmapped", T));
    volScalarField& Tmapped = tTmapped.ref();

    const fvMesh& nbrMesh = mesh_.time().lookupObject<fvMesh>(nbrRegionName_);

    const volScalarField& Tnbr =
        nbrMesh.lookupObject<volScalarField>(TNbrName_);

    interpolate
    (
        // Nbr temp relative to this side's reference
        (
            Tnbr.primitiveField()
          - basicThermo::TRefIfFound(nbrMesh).value()
          + basicThermo::TRefIfFound(Tmapped.db()).value()
        )(),
        Tmapped.primitiveFieldRef()
    );

    return tTmapped;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::interRegionHeatTransferModel::interRegionHeatTransferModel
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    interRegionOption
    (
        name,
        modelType,
        dict,
        obr
    ),
    nbrModelName_(coeffs_.lookup("nbrModel")),
    nbrModel_(nullptr),
    firstIter_(true),
    timeIndex_(-1),
    htc_
    (
        IOobject
        (
            type() + "Htc",
            mesh_.time().timeName(),
            obr_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar
        (
            "htc",
            dimEnergy/dimTime/dimTemperature/dimVolume,
            0.0
        ),
        zeroGradientFvPatchScalarField::typeName
    ),
    semiImplicit_(false),
    TName_(coeffs_.lookupOrDefault<word>("T", "T")),
    TNbrName_(coeffs_.lookupOrDefault<word>("TNbr", "T")),
    reportOn_(coeffs_.lookupOrDefault<Switch>("reportOn", false)),
    totalEnergy_(0.0),
    targetEnergy_(0.0),
    conservative_(coeffs_.lookupOrDefault<Switch>("conservative", true))
{
    if (active())
    {
        semiImplicit_ = coeffs_.lookup<bool>("semiImplicit");
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::interRegionHeatTransferModel::~interRegionHeatTransferModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::interRegionHeatTransferModel::sourceFields
(
    wordList& fieldNames
)
{
    if (active())
    {
        if (coeffs_.found("fields"))
        {
            fieldNames = coeffs_.lookup<wordList>("fields");
        }
        else
        {
            const basicThermo& thermo =
                obr_.lookupObject<basicThermo>(basicThermo::dictName);
            fieldNames = wordList(1, thermo.heT().name());
        }
    }
}


void Foam::fv::interRegionHeatTransferModel::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    setNbrModel();

    correct();
    timeIndex_ = mesh_.time().timeIndex();

    const volScalarField& heT = eqn.psi();

    const volScalarField& T = obr_.lookupObject<volScalarField>(TName_);
    tmp<volScalarField> tTmapped = getTmapped(T);
    volScalarField& Tmapped = tTmapped.ref();

    // scale total energy to ensure conservation
    scalar scalingFactor(1.0);

    if (conservative_)
    {
        totalEnergy_ = fvc::domainIntegrate(htc_*(T - Tmapped)).value();
        targetEnergy_ = getTargetEnergy();

        if
        (
            mag(nbrModel().totalEnergy()) > VSMALL
         && mag(totalEnergy()) > VSMALL
        ) // avoid scaling for first time step and if source is small
        {
            if (timeIndex_ > nbrModel().timeIndex_)
            {
                scalingFactor = mag(targetEnergy()/totalEnergy());
            }
            else
            {
                scalingFactor =
                    mag(nbrModel().targetEnergy()/totalEnergy());
            }
        }
        Info<< "scaling target energy by " << scalingFactor << endl;
    }

    if (reportOn_)
    {
        Info<< "Volumetric integral of htc: "
            << fvc::domainIntegrate(htc_).value()
            << endl;

        if (!conservative_) // avoid re-calculation
        {
            totalEnergy_ = fvc::domainIntegrate(htc_*(T - Tmapped)).value();
        }
        const fvMesh& nbrMesh = mesh_.time().lookupObject<fvMesh>(nbrRegionName_);

        Info<< "Energy exchange [J/s] from region " << nbrMesh.name()
            << " To " << mesh_.name() << " : " <<  scalingFactor*totalEnergy()
            << endl;

        if (mesh_.time().writeTime())
        {
            Tmapped.write();

            // Htc is made on construction which bypasses virtual dispatch
            // this corrects htc name for writing.
            if (htc_.name() != word(type() + "Htc"))
            {
                htc_.rename(type() + "Htc");
            }
            htc_.write();
        }
    }

    if (semiImplicit_)
    {
        if (heT.dimensions() == dimEnergy/dimMass)  //he
        {
            if (obr_.foundObject<basicThermo>(basicThermo::dictName))
            {
                const basicThermo& thermo =
                   obr_.lookupObject<basicThermo>(basicThermo::dictName);

                volScalarField htcByCpv(htc_/thermo.Cpv());

                eqn +=
                    scalingFactor*htc_*(Tmapped - T)
                  + htcByCpv*heT - fvm::Sp(htcByCpv, heT);
            }
            else
            {
                FatalErrorInFunction
                    << " on mesh " << mesh_.name()
                    << " could not find object basicThermo."
                    << " The available objects are: "
                    << mesh_.names()
                    << exit(FatalError);
            }
        }
        else if (heT.dimensions() == dimTemperature) //T
        {
            eqn +=
                scalingFactor*htc_*Tmapped - fvm::Sp(scalingFactor*htc_, heT);
        }
    }
    else
    {
        eqn += scalingFactor*htc_*(Tmapped - T);
    }
}


void Foam::fv::interRegionHeatTransferModel::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    addSup(eqn, fieldi);
}


Foam::scalar Foam::fv::interRegionHeatTransferModel::getTargetEnergy()
{
    // for first time step, take whatever side is evaluated first as target
    if
    (
        mesh_.time().restartTimeIndex() == 1
     && timeIndex_ > nbrModel().timeIndex_
    )
    {
        return totalEnergy();
    }

    return (totalEnergy() - nbrModel().totalEnergy())/2.0;
}


// ************************************************************************* //
