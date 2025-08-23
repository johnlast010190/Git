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
    (c) 2011-2023 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "singleStepCombustion/singleStepCombustion.H"
#include "finiteVolume/fvm/fvmSup.H"
#include "materialModels/baseModels/materialModels.H"
#include "materialModels/materialTables/materialTables.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::combustionModels::singleStepCombustion::calculateqFuel()
{
    const scalar Wu = mixture_.W(fuelIndex_);
    const word& phaseName = thermo().phaseName();
    const word& hfName = hfModel::typeName;

    forAll(reaction_.lhs(), i)
    {
        const label speciei = reaction_.lhs()[i].index;
        const scalar stoichCoeff = reaction_.lhs()[i].stoichCoeff;
        specieStoichCoeffs_[speciei] = -stoichCoeff;
        const word specieName = mixture_.Y(speciei).name();
        const scalar hf =
            (basicThermo::dictName == basicThermo::matDictName)
          ? thermo().materials()(hfName, phaseName, specieName)[0]
          : mixture_.hf(speciei);

        qFuel_.value() += hf*mixture_.W(speciei)*stoichCoeff/Wu;
    }

    forAll(reaction_.rhs(), i)
    {
        const label speciei = reaction_.rhs()[i].index;
        const scalar stoichCoeff = reaction_.rhs()[i].stoichCoeff;
        specieStoichCoeffs_[speciei] = stoichCoeff;
        const word specieName = mixture_.Y(speciei).name();
        const scalar hf =
            (basicThermo::dictName == basicThermo::matDictName)
          ? thermo().materials()(hfName, phaseName, specieName)[0]
          : mixture_.hf(speciei);

        qFuel_.value() -= hf*mixture_.W(speciei)*stoichCoeff/Wu;
        specieProd_[speciei] = -1;
    }

    Info<< "Fuel heat of combustion :" << qFuel_.value() << endl;
}


void Foam::combustionModels::singleStepCombustion::massAndAirStoichRatios()
{
    const label O2Index = mixture_.species()["O2"];
    const scalar Wu = mixture_.W(fuelIndex_);
    const label defaultSpecie = mixture_.defaultSpecie();
    stoicRatio_ =
        (
            mixture_.W(defaultSpecie)*specieStoichCoeffs_[defaultSpecie]
          + mixture_.W(O2Index)*mag(specieStoichCoeffs_[O2Index])
        )
       /(Wu*mag(specieStoichCoeffs_[fuelIndex_]));

    s_ =
        (mixture_.W(O2Index)*mag(specieStoichCoeffs_[O2Index]))
       /(Wu*mag(specieStoichCoeffs_[fuelIndex_]));

    Info<< "stoichiometric air-fuel ratio :" << stoicRatio_.value() << nl
        << "stoichiometric oxygen-fuel ratio :" << s_.value() << endl;
}


void Foam::combustionModels::singleStepCombustion::calculateMaxProducts()
{
    scalar Wm = 0.0;
    scalar totalMol = 0.0;
    forAll(reaction_.rhs(), i)
    {
        label speciei = reaction_.rhs()[i].index;
        totalMol += mag(specieStoichCoeffs_[speciei]);
    }

    scalarList Xi(reaction_.rhs().size());

    forAll(reaction_.rhs(), i)
    {
        const label speciei = reaction_.rhs()[i].index;
        Xi[i] = mag(specieStoichCoeffs_[speciei])/totalMol;
        Wm += Xi[i]*mixture_.W(speciei);
    }

    forAll(reaction_.rhs(), i)
    {
        const label speciei = reaction_.rhs()[i].index;
        Yprod0_[speciei] =  mixture_.W(speciei)/Wm*Xi[i];
    }

    Info<< "Maximum products mass concentrations:" << nl;
    forAll(Yprod0_, i)
    {
        if (Yprod0_[i] > 0)
        {
            Info<< "    " << mixture_.species()[i] << ": " << Yprod0_[i] << nl;
        }
    }

    // Normalize the stoichiometric coeff to mass
    forAll(specieStoichCoeffs_, i)
    {
        specieStoichCoeffs_[i] =
            specieStoichCoeffs_[i]*mixture_.W(i)
           /(mixture_.W(fuelIndex_)*mag(specieStoichCoeffs_[fuelIndex_]));
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::combustionModels::singleStepCombustion::singleStepCombustion
(
    const word& modelType,
    const fluidMulticomponentThermo& thermo,
    const compressibleTurbulenceModel& turb,
    const word& combustionProperties
)
:
    combustionModel(modelType, thermo, turb, combustionProperties),
    mixture_(dynamic_cast<const basicSpecieMixture&>(this->thermo())),
    reaction_(mixture_.species(), this->subDict("reaction")),
    stoicRatio_(dimensionedScalar("stoicRatio", dimless, 0)),
    s_(dimensionedScalar("s", dimless, 0)),
    qFuel_(dimensionedScalar("qFuel", sqr(dimVelocity), 0)),
    specieStoichCoeffs_(mixture_.species().size(), 0.0),
    Yprod0_(mixture_.species().size(), 0.0),
    fres_(Yprod0_.size()),
    fuelIndex_(mixture_.species()[thermo.properties().lookup("fuel")]),
    specieProd_(Yprod0_.size(), 1),
    wFuel_
    (
        IOobject
        (
            this->thermo().phasePropertyName("wFuel"),
            this->mesh().time().timeName(),
            this->thermo().db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar(dimMass/dimVolume/dimTime, 0)
    ),
    semiImplicit_(readBool(this->coeffs_.lookup("semiImplicit")))
{
    forAll(fres_, fresI)
    {
        fres_.set
        (
            fresI,
            new volScalarField
            (
                IOobject
                (
                    "fres_" + mixture_.species()[fresI],
                    this->mesh().time().timeName(),
                    this->mesh()
                ),
                this->mesh(),
                dimensionedScalar(dimless, 0)
            )
        );
    }

    calculateqFuel();

    massAndAirStoichRatios();

    calculateMaxProducts();

    if (semiImplicit_)
    {
        Info<< "Combustion mode: semi-implicit" << endl;
    }
    else
    {
        Info<< "Combustion mode: explicit" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::combustionModels::singleStepCombustion::~singleStepCombustion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::fvScalarMatrix> Foam::combustionModels::singleStepCombustion::R
(
    volScalarField& Y
) const
{
    const label specieI = this->thermo().composition().species()[Y.member()];

    volScalarField wSpecie(wFuel_*specieStoichCoeffs()[specieI]);

    if (semiImplicit_)
    {
        const label fNorm = specieProd()[specieI];
        const volScalarField fres(this->fres(specieI));
        wSpecie /= max(fNorm*(Y - fres), scalar(1e-2));

        return -fNorm*wSpecie*fres + fNorm*fvm::Sp(wSpecie, Y);
    }
    else
    {
        return wSpecie + fvm::Sp(0.0*wSpecie, Y);
    }
}


Foam::tmp<Foam::volScalarField>
Foam::combustionModels::singleStepCombustion::Qdot() const
{
    const label fuelI = fuelIndex();
    volScalarField& YFuel =
        const_cast<volScalarField&>(this->thermo().composition().Y(fuelI));

    return -qFuel()*(R(YFuel) & YFuel);
}


bool Foam::combustionModels::singleStepCombustion::read()
{
    if (combustionModel::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}


void Foam::combustionModels::singleStepCombustion::fresCorrect()
{
    const label O2Index = mixture_.species()["O2"];
    const volScalarField& YFuel = mixture_.Y()[fuelIndex_];
    const volScalarField& YO2 = mixture_.Y()[O2Index];

    // reactants
    forAll(reaction_.lhs(), i)
    {
        const label speciei = reaction_.lhs()[i].index;
        if (speciei == fuelIndex_)
        {
            fres_[speciei] = max(YFuel - YO2/s_, scalar(0));
        }
        else if (speciei == O2Index)
        {
            fres_[speciei] = max(YO2 - YFuel*s_, scalar(0));
        }
    }

    // products
    forAll(reaction_.rhs(), i)
    {
        const label speciei = reaction_.rhs()[i].index;
        if (speciei != mixture_.defaultSpecie())
        {
            forAll(fres_[speciei], celli)
            {
                if (fres_[fuelIndex_][celli] > 0.0)
                {
                    // rich mixture
                    fres_[speciei][celli] =
                        Yprod0_[speciei]
                      * (1.0 + YO2[celli]/s_.value() - YFuel[celli]);
                }
                else
                {
                    // lean mixture
                    fres_[speciei][celli] =
                        Yprod0_[speciei]
                      * (
                            1.0
                          - YO2[celli]/s_.value()*stoicRatio_.value()
                          + YFuel[celli]*stoicRatio_.value()
                        );
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
