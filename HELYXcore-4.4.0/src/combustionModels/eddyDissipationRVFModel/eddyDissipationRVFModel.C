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
    (c) 2011-2012 OpenFOAM Foundation
    (c) 2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "eddyDissipationRVFModel/eddyDissipationRVFModel.H"
#include "turbulenceModel.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "fields/volFields/volFields.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "cfdTools/general/include/fvCFD.H"
#include "global/etcFiles/etcFiles.H"
#include "db/IOstreams/Fstreams/IFstream.H"
#include "materialModels/materialTables/materialTables.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{
    defineTypeNameAndDebug(eddyDissipationRVFModel, 0);
    addToRunTimeSelectionTable
    (
        combustionModel,
        eddyDissipationRVFModel,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::combustionModels::eddyDissipationRVFModel::eddyDissipationRVFModel
(
    const word& modelType,
    const fluidMulticomponentThermo& thermo,
    const compressibleTurbulenceModel& turb,
    const word& combustionProperties
)
:
    singleStepCombustion(modelType, thermo, turb, combustionProperties),
    C_(this->coeffs().template lookup<scalar>("CEDC")),
    Cd_(this->coeffs().template lookup<scalar>("C_Diff")),
    Cstiff_(this->coeffs().template lookup<scalar>("CStiff")),
    tExt_
    (
        this->coeffs().template lookupOrDefault<scalar>("ExtinctionStart", 5.0)
    ),
    Cevap_(this->coeffs().template lookupOrDefault<scalar>("Cevap", 0.5)),
    ZN_(this->coeffs().template lookupOrDefault<scalar>("ZN", 10)),
    cKa_(this->coeffs().template lookupOrDefault<scalar>("cKa", 1)),
    cKapa_(this->coeffs().template lookupOrDefault<scalar>("cKapa", 1)),
    XrExt_(this->coeffs().template lookupOrDefault<scalar>("XrExt", 0)),
    TFuel_(this->coeffs().template lookupOrDefault<scalar>("TFuel", 293.15)),
    TAir_(this->coeffs().template lookupOrDefault<scalar>("TAir", 293.15)),
    TadAir_(this->coeffs().template lookupOrDefault<scalar>("TadAir", 2400)),
    SLC1_(this->coeffs().template lookupOrDefault<scalar>("SLC1", 56)),
    SLC2_(this->coeffs().template lookupOrDefault<scalar>("SLC2", 11.4)),
    fvSootAir_
    (
        this->coeffs().template lookupOrDefault<scalar>("fvSootAir", 2.3)
    ),   //- in PPM
    O2Soot_(this->coeffs().template lookupOrDefault<scalar>("O2Soot", 0.137)),
    fixedXr_(this->coeffs().template lookupOrDefault<scalar>("fixedXr", false)),
    RVFModelActivated_(false),
    variableXr_(false),
    Ka_
    (
        IOobject("Ka", this->mesh().time().timeName(), this->mesh()),
        this->mesh(),
        dimensionedScalar(dimless, 0)
    ),
    KaExt_
    (
        IOobject("KaExt", this->mesh().time().timeName(), this->mesh()),
        this->mesh(),
        dimensionedScalar(dimless, 0)
    ),
    KaMixed_
    (
        IOobject("KaMixed", this->mesh().time().timeName(), this->mesh()),
        this->mesh(),
        dimensionedScalar(dimless, 0)
    ),
    KaExtMixed_
    (
        IOobject("KaExtMixed", this->mesh().time().timeName(), this->mesh()),
        this->mesh(),
        dimensionedScalar(dimless, 0)
    ),
    RVF_
    (
        IOobject("RVF", this->mesh().time().timeName(), this->mesh()),
        this->mesh(),
        dimensionedScalar(dimless, 1.0)
    ),
    PV_
    (
        IOobject
        (
            "PV",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar(dimless, 0)
    ),
    FExt_
    (
        IOobject
        (
            "FExt",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar(dimless, 1.0)
    ),
    Fig_
    (
        IOobject
        (
            "Fig",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar(dimless, 0),
	    zeroGradientFvPatchScalarField::typeName
    ),
    epsSgs_
    (
        IOobject("epsSgs", this->mesh().time().timeName(), this->mesh()),
        this->mesh(),
        dimensionedScalar(dimLength*dimLength/dimTime/dimTime/dimTime, 0)
    ),
    epsG_
    (
        IOobject("epsilon", this->mesh().time().timeName(), this->mesh()),
        this->mesh(),
        dimensionedScalar(dimLength*dimLength/dimTime/dimTime/dimTime, 0)
    ),
    SL_
    (
        IOobject("SL", this->mesh().time().timeName(), this->mesh()),
        this->mesh(),
        dimensionedScalar(dimLength/dimTime, 0)
    ),
    SLMixed_
    (
        IOobject("SLMixed", this->mesh().time().timeName(), this->mesh()),
        this->mesh(),
        dimensionedScalar(dimLength/dimTime, 0)
    ),
    dQFuel_
    (
        IOobject("dQFuel", this->mesh().time().timeName(), this->mesh()),
        this->mesh(),
        dimensionedScalar(dimEnergy/dimTime/dimVolume, 0)
    ),
    dQFstar_
    (
        IOobject("dQFstar", this->mesh().time().timeName(), this->mesh()),
        this->mesh(),
        dimensionedScalar(dimEnergy/dimTime/dimVolume, 0)
    ),
    deltaFlame_
    (
        IOobject("deltaFlame", this->mesh().time().timeName(), this->mesh()),
        this->mesh(),
        dimensionedScalar(dimLength, 0)
    ),
    deltaFlameMixed_
    (
        IOobject
        (
            "deltaFlameMixed",
            this->mesh().time().timeName(),
            this->mesh()
        ),
        this->mesh(),
        dimensionedScalar(dimLength, 0)
    ),
    alphaLoss_
    (
        IOobject("alphaLoss", this->mesh().time().timeName(), this->mesh()),
        this->mesh(),
        dimensionedScalar(dimless, 0)
    ),
    ExpR_
    (
        IOobject
        (
            "ExpR",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar(dimless, 1)
    ),
    ExpRMixed_
    (
        IOobject
        (
            "ExpRMixed",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar(dimless, 1)
    ),
    ER_
    (
        IOobject
        (
            "ER",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar(dimless, 1)
    ),
    Beta_
    (
        IOobject("Beta", this->mesh().time().timeName(), this->mesh()),
        this->mesh(),
        dimensionedScalar(dimless, 10)
    ),
    Gama_
    (
        IOobject("Gama", this->mesh().time().timeName(), this->mesh()),
        this->mesh(),
        dimensionedScalar(dimless, 10)
    ),
    GamaMixed_
    (
        IOobject("GamaMixed", this->mesh().time().timeName(), this->mesh()),
        this->mesh(),
        dimensionedScalar(dimless, 10)
    ),
    ExtNumber_
    (
        IOobject("ExtNumber", this->mesh().time().timeName(), this->mesh()),
        this->mesh(),
        dimensionedScalar(dimless, 10)
    ),
    ExtNumberMixed_
    (
        IOobject
        (
            "ExtNumberMixed",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar(dimless, 10)
    ),
    Tad_
    (
        IOobject
        (
            "Tad",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar(dimTemperature, 293)
    ),
    TadMixed_
    (
        IOobject
        (
            "TadMixed",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar(dimTemperature, 293)
    ),
    XO2Local_
    (
        IOobject("XO2Local", this->mesh().time().timeName(), this->mesh()),
        this->mesh(),
        dimensionedScalar(dimless, 0.21)
    ),
    XrFlame_
    (
        IOobject("XrFlame", this->mesh().time().timeName(), this->mesh()),
        this->mesh(),
        dimensionedScalar(dimless, 0.2)
    ),
    fvSoot_
    (
        IOobject("fvSoot", this->mesh().time().timeName(), this->mesh()),
        this->mesh(),
        dimensionedScalar(dimless, 2.0)
    ),
    QdotRad_
    (
        IOobject("QdotRad", this->mesh().time().timeName(), this->mesh()),
        this->mesh(),
        dimensionedScalar(dimEnergy/pow3(dimLength)/dimTime, 0)
    ),
    WFstar_
    (
        IOobject("WFstar", this->mesh().time().timeName(), this->mesh()),
        this->mesh(),
        dimensionedScalar(dimMass/pow3(dimLength)/dimTime, 0)
    )
{}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::combustionModels::eddyDissipationRVFModel::~eddyDissipationRVFModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::combustionModels::eddyDissipationRVFModel::rtTurb() const
{
    return
        C_*this->turbulence().epsilon()
       /max
        (
            this->turbulence().k(),
            dimensionedScalar(dimVelocity*dimVelocity, SMALL)
        );
}


Foam::tmp<Foam::volScalarField>
Foam::combustionModels::eddyDissipationRVFModel::rtDiff() const
{
    const volScalarField& YO2 = this->thermo().composition().Y("O2");
    const compressible::LESModel& lesModel =
        YO2.db().lookupObject<compressible::LESModel>
	    (
	        turbulenceModel::propertiesName
	    );

    const dimensionedScalar Df("Df", dimLength, 0.01);
    return
        Cd_*(this->thermo().kappa()/this->thermo().Cp())
       /this->rho()/sqr(min(Df, lesModel.delta()));
}


Foam::tmp<Foam::volScalarField>
Foam::combustionModels::eddyDissipationRVFModel::rtSpread() const
{
    const volScalarField& YO2 = this->thermo().composition().Y("O2");
    const compressible::LESModel& lesModel =
        YO2.db().lookupObject<compressible::LESModel>
        (
            turbulenceModel::propertiesName
        );

    return SLMixed_/lesModel.delta();
}


void Foam::combustionModels::eddyDissipationRVFModel::correct()
{
    this->wFuel_.forceAssign
    (
        dimensionedScalar(dimMass/pow3(dimLength)/dimTime, 0)
    );

    this->fresCorrect();

    const label fuelI = this->fuelIndex();
    const volScalarField& YFuel = this->thermo().composition().Y()[fuelI];
    const volScalarField& YFstar = this->thermo().composition().Y("Fstar");
    const dimensionedScalar s = this->s().value();

    if (this->thermo().composition().contains("O2"))
    {
        const volScalarField& YO2 = this->thermo().composition().Y("O2");
        volScalarField rt(max(rtTurb(),rtDiff()));
	    calculateFlameTemperature();
	    calculateReactiveVolume();
        this->wFuel_.forceAssign
        (
            this->rho()*min(YFuel, YO2/s)
           /this->mesh().time().deltaT()/Cstiff_
           *(1 - exp(-Cstiff_*this->mesh().time().deltaT()*rt))
        );
        this->WFstar_.forceAssign
        (
            this->rho()*min(YFstar, YO2/s)
           /this->mesh().time().deltaT()/Cstiff_
           *(1 - exp(- Cstiff_*this->mesh().time().deltaT()*rt))
        );
    }
	dQFuel_ = this->QdotFuel();
	dQFstar_ = this->QdotFstar();
    QdotRad_ = XrFlame_*(dQFuel_ + dQFstar_);
}


void Foam::combustionModels::eddyDissipationRVFModel::calculateReactiveVolume()
{
    const label fuelI = this->fuelIndex();
    const volScalarField& YFuel = this->thermo().composition().Y()[fuelI];
    const volScalarField& YFstar = this->thermo().composition().Y("Fstar");
    const volScalarField& YO2 = this->thermo().composition().Y("O2");

    tmp<volScalarField> tnut(this->turbulence().nut());
    const volScalarField& nut = tnut();
    epsSgs_ = this->turbulence().epsilon();

    forAll(YO2, celli)
    {
        //- Initialize parameters for normal fuel
        SL_[celli] = 0;
        deltaFlame_[celli] = 0;
        Ka_[celli] = 0;
        KaExt_[celli] = 1;
        Gama_[celli] = 6.0;
        RVF_[celli] = 0;
        ExtNumber_[celli] = 0;
        alphaLoss_[celli] = 0;
        Beta_[celli] = 20.0;
        FExt_[celli] = 0;

        //- Initialize parameters for quenched fuel
        SLMixed_[celli] = 0;
        deltaFlameMixed_[celli] = 0;
        KaMixed_[celli] = 0;
        KaExtMixed_[celli] = 1;
        GamaMixed_[celli] = 6.0;
        ExtNumberMixed_[celli] = 0;
        Fig_[celli] = 0;

        //- Calculate radiative heat loss term
        scalar krad = 0.7;
        scalar kai = pow(Tad_[celli]/TAir_, 1.75)*1.4e-5/0.75;
        epsG_[celli] = epsSgs_[celli]*(kai + nut[celli])/(1.0e-8 + nut[celli]);
        scalar BetaFuel = 10;
        scalar BetaFstar = 10;

        //- Local flame radiant fraction
        XrFlame_[celli] = 0;
        scalar XrFuel = 0;
        scalar XrFstar = 0;

        if (YFuel[celli] > 1.0e-4 && YO2[celli] > 1.0e-4)
        {
            SL_[celli] = max(0.0001, SLC1_*exp(-1000.0*SLC2_/Tad_[celli]));
            deltaFlame_[celli] = min(0.1, kai/SL_[celli]);
            Gama_[celli] = 6.0*cKapa_*pow((Tad_[celli]/TAir_), 1.75);

            BetaFuel =
                max
                (
                    6.0,
                    min
                    (
                        20.0,
                        ZN_
                       *sqr(TadAir_/Tad_[celli])
                       *(Tad_[celli] - TAir_)
                       /(TadAir_ - TAir_)
                    )
                );
            Beta_[celli] = BetaFuel;
            scalar Xext =
                (
                    (BetaFuel + 0.667)
                  + sqrt(sqr(BetaFuel + 0.667) - 6.667*BetaFuel)
                )/(2*BetaFuel);

            scalar Cka = sqr(ExpR_[celli]*ZN_)/Gama_[celli]/exp(BetaFuel);
            scalar Te0 = Xext*Tad_[celli];
            scalar xo = XO2Local_[celli];
            scalar kradgas =
                233.67*xo*xo*xo - 110.62*xo*xo + 12.49*xo + 0.3849;

            // Calculate local radiant fraction
            scalar kradTotal = kradgas + 1.803e-3*fvSoot_[celli]*Te0;
            alphaLoss_[celli] =
                min
                (
                    1.0,
                    5.33*5.67e-8*kradTotal*sqr(deltaFlame_[celli])*pow(Te0, 3)
                   /(cKapa_*kai*(1.2*TAir_/Te0)*530*pow(Te0, 0.1131))
                );
            XrFuel = min(0.5, alphaLoss_[celli]/(1.0 + alphaLoss_[celli]));

            Ka_[celli] =
                min
                (
                    10.0,
                    cKa_*sqr(deltaFlame_[celli])
                   *sqrt(epsG_[celli])/pow(kai, 1.5)
                );
            KaExt_[celli] =
                pow
                (
                    max
                    (
                        0.0,
                        Cka*(1.0 - Xext)*exp(BetaFuel*Xext)/pow(Xext, 5.0/3.0)
                      - alphaLoss_[celli]*pow(Xext, 4.0/3.0)
                    ),
                    1.5
                );
            ExtNumber_[celli] = KaExt_[celli] - 1.0/Xext;

            // Purely thermal quenching
            if (ExtNumber_[celli] < 0)
            {
                RVF_[celli] = 0;
                FExt_[celli] = 1;
            }
            // Extinction due to high strain and thermal
            else
            {
                if (Ka_[celli] > KaExt_[celli])
                {
                    RVF_[celli] = 0;
                    FExt_[celli] = 1;
                }
                // Healthy flame
                else if (Ka_[celli] < 1.1)
                {
                    RVF_[celli] = 1.0;
                    FExt_[celli] = 0.0;
                }
                else
                {
                    // Close to quenching limit
                    if (KaExt_[celli] - Ka_[celli] < 0.05)
                    {
                   	    RVF_[celli] = 0.3;
                        FExt_[celli] = 0.7;
                    }
                    // Partial extinction
                    else
                    {
                        scalar X1 = Xext - 0.05;
                        scalar X2 = Xext + 0.05;
                        scalar X3 = Xext - 0.05;
                        scalar iterKa = 0;
                        while (iterKa < 10)
                        {
                            iterKa = iterKa + 1;
                            X1 =
                                log
                                (
                                    X1
                                   *(
                                        pow(Ka_[celli]*X1, 0.667)
                                      + alphaLoss_[celli]*sqr(X1)
                                    )/(Cka*(1-X1))
                                )/BetaFuel;

                            scalar YX2 =
                                exp(BetaFuel*X2)*Cka*pow(X2, 0.333)
                               /(
                                    pow(Ka_[celli], 0.667)
                                  + alphaLoss_[celli]*pow(X2, 1.333)
                                );

                            X2 = (sqrt(sqr(YX2) + 4.0*YX2) - YX2)/2.0;
                            X3 =
                                log
                                (
                                    X3*(1 + alphaLoss_[celli]*sqr(X3))
                                   /(Cka*(1-X3))
                                )
                               /BetaFuel;
                        }
                        RVF_[celli] =
                            max
                            (
                                0.0,
                                min
                                (
                                    1.0,
                                    (pow(X3/X1, 5.0) - pow(X3, 5.0))
                                   /(1.0 - pow(X3, 5.0))
                                )
                            );
                        FExt_[celli] = 1.0 - RVF_[celli];
                    }
                }
            }

            // Let pure fuel be flammable
            if (YO2[celli] < 1.0e-4)
            {
                RVF_[celli] = 1.0;
                FExt_[celli] = 0.0;
            }
        }

        if (YFstar[celli] > 1.0e-4)
        {
            scalar kaiMixed = pow(TadMixed_[celli]/TAir_, 1.75)*1.4e-5/0.75;
            BetaFstar =
                max
                (
                    6.0,
                    min
                    (
                        20.0,
                        ZN_*sqr(TadAir_/TadMixed_[celli])
                       *(TadMixed_[celli] - TAir_)
                       /(TadAir_ - TAir_)
                    )
                );
            SLMixed_[celli] =
                max(0.0001, SLC1_*exp(-1000.0*SLC2_/TadMixed_[celli]));
            deltaFlameMixed_[celli] = min(0.1, kaiMixed/SLMixed_[celli]);
            GamaMixed_[celli] = 6.0*cKapa_*pow((TadMixed_[celli]/TAir_), 1.75);

            scalar XextMixed =
                (
                    (BetaFstar + 0.667)
                  + sqrt(sqr(BetaFstar + 0.667) - 6.667*BetaFstar)
                )/(2*BetaFstar);

            scalar CkaMixed =
                sqr(ExpRMixed_[celli]*ZN_)/GamaMixed_[celli]/exp(BetaFstar);

            scalar Te0Mixed = XextMixed*TadMixed_[celli];
            scalar CalphaMixed =
                5.33*5.67e-8*krad*sqr(deltaFlameMixed_[celli])*pow(Te0Mixed, 3)
               /(
                    cKapa_*kaiMixed*(1.2*TAir_/Te0Mixed)
                   *530*pow(Te0Mixed, 0.1131)
                );

            KaMixed_[celli] =
                min
                (
                    10.,
                    cKa_*sqr(deltaFlameMixed_[celli])*sqrt(epsG_[celli])
                   /pow(kaiMixed, 1.5)
                );
            KaExtMixed_[celli] =
                pow
                (
                    max
                    (
                        0.0,
                        CkaMixed*(1.0 - XextMixed)*exp(BetaFstar*XextMixed)
                       /pow(XextMixed, 5.0/3.0)
                      - CalphaMixed*pow(XextMixed, 4.0/3.0)
                    ),
                    1.5
                );
            ExtNumberMixed_[celli] = KaExtMixed_[celli] - 1.0/XextMixed;

            //- Simplified treatment for re-ignition
            if (ExtNumberMixed_[celli] < 0)
            {
        	    Fig_[celli] = 0;
            }
            else
            {
                if (KaMixed_[celli] > KaExtMixed_[celli])
                {
                    Fig_[celli] = 0;
                }
                else
                {
                    Fig_[celli] = 1.0;
                }
            }
        }

        //- Local flame radiant fraction contributed from normal and quenched fuel
        XrFlame_[celli] =
            (
                XrFuel*YFuel[celli]*RVF_[celli]
              + XrFstar*YFstar[celli]*Fig_[celli]
            )
           /(
                YFuel[celli]*RVF_[celli]
              + YFstar[celli]*Fig_[celli] + 1.0e-8
            );

        // Activation flame extinction model after tExt
        if (!RVFModelActivated_)
        {
            RVF_[celli] = 1.0;
            Fig_[celli] = 1.0;
        }
    }

    // Activation flame extinction model after tExt
    if (this->mesh().time().value() < tExt_)
    {
        RVFModelActivated_ = false;
    }
    else
    {
        RVFModelActivated_ = true;
    }

    if (RVFModelActivated_ && !fixedXr_)
    {
        Info<< "Using RVF soot emission model!" << endl;
        variableXr_ = true;
    }
    else
    {
        Info<< "Using constant emission model!" << endl;
        variableXr_ = false;
    }
}


void Foam::combustionModels::eddyDissipationRVFModel::
calculateFlameTemperature()
{
    //- Get species mass fraction
    const label fuelI = this->fuelIndex();
    const volScalarField& YFuel = this->thermo().composition().Y()[fuelI];
    const volScalarField& YO2 = this->thermo().composition().Y("O2");
    const volScalarField& YN2 = this->thermo().composition().Y("N2");
    const volScalarField& YCO2 = this->thermo().composition().Y("CO2");
    const volScalarField& YH2O = this->thermo().composition().Y("H2O");
    const volScalarField& YFstar = this->thermo().composition().Y("Fstar");

    const volScalarField& TCellRef = this->thermo().T();
    const volScalarField& pCellRef = this->thermo().p();
    const dimensionedScalar qF(this->qFuel());
    const scalar s = this->s().value();
    const volScalarField rhoCell(this->thermo().rho());


    // Prepare material models
    const label O2Index = mixture_.species()["O2"];
    const label N2Index = mixture_.species()["N2"];
    const label CO2Index = mixture_.species()["CO2"];
    const label H2OIndex = mixture_.species()["H2O"];
    const word fuelName(YFuel.name());
    // fuelI
    const word& phaseName = this->thermo().phaseName();

    // Sensible enthalphy models
    const baseModels<scalar>& O2Hs =
        thermo().materials()(hsModel::typeName, phaseName, "O2");
    const baseModels<scalar>& N2Hs =
        thermo().materials()(hsModel::typeName, phaseName, "N2");
    const baseModels<scalar>& CO2Hs =
        thermo().materials()(hsModel::typeName, phaseName, "CO2");
    const baseModels<scalar>& H2OHs =
        thermo().materials()(hsModel::typeName, phaseName, "H2O");
    const baseModels<scalar>& fuelHs =
        thermo().materials()(hsModel::typeName, phaseName, fuelName);

    // Cp models
    const baseModels<scalar>& O2Cp =
        thermo().materials()(CpModel::typeName, phaseName, "O2");
    const baseModels<scalar>& N2Cp =
        thermo().materials()(CpModel::typeName, phaseName, "N2");
    const baseModels<scalar>& CO2Cp =
        thermo().materials()(CpModel::typeName, phaseName, "CO2");
    const baseModels<scalar>& H2OCp =
        thermo().materials()(CpModel::typeName, phaseName, "H2O");
    const baseModels<scalar>& fuelCp =
        thermo().materials()(CpModel::typeName, phaseName, fuelName);

    //- Get Mspecies/Mfuel from reaction equation
    scalar rCO2 =
        this->specieStoichCoeffs()
        [
            this->thermo().composition().species()["CO2"]
        ];

    scalar rH2O =
        this->specieStoichCoeffs()
        [
            this->thermo().composition().species()["H2O"]
        ];

    scalar rN2 =
        this->specieStoichCoeffs()
        [
            this->thermo().composition().species()["N2"]
        ];

    // Get Spray Info
    volScalarField sprayDensity
    (
        this->mesh().template lookupObject<volScalarField>("rhoSpray")
    );

    if (this->mesh().template foundObject<volScalarField>("rhoSprayMean_MA"))
    {
        Info<< "Using moving averaged spray density." << endl;
        sprayDensity =
            this->mesh().template lookupObject<volScalarField>
            (
                "rhoSprayMean_MA"
            );
    }

    PV_ =
        (YCO2*(1.0 + rH2O/rCO2) + SMALL)
       /(YCO2*(1.0 + rH2O/rCO2) + SMALL + min(YFuel, YO2/s)*(1.0 + s));

    forAll(YFuel, celli)
    {
	    ER_[celli] = 0.0;
	    ExpR_[celli] = 1.0;
	    ExpRMixed_[celli] = 1.0;
	    scalar pValue = pCellRef[celli];
	    scalar TValue = TCellRef[celli];

	    //- Calculate local original O2 mole fraction
        //  (convert combustion product to original species)
        scalar O2Total =
            YCO2[celli]/mixture_.W(CO2Index)
          + 0.5*YCO2[celli]*(rH2O/rCO2)/mixture_.W(H2OIndex)
          + YO2[celli]/mixture_.W(O2Index);

        XO2Local_[celli] =
            min
            (
                0.25,
                max
                (
                    O2Total
                  /(O2Total + YN2[celli]/mixture_.W(N2Index) + 1.0e-6),
                    1.0e-6
                )
            );
        scalar XO2Air(XO2Local_[celli]);
        rN2 = s*(1 - XO2Air)*mixture_.W(N2Index)/XO2Air/mixture_.W(O2Index);
	    scalar Yspray =
            max(min(Cevap_*sprayDensity[celli]/rhoCell[celli],1.0), 0.0);

	    //- Calculate adiabatic flame temperature for normal fuel
	    if ((YFuel[celli] > 1.0e-4) && (YO2[celli] > 1.0e-4))
	    {
	        scalar YN2R =
                min
                (
                    YN2[celli],
                    YO2[celli]*mixture_.W(N2Index)*(1.0 - XO2Air)
                   /(mixture_.W(O2Index)*XO2Air)
                );

	        scalar YN2P(max(0.0,YN2[celli] - YN2R));
	        scalar YFuelR(min(YFuel[celli],YO2[celli]/s));
	        scalar MFO = mixture_.W(fuelI)/mixture_.W(O2Index);
	        scalar MFN = mixture_.W(fuelI)/mixture_.W(N2Index);

	        ER_[celli] =
                YFuelR*(1.0 + s*MFO + s*MFN*YN2R/YO2[celli])
	          /(YFuel[celli] + MFO*YO2[celli] + MFN*YN2[celli]);
	        scalar ERValue = ER_[celli];
	        scalar CoffN2 = ERValue*YN2P + rN2*YFuelR;
	        scalar CoffCO2 = ERValue*YCO2[celli] + rCO2*YFuelR;
	        scalar CoffH2O = ERValue*(YH2O[celli] + Yspray) + rH2O*YFuelR;
	        scalar CoffFstar = ERValue*YFstar[celli];
	        scalar CoffO2 = 0;
	        scalar RHS1 = YFuelR*qF.value()*(1.0 - XrExt_);
	        scalar RHS2 = min(RHS1, ERValue*Yspray*3.0e6);
	        scalar RHS3 =
                ERValue
               *(
                    YO2[celli]*O2Hs.value(pValue, TValue)
	    	      + YN2[celli]*N2Hs.value(pValue, TValue)
	    	      + YCO2[celli]*CO2Hs.value(pValue, TValue)
	    	      + YH2O[celli]*H2OHs.value(pValue, TValue)
	    	      + (YFuel[celli] + YFstar[celli])*fuelHs.value(pValue,TValue)
	    	    );

	        scalar RHS4 =
	    	    (YFuelR - ERValue*YFuel[celli])*fuelHs.value(pValue, TFuel_)
	    	  + (s*YFuelR - ERValue*YO2[celli])*O2Hs.value(pValue, TAir_)
	    	  + (rN2*YFuelR - ERValue*YN2R)*N2Hs.value(pValue, TAir_);

	        scalar RHSall = RHS1 - RHS2 + RHS3 + RHS4;

	        scalar Test = TValue;
	        scalar Tnew = TValue;
	        int iter = 0;
	        do
	        {
	            Test = Tnew;
	            scalar CpEff =
	                CoffN2*N2Cp.value(pValue, Test)
	              + CoffCO2*CO2Cp.value(pValue, Test)
	              + CoffH2O*H2OCp.value(pValue, Test)
	              + CoffFstar*fuelCp.value(pValue, Test)
	              + CoffO2*O2Cp.value(pValue, Test);

	            scalar LHSest =
	                CoffN2*N2Hs.value(pValue, Test)
	              + CoffCO2*CO2Hs.value(pValue, Test)
	              + CoffH2O*H2OHs.value(pValue, Test)
	              + CoffFstar*fuelHs.value(pValue, Test)
	              + CoffO2*O2Hs.value(pValue, Test);

	            Tnew = Test + (RHSall - LHSest)/CpEff;
	            if (iter++ > 10)
	            {
	                FatalErrorInFunction
	                    << "Maximum number of iteration exceeded in Tad model"
                        << Tnew << nl << "Local XO2: " << XO2Air << tab
                        << celli << tab<< iter << nl << "Coeffs: " << CoffN2
                        << tab << CoffH2O << tab << CoffCO2 << tab << CoffO2
                        << nl << "RHSall: " << RHSall << tab << LHSest << tab
                        << CpEff << nl
	                    << abort(FatalError);
	            }
	            if ((Tnew > 5000) || (Tnew < 200))
	            {
	                FatalErrorInFunction
	                    << "Tad exceed range of [200, 5000]: " << Tnew << nl
	                    << "Local XO2: " << XO2Air << tab << celli << tab
                        << iter << nl << "O2 (total 1) & N2: " << O2Total
                        << tab << YN2[celli] << nl << "Coeffs: " << CoffN2
                        << tab << CoffH2O << tab << CoffCO2 << tab
                        << CoffO2 << nl << "RHSall: " << RHSall << tab
                        << LHSest << tab << CpEff << nl << "Yspray: "
                        << Yspray << nl << "ERValue: " << ERValue << nl
	                    << abort(FatalError);
	            }
	        } while (mag(Tnew - Test) > 5.0);

	        Tad_[celli] = Tnew;
	        ExpR_[celli] =
                (Tnew/TAir_)
	           *(
                    (
                        rCO2/mixture_.W(CO2Index)
                      + rH2O/mixture_.W(H2OIndex)
                      + rN2/mixture_.W(N2Index)
                    )
                   /(
                        1.0/mixture_.W(fuelI)
                      + s/mixture_.W(O2Index)
                      + rN2/mixture_.W(N2Index)
                    )
	            );
	    }
	    else
	    {
	        Tad_[celli] = TValue;
	    }

	    //- Calculate adiabatic flame temperature for mixed fuel
	    if ((YFstar[celli] > 1.0e-4) && (YO2[celli] > 1.0e-4))
	    {
	        scalar YFuelR = min(YFstar[celli], YO2[celli]/s);
	        scalar CoffN2 = YN2[celli];
	        scalar CoffH2O = YH2O[celli] + rH2O*YFuelR + Yspray;
	        scalar CoffCO2 = YCO2[celli] + rCO2*YFuelR;
	        scalar CoffFstar = YFstar[celli] - YFuelR;
	        scalar CoffO2 = YO2[celli] - s*YFuelR;

	        scalar RHS1 = YFuelR*qF.value()*(1.0 - XrExt_);
	        scalar RHS2 = min(RHS1, Yspray*3.0e6);
	        scalar RHS3 =
                YO2[celli]*O2Hs.value(pValue, TValue)
	          + YN2[celli]*N2Hs.value(pValue, TValue)
	          + YCO2[celli]*CO2Hs.value(pValue, TValue)
	          + YH2O[celli]*H2OHs.value(pValue, TValue)
	          + YFstar[celli]*fuelHs.value(pValue, TValue);

	        scalar RHS4 = 0;

	        scalar RHSall = RHS1 - RHS2 + RHS3 + RHS4;

	        scalar Test = TValue;
	        scalar Tnew = TValue;
	        int iter(0);
	        do
	        {
	    	    Test = Tnew;
	    	    scalar CpEff =
	    	        CoffN2*N2Cp.value(pValue, Test)
	    	      + CoffCO2*CO2Cp.value(pValue, Test)
	    	      + CoffH2O*H2OCp.value(pValue, Test)
	    	      + CoffFstar*fuelCp.value(pValue, Test)
	    	      + CoffO2*O2Cp.value(pValue, Test);

	    	    scalar LHSest =
	    	        CoffN2*N2Hs.value(pValue, Test)
	    	      + CoffCO2*CO2Hs.value(pValue, Test)
	    	      + CoffH2O*H2OHs.value(pValue, Test)
	    	      + CoffFstar*fuelHs.value(pValue, Test)
	    	      + CoffO2*O2Hs.value(pValue, Test);

	    	    Tnew = Test + (RHSall - LHSest)/CpEff;
	    	    if (iter++ > 10)
	    	    {
	    	        FatalErrorInFunction
	    	            << "Maximum number of iteration exceeded in Tad model"
	    	   	        << abort(FatalError);
	    	    }
	    	    if ((Tnew > 5000) || (Tnew < 200))
	    	    {
	    	        FatalErrorInFunction
	    	            << "Tad exceed range of [200, 5000]: " << Tnew << nl
	    	            << "Local XO2: " << XO2Air << tab << celli << tab
                        << iter << nl << "O2 (total 2) & N2: " << O2Total
                        << tab << YN2[celli] << nl << "Coeffs: " << CoffN2
                        << tab << CoffH2O << tab << CoffCO2 << tab << CoffO2
                        << nl << "RHSall: " << RHSall << tab << LHSest << tab
                        << CpEff << nl << "Yspray: " << Yspray << nl
	    	            << abort(FatalError);
	    	    }
	        } while (mag(Tnew - Test) > 5.0);

	        TadMixed_[celli] = Tnew;
	        ExpRMixed_[celli] =
                (Tnew/TAir_)
	    	   *(
                    (
                        rCO2/mixture_.W(CO2Index)
                      + rH2O/mixture_.W(H2OIndex)
                      + rN2/mixture_.W(N2Index)
                    )
                   /(
                        1.0/mixture_.W(fuelI)
                      + s/mixture_.W(O2Index)
                      + rN2/mixture_.W(N2Index)
                    )
	    	    );
	    }
	    else
	    {
	        TadMixed_[celli] = TValue;
	    }

        //- Local soot volume fraction near flame sheet
        scalar XO2Eff(0);
        if (YFuel[celli] > 1.0e-4)
        {
            XO2Eff = max(0, 0.21 - (TadAir_ - Tad_[celli])/8055);
        }
        else if (YFstar[celli] > 1.0e4)
        {
            XO2Eff = max(0, 0.21 - (TadAir_ - TadMixed_[celli])/8055);
        }
        else
        {
            XO2Eff = 0;
        }
        fvSoot_[celli] = max(0, fvSootAir_*(XO2Eff - O2Soot_)/(0.21-O2Soot_));
    }
}


Foam::tmp<Foam::fvScalarMatrix>
Foam::combustionModels::eddyDissipationRVFModel::R(volScalarField& Y) const
{
    const label specieI = this->thermo().composition().species()[Y.name()];
    const label fuelI = this->fuelIndex();
    if (specieI == fuelI)
    {
        volScalarField wSpecie
        (
            this->wFuel_*this->specieStoichCoeffs()[specieI]
        );
        return wSpecie + fvm::Sp(0.0*wSpecie, Y);
    }
    else if (Y.name() == "Fstar")
    {
        volScalarField wSpecie
        (
           (1-RVF_)*this->wFuel_ - Fig_*this->WFstar_
        );
        return wSpecie + fvm::Sp(0.0*wSpecie, Y);
    }
    else
    {
        volScalarField wSpecie
        (
           (RVF_*this->wFuel_ + Fig_*this->WFstar_)
	       *this->specieStoichCoeffs()[specieI]
        );
        return wSpecie + fvm::Sp(0.0*wSpecie, Y);
    }
}


Foam::tmp<Foam::volScalarField>
Foam::combustionModels::eddyDissipationRVFModel::Qdot() const
{
    const label fuelI = this->fuelIndex();
    volScalarField& YFuel =
        const_cast<volScalarField&>(this->thermo().composition().Y(fuelI));

    const label indexFstar = this->thermo().composition().species()["Fstar"];
    volScalarField& YFstar =
        const_cast<volScalarField&>(this->thermo().composition().Y(indexFstar));

    return -this->qFuel()*((R(YFuel) & YFuel) + (R(YFstar) & YFstar));
}


Foam::tmp<Foam::volScalarField>
Foam::combustionModels::eddyDissipationRVFModel::QdotFuel() const
{
    const label fuelI = this->fuelIndex();
    volScalarField& YFuel =
        const_cast<volScalarField&>(this->thermo().composition().Y(fuelI));

    return -this->qFuel()*(R(YFuel) & YFuel)*RVF_;
}


Foam::tmp<Foam::volScalarField>
Foam::combustionModels::eddyDissipationRVFModel::QdotFstar() const
{
    const label fuelI = this->fuelIndex();
    volScalarField& YFuel =
        const_cast<volScalarField&>(this->thermo().composition().Y(fuelI));

    const label indexFstar(this->thermo().composition().species()["Fstar"]);
    volScalarField& YFstar =
        const_cast<volScalarField&>
        (
            this->thermo().composition().Y(indexFstar)
        );

    return
        -this->qFuel()
        *((R(YFuel) & YFuel)*(1-RVF_)+(R(YFstar) & YFstar));
}


bool Foam::combustionModels::eddyDissipationRVFModel::variableXr() const
{
    return variableXr_;
}


Foam::scalar
Foam::combustionModels::eddyDissipationRVFModel::RVFactivationTime() const
{
    return tExt_;
}


bool Foam::combustionModels::eddyDissipationRVFModel::read()
{
    if (singleStepCombustion::read())
    {
        C_ = this->coeffs().template lookup<scalar>("C");
        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
