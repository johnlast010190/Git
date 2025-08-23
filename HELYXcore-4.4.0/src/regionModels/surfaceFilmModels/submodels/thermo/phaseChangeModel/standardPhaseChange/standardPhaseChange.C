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

#include "surfaceFilmModels/submodels/thermo/phaseChangeModel/standardPhaseChange/standardPhaseChange.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "surfaceFilmModels/thermoSingleLayer/thermoSingleLayer.H"
#include "mixtures/basicSpecieMixture/basicSpecieMixture.H"
#include "fields/Fields/zeroField/zeroField.H"
#include "materialModels/materialTables/materialTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(standardPhaseChange, 0);
addToRunTimeSelectionTable(phaseChangeModel, standardPhaseChange, dictionary);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

scalar standardPhaseChange::Sh
(
    const scalar Re,
    const scalar Sc
) const
{
    if (Re < 5.0E+05)
    {
        return 0.664*sqrt(Re)*cbrt(Sc);
    }
    else
    {
        return 0.037*pow(Re, 0.8)*cbrt(Sc);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

standardPhaseChange::standardPhaseChange
(
    surfaceFilmRegionModel& film,
    const dictionary& dict
)
:
    speciePhaseChange(typeName, film, dict),
    deltaMin_(coeffDict_.lookup<scalar>("deltaMin")),
    L_(coeffDict_.lookup<scalar>("L")),
    TbFactor_(coeffDict_.lookupOrDefault<scalar>("TbFactor", 1.1)),
    YInfZero_(coeffDict_.lookupOrDefault<Switch>("YInfZero", false))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

standardPhaseChange::~standardPhaseChange()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class YInfType>
void standardPhaseChange::correctModel
(
    const scalar dt,
    scalarField& availableMass,
    scalarField& dMass,
    scalarField& dEnergy,
    YInfType YInf
)
{
    const thermoSingleLayer& film = filmType<thermoSingleLayer>();

    // Set local liquidThermo properties
    const rhoThermo& liquidThermo = film.thermo();
    const materialTables& liquidTables = liquidThermo.materials();

    const basicSpecieMixture& primarySpecieThermo =
        refCast<const basicSpecieMixture>(film.primaryThermo());

    // Retrieve fields from film model
    const scalarField& delta = film.delta();
    const scalarField T(film.thermo().TAbs());
    const scalarField pInf(film.thermo().pAbs());
    const scalarField& rho = film.rho();
    const scalarField& rhoInf = film.rhoPrimary();
    const scalarField& muInf = film.muPrimary();
    const scalarField& magSf = film.magSf();
    const vectorField dU(film.UPrimary() - film.Us());
    const scalarField limMass
    (
        max(scalar(0), availableMass - deltaMin_*rho*magSf)
    );

    // molecular weight of vapour [kg/kmol]
    const scalar Wvap = this->Wvap();

    // molecular weight of liquid [kg/kmol]
    const scalar Wliq =
        liquidTables(WModel::typeName, film.thermo().phaseName(), word::null)[0];

    forAll(dMass, celli)
    {
        scalar dm = 0;

        if (delta[celli] > deltaMin_)
        {
            // Cell pressure [Pa]
            const scalar pc = pInf[celli];

            // Calculate the boiling temperature
            const scalar Tb = liquidTables("pvInvert")[celli];

            // Local temperature - impose lower limit of 200 K for stability
            const scalar Tloc = min(TbFactor_*Tb, max(200.0, T[celli]));

            const scalar Tback = T[celli];

            const_cast<volScalarField&>(liquidThermo.T())[celli] =
                (Tloc - liquidThermo.TRefValue());

            // Saturation pressure [Pa]
            const scalar pSat = liquidTables(pvModel::typeName)[celli];

            // Latent heat [J/kg]
            const scalar hVap = liquidTables(hlModel::typeName)[celli];

            // Calculate mass transfer
            if (pSat >= 0.95*pc)
            {
                // Boiling
                const scalar Cp = liquidTables(CpModel::typeName)[celli];
                const scalar Tcorr = max(0.0, T[celli] - Tb);
                const scalar qCorr = limMass[celli]*Cp*(Tcorr);
                dm = qCorr/hVap;
            }
            else
            {
                // Primary region density [kg/m3]
                const scalar rhoInfc = rhoInf[celli];

                // Primary region viscosity [Pa.s]
                const scalar muInfc = muInf[celli];

                // Reynolds number
                const scalar Re = rhoInfc*mag(dU[celli])*L_/muInfc;

                // vapour mass fraction at interface
                const scalar Ys = Wliq*pSat/(Wliq*pSat + Wvap*(pc - pSat));

                // vapour diffusivity [m2/s]
                const scalar Dab = liquidTables(DModel::typeName)[celli];

                // Schmidt number
                const scalar Sc = muInfc/(rhoInfc*(Dab + ROOTVSMALL));

                // Sherwood number
                const scalar Sh = this->Sh(Re, Sc);

                // mass transfer coefficient [m/s]
                const scalar hm = Sh*Dab/(L_ + ROOTVSMALL);

                // add mass contribution to source
                dm = dt*magSf[celli]*rhoInfc*hm*(Ys - YInf[celli])/(1.0 - Ys);
            }

            dm = min(limMass[celli], max(dm, 0));

            dMass[celli] += dm;

            // Assume that the vapour transferred to the primary region is
            // already at temperature Tloc so that all heat required for
            // the phase-change is provided by the film

            // Needs to be tested!!!
            volScalarField& pField =
                const_cast<volScalarField&>(film.primaryThermo().p());
            volScalarField& TField =
                const_cast<volScalarField&>(film.primaryThermo().T());
            const scalar pBack = pField[0];
            const scalar TBack = TField[0];
            pField[0] = (pc - film.primaryThermo().pRefValue());
            TField[0] = (Tloc - film.primaryThermo().TRefValue());
            const scalar hs =
                film.primaryThermo().materials()
                (
                    hsModel::typeName,
                    film.primaryThermo().phaseName(),
                    primarySpecieThermo.Y(vapId()).name()
                )[0];
            pField[0] = pBack;
            TField[0] = TBack;
            dEnergy[celli] += dm*hs;
            const_cast<volScalarField&>(liquidThermo.T())[celli] = Tback;
        }
    }
}


void standardPhaseChange::correctModel
(
    const scalar dt,
    scalarField& availableMass,
    scalarField& dMass,
    scalarField& dEnergy
)
{
    if (YInfZero_)
    {
        correctModel(dt, availableMass, dMass, dEnergy, zeroField());
    }
    else
    {
        const thermoSingleLayer& film = filmType<thermoSingleLayer>();
        const scalarField& YInf = film.YPrimary()[vapId()];

        correctModel(dt, availableMass, dMass, dEnergy, YInf);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
