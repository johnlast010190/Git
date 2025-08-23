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
    (c) 2021-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "materialModels.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    declareMaterialModelScalarDims(W, dimMass/dimMoles);
    declareMaterialModelScalarDims(R, dimEnergy/dimMass/dimTemperature);
    declareMaterialModelScalarDims(Y, dimless);
    declareMaterialModelScalarDims(limit, dimless);
    declareMaterialModelScalarDims(rho, dimMass/dimVolume);
    declareMaterialModelScalarDims(hContribution, dimEnergy/dimMass);
    declareMaterialModelScalarDims(CpContribution, dimEnergy/dimMass/dimTemperature);
    declareMaterialModelScalarDims(eContribution, dimEnergy/dimMass);
    declareMaterialModelScalarDims(CvContribution, dimEnergy/dimMass/dimTemperature);
    declareMaterialModelScalarDims(spContribution, dimEnergy/dimMass/dimTemperature);
    declareMaterialModelScalarDims(svContribution, dimEnergy/dimMass/dimTemperature);
    declareMaterialModelScalarDims(psi, sqr(dimTime)/dimArea);
    declareMaterialModelScalarDims(Z, dimless);
    declareMaterialModelScalarDims(CpMCv, dimEnergy/dimMass/dimTemperature);
    declareMaterialModelScalarDims(hf, dimEnergy/dimMass);
    declareMaterialModelScalarDims(a, dimEnergy/dimMass);
    declareMaterialModelScalarDims(Cp, dimEnergy/dimMass/dimTemperature);
    declareMaterialModelScalarDims(Cpv, dimEnergy/dimMass/dimTemperature);
    declareMaterialModelScalarDims(Cv, dimEnergy/dimMass/dimTemperature);
    declareMaterialModelScalarDims(ea, dimEnergy/dimMass);
    declareMaterialModelScalarDims(es, dimEnergy/dimMass);
    declareMaterialModelScalarDims(gamma, dimless);
    declareMaterialModelScalarDims(g, dimEnergy/dimMass);
    declareMaterialModelScalarDims(ha,dimEnergy/dimMass);
    declareMaterialModelScalarDims(he, dimEnergy/dimMass);
    declareMaterialModelScalarDims(hs, dimEnergy/dimMass);
    declareMaterialModelScalarDims(s, dimEnergy/dimMass/dimTemperature);
    declareMaterialModelScalarDims(strainRate, dimless/dimTime);
    declareMaterialModelScalarDims(mu, dimDynamicViscosity);
    declareMaterialModelScalarDims(kappa, dimPower/(dimLength*dimTemperature));
    declareMaterialModelVectorDims(vKappa, dimPower/(dimLength*dimTemperature));

    // Mixture models
    declareMaterialModelScalarMixture(weightedAverageMixture, dimless);
    declareMaterialModelVectorMixture(vectorWeightedAverageMixture, dimless);
    declareMaterialModelTensorMixture(tensorWeightedAverageMixture, dimless);
    declareMaterialModelScalarMixture(clippedWeightedAverageMixture, dimless);
    declareMaterialModelVectorMixture(vectorClippedWeightedAverageMixture, dimless);
    declareMaterialModelTensorMixture(tensorClippedWeightedAverageMixture, dimless);
    declareMaterialModelScalarMixture(harmonicMixture, dimless);

    // Combustion mixtures:
    declareMaterialModelScalarMixture(egrReactantMixture, dimless);
    declareMaterialModelScalarMixture(egrProductMixture, dimless);
    declareMaterialModelScalarMixture(egrThermoMixture, dimless);
    declareMaterialModelScalarMixture(homogeneousReactantMixture, dimless);
    declareMaterialModelScalarMixture(homogeneousProductMixture, dimless);
    declareMaterialModelScalarMixture(homogeneousThermoMixture, dimless);
    declareMaterialModelScalarMixture(inhomogeneousReactantMixture, dimless);
    declareMaterialModelScalarMixture(inhomogeneousProductMixture, dimless);
    declareMaterialModelScalarMixture(inhomogeneousThermoMixture, dimless);
    declareMaterialModelScalarMixture(veryInhomogeneousReactantMixture, dimless);
    declareMaterialModelScalarMixture(veryInhomogeneousProductMixture, dimless);
    declareMaterialModelScalarMixture(veryInhomogeneousThermoMixture, dimless);

    declareMaterialModelScalarDims(T, dimTemperature);
    declareMaterialModelScalarDims(Tha, dimTemperature);
    declareMaterialModelScalarDims(Ths, dimTemperature);
    declareMaterialModelScalarDims(Tea, dimTemperature);
    declareMaterialModelScalarDims(Tes, dimTemperature);
    declareMaterialModelScalarDims(p, dimPressure);
    declareMaterialModelScalarDims(dCpdT, dimEnergy/dimMass);
    declareMaterialModelScalarDims(alphav, dimless/dimTemperature);

    // Liquid properties
    declareMaterialModelScalarDims(pv, dimPressure);
    declareMaterialModelScalarDims(hl, dimEnergy/dimMass);
    declareMaterialModelScalarDims(Cpg, dimEnergy/dimMass/dimTemperature);
    declareMaterialModelScalarDims(B, dimVolume/dimMass);
    declareMaterialModelScalarDims(kappag, dimEnergy/(dimLength*dimTemperature));
    declareMaterialModelScalarDims(mug, dimPressure*dimTime);
    declareMaterialModelScalarDims(sigma, dimForce/dimLength);
    declareMaterialModelScalarDims(D, dimArea/dimTime);
    declareMaterialModelScalarDims(pvInvert, dimTemperature);
    declareMaterialModelScalarDims(c0, dimLength/dimTime);

    // Sturation models
    declareMaterialModelScalarDims(pSat, dimPressure);
    declareMaterialModelScalarDims(pSatPrime, dimPressure/dimTemperature);
    declareMaterialModelScalarDims(lnPSat, dimless);
    declareMaterialModelScalarDims(Tsat, dimTemperature);

    // Solid properties
    declareMaterialModelScalarDims(E, dimPressure);
    declareMaterialModelScalarDims(nu, dimless);
}


// ************************************************************************* //
