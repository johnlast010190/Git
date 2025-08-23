#[[---------------------------------------------------------------------------]
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.4.0
|    o     o     |  ENGYS Ltd. <http://engys.com/>
|       o        |
[-----------------------------------------------------------------------------]
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
    (c) 2022-2024 Engys Ltd.

[----------------------------------------------------------------------------]]

set(TARGET_NAME "materialModels")

add_helyx_library(${TARGET_NAME} SHARED

# Selection and collection of models
materialModels/baseModels/materialModels.C
materialModels/baseModels/generalModel.C
materialModels/materialModel/materialModel.C
materialModels/materialModel/materialModelNew.C
materialModels/materialTables/materialTables.C

# Species models
materialsSpecie/specieModels/specieModels.C
materialsSpecie/materialsSpecie.C

# Generic models
general/scalarModelFunction1/scalarModelFunction1.C
general/tabulated/tabulated.C
general/logPolynomial/logPolynomial.C
general/fieldLimiter/fieldLimiter.C

# Water-steam model
general/IAPWSIF97EOS/IAPWSIF97EOS.C
general/IAPWSIF97EOS/region1IAPWSIF97EOS.C
general/IAPWSIF97EOS/region2IAPWSIF97EOS.C
general/IAPWSIF97EOS/region3IAPWSIF97EOS.C
general/IAPWSIF97EOS/region5IAPWSIF97EOS.C
general/IAPWSIF97EOS/IAPWSIF97EOSModels.C

# Equation of state models
equationOfState/equationOfState/equationOfState.C
equationOfState/perfectGasLaw/perfectGasLaw.C
equationOfState/hybridConstPerfectGas/hybridConstPerfectGas.C
equationOfState/adiabaticPerfectFluidLaw/adiabaticPerfectFluidLaw.C
equationOfState/BoussinesqLaw/BoussinesqLaw.C
equationOfState/icoPolynomialLaw/icoPolynomialLaw.C
equationOfState/incompressiblePerfectGasLaw/incompressiblePerfectGasLaw.C
equationOfState/linearLaw/linearLaw.C
equationOfState/PengRobinsonGasLaw/PengRobinsonGasLaw.C
equationOfState/SoaveRedlichKwongGas/SoaveRedlichKwongGas.C
equationOfState/RedlichKwongPengRobinsonGas/RedlichKwongPengRobinsonGas.C
equationOfState/rhoConstLaw/rhoConstLaw.C
equationOfState/rhoFunction1/rhoFunction1.C
equationOfState/rPolynomialLaw/rPolynomialLaw.C

# Thermo related models
thermo/thermo/thermodynamics/thermodynamics.C
thermo/thermo/hConst/hConst.C
thermo/thermo/eConst/eConst.C
thermo/thermo/ePower/ePower.C
thermo/thermo/ePolynomial/ePolynomial.C
thermo/thermo/hFunctionT/hFunctionT.C
thermo/thermo/hPolynomial/hPolynomial.C
thermo/thermo/hPower/hPower.C
thermo/thermo/hTabulated/hTabulated.C
thermo/thermo/janaf/janaf.C

thermo/energyConversion/energyConversion/energyConversion.C
thermo/energyConversion/standardThermo/standardThermo.C

thermo/energy/energy/energy.C
thermo/energy/matSensibleInternalEnergy/matSensibleInternalEnergy.C
thermo/energy/matSensibleEnthalpy/matSensibleEnthalpy.C
thermo/energy/matAbsoluteEnthalpy/matAbsoluteEnthalpy.C
thermo/energy/matAbsoluteInternalEnergy/matAbsoluteInternalEnergy.C
thermo/energy/energyInverseTemperature/energyInverseTemperature.C

# Transport models
transport/muSutherland/muSutherland.C
transport/kappaSutherland/kappaSutherland.C
transport/kappaPrandtl/kappaPrandtl.C

# Solid transport models
transport/constSolid/constSolid.C
transport/exponentialSolid/exponentialSolid.C
transport/functionTIso/functionTIso.C
transport/functionTAnIso/functionTAnIso.C
transport/solidPolynomial/solidPolynomial.C

# Transport mu models
transport/strainRateLaw/strainRateLaw.C
transport/approachingMuLaw/approachingMuLaw.C
transport/BirdCarreauMu/BirdCarreauMu.C
transport/BokilBewtraMu/BokilBewtraMu.C
transport/CassonMu/CassonMu.C
transport/CrossPowerLawMu/CrossPowerLawMu.C
transport/HerschelBulkleyMu/HerschelBulkleyMu.C
transport/powerLawMu/powerLawMu.C
transport/strainRateFunctionMu/strainRateFunctionMu.C
transport/WLFMu/WLFMu.C
transport/vapourDiffusivity/vapourDiffusivity.C
transport/pvInvertFunc/pvInvertFunc.C
transport/speedSoundGas/speedSoundGas.C

# Mixture models
mixtures/weightedAverageMixture/weightedAverageMixture.C
mixtures/clippedWeightedAverageMixture/clippedWeightedAverageMixture.C
mixtures/harmonicMixture/harmonicMixture.C
mixtures/coefficientMixture/coefficientMixture.C
mixtures/coefficientMixture/coefficientMixtureNew.C
mixtures/egrMixture/egrThermoMixture.C
mixtures/homogeneousMixture/homogeneousThermoMixture.C
mixtures/veryInhomogeneousMixture/veryInhomogeneousThermoMixture.C
mixtures/inhomogeneousMixture/inhomogeneousThermoMixture.C

# Saturation models
saturationModels/Antoine/Antoine.C
saturationModels/AntoineExtended/AntoineExtended.C
saturationModels/ArdenBuck/ArdenBuck.C
saturationModels/constantPressureSaturation/constantPressureSaturation.C
)

helyx_link_libraries(${TARGET_NAME}
    PRIVATE finiteVolume
    PUBLIC thermophysicalMixtures
)
