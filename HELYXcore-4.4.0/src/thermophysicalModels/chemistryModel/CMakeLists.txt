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
    (c) 2019-2025 Engys Ltd.

[----------------------------------------------------------------------------]]

set(TARGET_NAME "chemistryModel")

add_helyx_library(${TARGET_NAME} SHARED
./basicChemistryModel/basicChemistryModel.C
./basicChemistryModel/basicChemistryModelNew.C

./chemistrySolver/chemistrySolver/chemistrySolvers.C
./chemistrySolver/noChemistrySolver/noChemistrySolvers.C
./chemistrySolver/EulerImplicit/EulerImplicitChemistrySolvers.C
./chemistrySolver/ode/odeChemistrySolvers.C

./odeChemistryModel/odeChemistryModel.C

./chemistryModel/reduction/chemistryReductionMethod/chemistryReductionMethods.C
./chemistryModel/reduction/noChemistryReduction/noChemistryReductionMethods.C
./chemistryModel/reduction/DAC/DACChemistryReductionMethods.C
./chemistryModel/reduction/DRG/DRGChemistryReductionMethods.C
./chemistryModel/reduction/DRGEP/DRGEPChemistryReductionMethods.C
./chemistryModel/reduction/EFA/EFAChemistryReductionMethods.C
./chemistryModel/reduction/PFA/PFAChemistryReductionMethods.C

./chemistryModel/tabulation/chemistryTabulationMethod/chemistryTabulationMethod.C
./chemistryModel/tabulation/chemistryTabulationMethod/chemistryTabulationMethodNew.C
./chemistryModel/tabulation/noChemistryTabulation/noChemistryTabulation.C
./chemistryModel/tabulation/ISAT/ISAT.C
./chemistryModel/tabulation/ISAT/chemPointISAT/chemPointISAT.C
./chemistryModel/tabulation/ISAT/binaryNode/binaryNode.C
./chemistryModel/tabulation/ISAT/binaryTree/binaryTree.C

./reaction/makeMaterialReactions.C
./reaction/makeGasesReactions.C
./reaction/makeLiquidReactions.C
./reaction/makePolynomialReactions.C
./reaction/makeLangmuirHinshelwoodReactions.C
./reaction/makeMichaelisMentenReactions.C
./reaction/makeSolidReactions.C
./reaction/makefluxLimitedLangmuirHinshelwoodReactions.C
./reaction/makesurfaceArrheniusReactions.C
./functionObjects/specieReactionRates/specieReactionRates.C

./basicSolidChemistryModel/basicSolidChemistryModel.C
./basicSolidChemistryModel/basicSolidChemistryModelNew.C
./basicSolidChemistryModel/basicSolidChemistryModels.C
)

helyx_link_libraries(${TARGET_NAME}
    PUBLIC fluidThermophysicalModels
    PRIVATE materialModels
    PRIVATE multiphaseMaterials
    PRIVATE specie
    PRIVATE solidThermo
    PRIVATE ODE
    PRIVATE finiteVolume
    PRIVATE meshTools
)


set_source_files_properties(
    ${CMAKE_CURRENT_LIST_DIR}/chemistryModel/tabulation/chemistryTabulationMethod/chemistryTabulationMethods.C
    ${CMAKE_CURRENT_LIST_DIR}/chemistryModel/tabulation/noChemistryTabulation/noChemistryTabulationMethods.C
    ${CMAKE_CURRENT_LIST_DIR}/chemistryModel/tabulation/ISAT/ISATChemistryTabulationMethods.C
    ${CMAKE_CURRENT_LIST_DIR}/chemistryModel/basicSolidChemistryModel/basicSolidChemistryModels.C
    PROPERTIES SKIP_UNITY_BUILD_INCLUSION ON
)
