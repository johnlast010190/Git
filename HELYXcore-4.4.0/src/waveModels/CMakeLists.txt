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
    (c) 2019 Engys Ltd.

[----------------------------------------------------------------------------]]

set(TARGET_NAME "waveModels")

add_helyx_library(${TARGET_NAME} SHARED
./waveModel/waveModel.C
./waveModel/waveModelNew.C
./waveGenerationModels/base/waveGenerationModel/waveGenerationModel.C
./waveGenerationModels/base/regularWaveModel/regularWaveModel.C
./waveGenerationModels/base/solitaryWaveModel/solitaryWaveModel.C
./waveGenerationModels/derived/Boussinesq/BoussinesqWaveModel.C
./waveGenerationModels/derived/cnoidal/cnoidalWaveModel.C
./waveGenerationModels/derived/Grimshaw/GrimshawWaveModel.C
./waveGenerationModels/derived/McCowan/McCowanWaveModel.C
./waveGenerationModels/derived/StokesII/StokesIIWaveModel.C
./waveGenerationModels/derived/StokesI/StokesIWaveModel.C
./waveGenerationModels/derived/StokesV/StokesVWaveModel.C
./waveAbsorptionModels/base/waveAbsorptionModel/waveAbsorptionModel.C
./waveAbsorptionModels/derived/shallowWaterAbsorption/shallowWaterAbsorption.C
./derivedFvPatchFields/waveVelocity/waveVelocityFvPatchVectorField.C
./derivedFvPatchFields/waveAlpha/waveAlphaFvPatchScalarField.C
)



helyx_link_libraries(${TARGET_NAME}
    PRIVATE finiteVolume
)
