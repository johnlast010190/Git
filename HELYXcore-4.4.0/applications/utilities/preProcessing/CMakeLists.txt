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
    (c) 2020-2025 Engys Ltd.

[----------------------------------------------------------------------------]]

add_subdirectory(applyBoundaryLayer)

if(${FFTW_FOUND})
    add_subdirectory(boxTurb)
else()
    add_dummy_target(boxTurb-exe FFTW)
endif()

add_subdirectory(caseSetup)

add_subdirectory(changeDictionary)

add_subdirectory(compareFields)

add_subdirectory(createCase)

add_subdirectory(createExternalCoupledPatchGeometry)

add_subdirectory(createZeroDirectory)

add_subdirectory(dsmcInitialise)

add_subdirectory(engineSwirl)

add_subdirectory(exportObservers)

add_subdirectory(faceAgglomerate)

add_subdirectory(frontalArea)

add_subdirectory(helyxMeshSetup)

add_subdirectory(mapFields)

add_subdirectory(mapFieldsPar)

add_subdirectory(mdInitialise)

add_subdirectory(meshUpdate)

add_subdirectory(modifyGIBBoundary)

add_subdirectory(perturbU)

add_subdirectory(scalePressure)

add_subdirectory(scaleUprime)

add_subdirectory(setAlphaField)

add_subdirectory(setFields)

add_subdirectory(viewFactorsGen)

add_subdirectory(helyxMap)

#if (${OPENVDB_FOUND})
#    add_subdirectory(helyxVDBMap)
#else ()
#    add_dummy_target(helyxVDBMap OpenVDB)
#endif ()



