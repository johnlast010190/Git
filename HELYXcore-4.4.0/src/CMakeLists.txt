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
    (c) 2019-2024 Engys Ltd.

[----------------------------------------------------------------------------]]

if(${TBB_FOUND})
    add_compile_definitions(HELYX_USE_TBB)
endif()

add_subdirectory(OSspecific)

add_subdirectory(OpenFOAM)

# All libraries from this point on should link OpenFOAM
link_libraries(OpenFOAM)

add_subdirectory(ODE)

add_subdirectory(TurbulenceModels)

add_subdirectory(combustionModels)

add_subdirectory(conversion)

add_subdirectory(dynamicFvMesh)

add_subdirectory(deformingBodyMotionFunctions)

add_subdirectory(dynamicMesh)

add_subdirectory(engine)

add_subdirectory(fileFormats)

add_subdirectory(finiteArea)

add_subdirectory(finiteVolume)

add_subdirectory(functionObjects)

add_subdirectory(fvAgglomerationMethods)

add_subdirectory(fvMeshMovers)

add_subdirectory(fvMeshTopoChangers)

add_subdirectory(fvMeshDistributors)

add_subdirectory(fvMeshGIBChangers)

add_subdirectory(fvMeshStitchers)

add_subdirectory(fvMotionSolver)

add_subdirectory(fvOptions)

add_subdirectory(genericPatchFields)

add_subdirectory(helyxMap)

add_subdirectory(helyxSolve)

add_subdirectory(lagrangian)

add_subdirectory(mesh)

add_subdirectory(meshTools)

add_subdirectory(overset)

add_subdirectory(parallel)

add_subdirectory(radiationModels)

if(${FFTW_FOUND})
    add_subdirectory(randomProcesses)
else()
    add_dummy_target(randomProcesses FFTW)
endif()

add_subdirectory(regionModels)

add_subdirectory(renumber)

add_subdirectory(rigidBodyDynamics)

add_subdirectory(rigidBodyMeshMotion)

add_subdirectory(sampling)

add_subdirectory(sixDoFRigidBodyMotion)

add_subdirectory(surfMesh)

add_subdirectory(thermophysicalModels)

add_subdirectory(materialModels)

add_subdirectory(topoChangerFvMesh)

add_subdirectory(transportModels)

add_subdirectory(waveModels)
