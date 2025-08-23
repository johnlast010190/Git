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

set(TARGET_NAME "ODE")

add_helyx_library(${TARGET_NAME} SHARED
    ./findroot/findRoot.C
    ./ODESolvers/ODESolver/ODESolver.C
    ./ODESolvers/ODESolver/ODESolverNew.C
    ./ODESolvers/adaptiveSolver/adaptiveSolver.C
    ./ODESolvers/Euler/Euler.C
    ./ODESolvers/EulerSI/EulerSI.C
    ./ODESolvers/Trapezoid/Trapezoid.C
    ./ODESolvers/RKF45/RKF45.C
    ./ODESolvers/RKCK45/RKCK45.C
    ./ODESolvers/RKDP45/RKDP45.C
    ./ODESolvers/Rosenbrock12/Rosenbrock12.C
    ./ODESolvers/Rosenbrock23/Rosenbrock23.C
    ./ODESolvers/Rosenbrock34/Rosenbrock34.C
    ./ODESolvers/rodas23/rodas23.C
    ./ODESolvers/rodas34/rodas34.C
    ./ODESolvers/SIBS/SIBS.C
    ./ODESolvers/SIBS/SIMPR.C
    ./ODESolvers/SIBS/polyExtrapolate.C
    ./ODESolvers/seulex/seulex.C
    ./ODESystem/ODESystem.C
)

helyx_additional_includes(${TARGET_NAME}
    PRIVATE
        meshTools-includes
)

helyx_link_libraries(${TARGET_NAME} PRIVATE
    finiteVolume
)
