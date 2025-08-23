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
    (c) 2011-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "matrices/LduMatrix/Solvers/PCICG/PCICG.H"
#include "matrices/LduMatrix/Solvers/PBiCCCG/PBiCCCG.H"
#include "matrices/LduMatrix/Solvers/PBiCICG/PBiCICG.H"
#include "matrices/LduMatrix/Solvers/SmoothSolver/SmoothSolver.H"
#include "fields/Fields/fieldTypes.H"

#define makeLduSolvers(Type, DType, LUType)                                    \
                                                                               \
    makeLduSolver(DiagonalSolver, Type, DType, LUType);                        \
    makeLduSymSolver(DiagonalSolver, Type, DType, LUType);                     \
    makeLduAsymSolver(DiagonalSolver, Type, DType, LUType);                    \
                                                                               \
    makeLduSolver(PCICG, Type, DType, LUType);                                 \
    makeLduSymSolver(PCICG, Type, DType, LUType);                              \
                                                                               \
    makeLduSolver(PBiCCCG, Type, DType, LUType);                               \
    makeLduAsymSolver(PBiCCCG, Type, DType, LUType);                           \
                                                                               \
    makeLduSolver(PBiCICG, Type, DType, LUType);                               \
    makeLduAsymSolver(PBiCICG, Type, DType, LUType);                           \
                                                                               \
    makeLduSolver(SmoothSolver, Type, DType, LUType);                          \
    makeLduSymSolver(SmoothSolver, Type, DType, LUType);                       \
    makeLduAsymSolver(SmoothSolver, Type, DType, LUType);

namespace Foam
{
    makeLduSolvers(scalar, scalar, scalar);
    makeLduSolvers(vector, scalar, scalar);
    makeLduSolvers(sphericalTensor, scalar, scalar);
    makeLduSolvers(symmTensor, scalar, scalar);
    makeLduSolvers(tensor, scalar, scalar);
};


// ************************************************************************* //
