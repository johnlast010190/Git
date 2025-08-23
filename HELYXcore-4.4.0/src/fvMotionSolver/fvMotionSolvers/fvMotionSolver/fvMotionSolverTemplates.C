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
    (c) 2012-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "fvMotionSolvers/fvMotionSolver/fvMotionSolver.H"
#include "fields/pointPatchFields/basic/fixedValue/fixedValuePointPatchFields.H"
#include "fvPatchFields/derived/cellMotion/cellMotionFvPatchFields.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::wordList Foam::fvMotionSolver::cellMotionBoundaryTypes
(
    const typename PointField<Type>::
    Boundary& pmUbf
) const
{
    wordList cmUbf = pmUbf.types();

    // Remove global patches from the end of the list
    cmUbf.setSize(fvMesh_.boundary().size());

    forAll(cmUbf, patchi)
    {
        if (isA<fixedValuePointPatchField<Type>>(pmUbf[patchi]))
        {
            cmUbf[patchi] = cellMotionFvPatchField<Type>::typeName;
        }

        if (debug)
        {
            Pout<< "Patch:" << fvMesh_.boundary()[patchi].patch().name()
                << " pointType:" << pmUbf.types()[patchi]
                << " cellType:" << cmUbf[patchi] << endl;
        }
    }

    return cmUbf;
}


// ************************************************************************* //
