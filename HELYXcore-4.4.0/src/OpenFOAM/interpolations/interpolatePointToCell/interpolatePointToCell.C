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

#include "interpolations/interpolatePointToCell/interpolatePointToCell.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::interpolatePointToCell
(
    const PointField<Type>& ptf,
    const label celli
)
{
    const primitiveMesh& mesh = ptf.mesh()();

    const cell& cFaces = mesh.cells()[celli];

    labelHashSet pointHad(10*cFaces.size());

    Type sum = Zero;

    forAll(cFaces, i)
    {
        const face& f = mesh.faces()[cFaces[i]];

        forAll(f, fp)
        {
            label v = f[fp];

            if (pointHad.insert(v))
            {
                sum += ptf[v];
            }
        }
    }

    return sum/pointHad.size();
}


// ************************************************************************* //
