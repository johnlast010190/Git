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
    (c) 2011-2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "meshes/treeBoundBox/treeBoundBox.H"
#include "containers/Lists/FixedList/FixedList.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


template<unsigned Size>
Foam::treeBoundBox::treeBoundBox
(
    const UList<point>& points,
    const FixedList<label, Size>& indices
)
:
    boundBox(points, indices, false)
{
    // points may be empty, but a FixedList is never empty
    if (points.empty())
    {
        WarningInFunction
            << "cannot find bounding box for zero-sized pointField, "
            << "returning zero" << endl;

        return;
    }
}

// ************************************************************************* //
