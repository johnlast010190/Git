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
    (c) 2017 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "vtk/part/foamVtkMeshMaps.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::foamVtkMeshMaps::foamVtkMeshMaps(const label size)
:
    cellMap_(size),
    pointMap_(size),
    additionalIds_(size)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline void Foam::foamVtkMeshMaps::clear()
{
    cellMap_.clear();
    pointMap_.clear();
    additionalIds_.clear();
}


inline const Foam::labelList&
Foam::foamVtkMeshMaps::cellMap() const
{
    return cellMap_;
}


inline Foam::DynamicList<Foam::label>&
Foam::foamVtkMeshMaps::cellMap()
{
    return cellMap_;
}


inline const Foam::labelList&
Foam::foamVtkMeshMaps::pointMap() const
{
    return pointMap_;
}


inline Foam::DynamicList<Foam::label>&
Foam::foamVtkMeshMaps::pointMap()
{
    return pointMap_;
}


inline const Foam::labelList&
Foam::foamVtkMeshMaps::additionalIds() const
{
    return additionalIds_;
}


inline Foam::DynamicList<Foam::label>&
Foam::foamVtkMeshMaps::additionalIds()
{
    return additionalIds_;
}


// ************************************************************************* //
