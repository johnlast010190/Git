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
    (c) 2016-2019 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "vtk/part/foamVtuCells.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


inline enum Foam::vtk::vtuCells::contentType
Foam::vtk::vtuCells::content() const
{
    return output_;
}


inline bool Foam::vtk::vtuCells::decomposeRequested() const
{
    return decomposeRequest_;
}


inline bool Foam::vtk::vtuCells::empty() const
{
    return cellTypes_.empty();
}


inline Foam::label Foam::vtk::vtuCells::size() const
{
    return cellTypes_.size();
}


inline const Foam::List<uint8_t>&
Foam::vtk::vtuCells::cellTypes() const
{
    return cellTypes_;
}


inline const Foam::labelList&
Foam::vtk::vtuCells::vertLabels() const
{
    return vertLabels_;
}


inline const Foam::labelList&
Foam::vtk::vtuCells::vertOffsets() const
{
    return vertOffset_;
}


inline const Foam::labelList&
Foam::vtk::vtuCells::faceLabels() const
{
    return faceLabels_;
}


inline const Foam::labelList&
Foam::vtk::vtuCells::faceOffsets() const
{
    return faceOffset_;
}


inline const Foam::labelList&
Foam::vtk::vtuCells::addPointCellLabels() const
{
    return maps_.additionalIds();
}


inline const Foam::labelList&
Foam::vtk::vtuCells::cellMap() const
{
    return maps_.cellMap();
}


// ************************************************************************* //
