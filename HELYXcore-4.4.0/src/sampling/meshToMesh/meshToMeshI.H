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
    (c) 2016 OpenCFD Ltd.
    (c) 2012-2018 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "meshToMesh/meshToMesh.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline const Foam::polyMesh& Foam::meshToMesh::srcRegion() const
{
    return srcRegion_;
}


inline const Foam::polyMesh& Foam::meshToMesh::tgtRegion() const
{
    return tgtRegion_;
}


inline const Foam::labelListList& Foam::meshToMesh::srcToTgtCellAddr() const
{
    return srcToTgtCellAddr_;
}


inline const Foam::labelListList& Foam::meshToMesh::tgtToSrcCellAddr() const
{
    return tgtToSrcCellAddr_;
}


inline const Foam::scalarListList& Foam::meshToMesh::srcToTgtCellWght() const
{
    return srcToTgtCellWght_;
}


inline const Foam::scalarListList& Foam::meshToMesh::tgtToSrcCellWght() const
{
    return tgtToSrcCellWght_;
}


inline const Foam::pointListList& Foam::meshToMesh::srcToTgtCellVec() const
{
    return srcToTgtCellVec_;
}


inline const Foam::pointListList& Foam::meshToMesh::tgtToSrcCellVec() const
{
    return tgtToSrcCellVec_;
}


inline Foam::scalar Foam::meshToMesh::V() const
{
    return V_;
}


inline const Foam::autoPtr<Foam::distributionMap>&
Foam::meshToMesh::srcMap() const
{
    return srcMapPtr_;
}


inline const Foam::autoPtr<Foam::distributionMap>&
Foam::meshToMesh::tgtMap() const
{
    return tgtMapPtr_;
}


inline const Foam::PtrList<Foam::AMIInterpolation>&
Foam::meshToMesh::patchAMIs() const
{
    return patchAMIs_;
}


// ************************************************************************* //
