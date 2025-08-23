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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline const Foam::fvMesh& Foam::ensightMesh::mesh() const
{
    return mesh_;
}


inline const Foam::ensightMesh::options& Foam::ensightMesh::option() const
{
    return *options_;
}


inline Foam::IOstream::streamFormat Foam::ensightMesh::format() const
{
    return options_->format();
}


inline bool Foam::ensightMesh::deprecatedOrder() const
{
    return options_->deprecatedOrder();
}


inline const Foam::ensightCells& Foam::ensightMesh::meshCells() const
{
    return meshCells_;
}


inline const Foam::Map<Foam::word>& Foam::ensightMesh::patches() const
{
    return patchLookup_;
}


inline const Foam::HashTable<Foam::ensightFaces>&
Foam::ensightMesh::boundaryPatchFaces() const
{
    return boundaryPatchFaces_;
}


inline const Foam::HashTable<Foam::ensightFaces>&
Foam::ensightMesh::faceZoneFaces() const
{
    return faceZoneFaces_;
}


inline void Foam::ensightMesh::write(autoPtr<ensightGeoFile>& os) const
{
    write(os.rawRef());
}


// ************************************************************************* //
