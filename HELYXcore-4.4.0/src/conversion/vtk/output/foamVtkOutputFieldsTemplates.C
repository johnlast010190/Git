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
    (c) 2016-2107 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


template<class Type>
void Foam::vtk::writeField
(
    vtk::formatter& fmt,
    const DimensionedField<Type, volMesh>& fld
)
{
    const uint64_t payLoad =
    (
        fld.size() * pTraits<Type>::nComponents * sizeof(float)
    );

    fmt.writeSize(payLoad);
    writeList(fmt, fld);

    fmt.flush();
}


template<class Type>
void Foam::vtk::writeField
(
    vtk::formatter& fmt,
    const DimensionedField<Type, volMesh>& fld,
    const labelUList& cellMap
)
{
    const uint64_t payLoad =
    (
        cellMap.size() * pTraits<Type>::nComponents * sizeof(float)
    );

    fmt.writeSize(payLoad);
    writeList(fmt, fld, cellMap);

    fmt.flush();
}


template<class Type>
void Foam::vtk::writeField
(
    vtk::formatter& fmt,
    const VolField<Type>& fld
)
{
    const uint64_t payLoad =
    (
        fld.size() * pTraits<Type>::nComponents * sizeof(float)
    );

    fmt.writeSize(payLoad);
    writeList(fmt, fld.internalField());

    fmt.flush();
}


template<class Type>
void Foam::vtk::writeField
(
    vtk::formatter& fmt,
    const VolField<Type>& fld,
    const labelUList& cellMap
)
{
    const uint64_t payLoad =
    (
        cellMap.size() * pTraits<Type>::nComponents * sizeof(float)
    );

    fmt.writeSize(payLoad);
    writeList(fmt, fld.internalField(), cellMap);

    fmt.flush();
}


// ************************************************************************* //
