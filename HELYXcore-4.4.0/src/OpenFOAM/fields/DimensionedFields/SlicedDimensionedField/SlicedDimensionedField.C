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
    (c) 2011-2022 OpenFOAM Foundation
    (c) 2022 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/DimensionedFields/SlicedDimensionedField/SlicedDimensionedField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, class GeoMesh>
Foam::SlicedDimensionedField<Type, GeoMesh>::SlicedDimensionedField
(
    const IOobject& io,
    const Mesh& mesh,
    const dimensionSet& ds,
    const Field<Type>& iField
)
:
    DimensionedField<Type, GeoMesh>(io, mesh, ds, Field<Type>())
{
    // Set the internalField to the slice of the complete field
    UList<Type>::shallowCopy
    (
        typename Field<Type>::subField(iField, GeoMesh::size(mesh))
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type, class GeoMesh>
Foam::SlicedDimensionedField<Type, GeoMesh>::~SlicedDimensionedField()
{
    // Set the internalField storage pointer to nullptr before its destruction
    // to protect the field it a slice of.
    UList<Type>::shallowCopy(UList<Type>(nullptr, 0));
}


// ************************************************************************* //
