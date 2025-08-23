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
    (c) 2012-2023 OpenFOAM Foundation
    (c) 2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "fvMeshTools/fvMeshTools.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class GeoField>
void Foam::fvMeshTools::setPatchFields
(
    typename GeoField::Mesh& mesh,
    const label patchi,
    const dictionary& patchFieldDict
)
{
    objectRegistry& obr = const_cast<objectRegistry&>(mesh.thisDb());

    HashTable<GeoField*> fields(obr.lookupClass<GeoField>());

    forAllIter(typename HashTable<GeoField*>, fields, iter)
    {
        GeoField& field = *iter();

        if (GeoField::Mesh::geometryFields.found(field.name())) continue;

        typename GeoField::Boundary& bfield = field.boundaryFieldRef();

        if (patchFieldDict.found(field.name()))
        {
            bfield.set
            (
                patchi,
                GeoField::Patch::New
                (
                    mesh.boundary()[patchi],
                    field(),
                    patchFieldDict.subDict(field.name())
                )
            );
        }
    }
}


template<class GeoField>
void Foam::fvMeshTools::setPatchFields
(
    typename GeoField::Mesh& mesh,
    const label patchi,
    const typename GeoField::value_type& value
)
{
    objectRegistry& obr = const_cast<objectRegistry&>(mesh.thisDb());

    HashTable<GeoField*> fields(obr.lookupClass<GeoField>());

    forAllIter(typename HashTable<GeoField*>, fields, iter)
    {
        GeoField& field = *iter();

        if (GeoField::Mesh::geometryFields.found(field.name())) continue;

        typename GeoField::Boundary& bfield = field.boundaryFieldRef();

        bfield[patchi].forceAssign(value);
    }
}


// ************************************************************************* //
