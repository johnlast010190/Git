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
    (c) 2015-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/fvPatchField/fieldMappers/gibFvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type, template<class> class PatchField, class GeoMesh>
void GIBMapping::StoreOldFieldsToPatch
(
    List<Tuple2<word, Field<Type>>>& ptrFields,
    const label& pI
)
{
    typedef GeometricField<Type, PatchField, GeoMesh> GeoField;

    const UPtrList<GeoField> fields(GIBChanger_.mesh().fields<GeoField>());

    ptrFields.resize(fields.size());

    forAll(fields, i)
    {
        const GeoField& field = fields[i];

        const List<Field<Type>>& bfield = field.boundaryField().oldTimeField();
        Type defaultValue = field.boundaryField()[pI].defaultGIBValue();

        ptrFields[i].first() = field.name();
        if (pI == GIBChanger_.masterId())
        {
            ptrFields[i].second().setSize(sourcePatchm_->size());
        }
        else if (pI == GIBChanger_.slaveId())
        {
            ptrFields[i].second().setSize(sourcePatchs_->size());
        }
        else
        {
            FatalErrorInFunction
                << abort(FatalError);
        }

        forAll(ptrFields[i].second(), fI)
        {
            if (fI < bfield[pI].size())
            {
                ptrFields[i].second()[fI] = bfield[pI][fI];
            }
            else
            {
                ptrFields[i].second()[fI] = defaultValue;
            }
        }
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void GIBMapping::MapGIBField
(
    List<Tuple2<word, Field<Type>>>& ptrFields,
    const label& pI,
    const AMIInterpolation& ipToP
)
{
    typedef GeometricField<Type, PatchField, GeoMesh> GeoField;

    UPtrList<GeoField> fields(GIBChanger_.mesh().fields<GeoField>());

    forAll(fields, i)
    {
        GeoField& field = fields[i];

        typename GeoField::Boundary& bfield = field.boundaryFieldRef();

        Field<Type>& f = bfield[pI];

        forAll(ptrFields, j)
        {
            if (ptrFields[j].first() == field.name())
            {
                const Field<Type>& oldFaceField = ptrFields[j].second();

                tmp<Field<Type>> fnewt =
                    ipToP.interpolateToTarget(oldFaceField);
                f = fnewt();

                // Ask BC to map any stored fields
                gibFvPatchFieldMapper mapper(ipToP);
                bfield[pI].autoMapGIB(mapper);
            }
        }
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void GIBMapping::Map1to1GIBField
(
    const List<Tuple2<word, Field<Type>>>& ptrmFields,
    const List<Tuple2<word, Field<Type>>>& ptrsFields
)
{
    typedef GeometricField<Type, PatchField, GeoMesh> GeoField;

    UPtrList<GeoField> fields(GIBChanger_.mesh().fields<GeoField>());

    const label& mId = GIBChanger_.masterId();
    const label& sId = GIBChanger_.slaveId();

    forAll(fields, i)
    {
        GeoField& field = fields[i];

        typename GeoField::Boundary& bfield = field.boundaryFieldRef();

        forAll(ptrmFields, j)
        {
            if (ptrmFields[j].first() == field.name())
            {
                Field<Type>& fm = bfield[mId];
                Field<Type>& fs = bfield[sId];

                const Field<Type>& moldFaceField = ptrmFields[j].second();
                const Field<Type>& soldFaceField = ptrsFields[j].second();

                fm = moldFaceField;
                fs = soldFaceField;
            }
        }
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void GIBMapping::MapNonOverlapFaces
(
    const List<Tuple2<word, Field<Type>>>& ptrmFields,
    const label mId,
    const labelList& dNonOverFaces,
    const labelList& nonOverFaces
)
{
    typedef GeometricField<Type, PatchField, GeoMesh> GeoField;

    UPtrList<GeoField> fields(GIBChanger_.mesh().fields<GeoField>());

    forAll(fields, i)
    {
        GeoField& field = fields[i];

        typename GeoField::Boundary& bfield = field.boundaryFieldRef();

        forAll(ptrmFields, j)
        {
            if (ptrmFields[j].first() == field.name())
            {
                Field<Type>& fm = bfield[mId];

                const Field<Type>& moldFaceField = ptrmFields[j].second();

                Type value = gAverage(moldFaceField);

                forAll(nonOverFaces, fI)
                {
                    label dfII = dNonOverFaces[fI];
                    label fII = nonOverFaces[fI];

                    if (dfII != -1)
                    {
                        if (fII != -1)
                        {
                            fm[dfII] = moldFaceField[fII];
                        }
                        else
                        {
                            fm[dfII] = value;
                        }
                    }
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
