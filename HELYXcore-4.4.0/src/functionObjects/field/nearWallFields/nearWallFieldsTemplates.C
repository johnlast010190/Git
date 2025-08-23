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
    (c) 2015 OpenCFD Ltd.
    (c) 2011-2023 OpenFOAM Foundation
    (c) 2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "nearWallFields/nearWallFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::functionObjects::nearWallFields::createFields
(
    PtrList<VolField<Type>>& sfields
) const
{
    const UPtrList<VolField<Type>> fields(mesh_.fields<VolField<Type>>());

    forAll(fields, i)
    {
        const VolField<Type>& field = fields[i];

        if (fieldMap_.found(field.name()))
        {
            const word& sampleFieldName = fieldMap_[field.name()];

            if (obr_.found(sampleFieldName))
            {
                WarningInFunction
                    << "    a field named " << sampleFieldName
                    << " already exists on the mesh"
                    << endl;
            }
            else
            {
                label sz = sfields.size();
                sfields.setSize(sz + 1);

                IOobject io(field);
                io.readOpt() = IOobject::NO_READ;
                io.writeOpt() = IOobject::NO_WRITE;
                io.rename(sampleFieldName);

                sfields.set(sz, new VolField<Type>(io, field));

                Log << "    created " << sfields[sz].name()
                    << " to sample " << field.name() << endl;
            }
        }
    }
}


template<class Type>
void Foam::functionObjects::nearWallFields::sampleBoundaryField
(
    const interpolationCellPoint<Type>& interpolator,
    VolField<Type>& fld
) const
{
    // Construct flat fields for all patch faces to be sampled
    Field<Type> sampledValues(getPatchDataMapPtr_().constructSize());

    forAll(cellToWalls_, celli)
    {
        const labelList& cData = cellToWalls_[celli];

        forAll(cData, i)
        {
            const point& samplePt = cellToSamples_[celli][i];
            sampledValues[cData[i]] = interpolator.interpolate(samplePt, celli);
        }
    }

    // Send back sampled values to patch faces
    getPatchDataMapPtr_().reverseDistribute
    (
        getPatchDataMapPtr_().constructSize(),
        sampledValues
    );

    typename VolField<Type>::
        Boundary& fldBf = fld.boundaryFieldRef();

    // Pick up data
    label nPatchFaces = 0;
    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        label patchi = iter.key();

        fvPatchField<Type>& pfld = fldBf[patchi];

        Field<Type> newFld(pfld.size());
        forAll(pfld, i)
        {
            newFld[i] = sampledValues[nPatchFaces++];
        }

        pfld.forceAssign(newFld);
    }
}


template<class Type>
void Foam::functionObjects::nearWallFields::sampleFields
(
    PtrList<VolField<Type>>& sflds
) const
{
    typedef VolField<Type> VolFieldType;

    forAll(sflds, i)
    {
        const word& fldName = reverseFieldMap_[sflds[i].name()];
        const VolFieldType& fld = lookupObject<VolFieldType>(fldName);

        // Take over internal and boundary values
        sflds[i].forceAssign(fld);

        // Construct interpolation method
        interpolationCellPoint<Type> interpolator(fld);

        // Override sampled values
        sampleBoundaryField(interpolator, sflds[i]);
    }
}


// ************************************************************************* //
