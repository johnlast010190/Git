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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2016 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "fieldCoordinateSystemTransform/fieldCoordinateSystemTransform.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fields/GeometricFields/transformGeometricField/transformGeometricField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class FieldType>
void Foam::functionObjects::fieldCoordinateSystemTransform::transformField
(
    const FieldType& field
)
{
    word transFieldName(transformFieldName(field.name()));

    store
    (
        transFieldName,
        Foam::transform(dimensionedTensor(coordSysPtr_->R()), field)
    );
}


template<class Type>
void Foam::functionObjects::fieldCoordinateSystemTransform::transform
(
    const word& fieldName
)
{
    typedef VolField<Type> VolFieldType;
    typedef SurfaceField<Type> SurfaceFieldType;

    if (foundObject<VolFieldType>(fieldName))
    {
        DebugInformation
            << type() << ": Field " << fieldName << " already in database"
            << endl;

        transformField<VolFieldType>(lookupObject<VolFieldType>(fieldName));
    }
    else if (foundObject<SurfaceFieldType>(fieldName))
    {
        DebugInformation
            << type() << ": Field " << fieldName << " already in database"
            << endl;

        transformField<SurfaceFieldType>
        (
            lookupObject<SurfaceFieldType>(fieldName)
        );
    }
    else
    {
        IOobject fieldHeader
        (
            fieldName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        if (fieldHeader.typeHeaderOk<VolFieldType>(true, true, false))
        {
            DebugInformation
                << type() << ": Field " << fieldName << " read from file"
                << endl;

            transformField<VolFieldType>
            (
                lookupObject<VolFieldType>(fieldName)
            );
        }
        else if (fieldHeader.typeHeaderOk<SurfaceFieldType>(true, true, false))
        {
            DebugInformation
                << type() << ": Field " << fieldName << " read from file"
                << endl;

            transformField<SurfaceFieldType>
            (
                lookupObject<SurfaceFieldType>(fieldName)
            );
        }
    }
}


// ************************************************************************* //
