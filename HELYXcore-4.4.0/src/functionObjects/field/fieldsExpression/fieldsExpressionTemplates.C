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
    (c) 2016 OpenFOAM Foundation
    (c) 2016 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "surfFields/surfFields/surfFields.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type, class FOType>
bool Foam::functionObjects::fieldsExpression::calcFieldTypes(FOType& fo)
{
    typedef VolField<Type> VolFieldType;
    typedef SurfaceField<Type> SurfaceFieldType;
    typedef DimensionedField<Type, surfGeoMesh> SurfFieldType;

    if (foundObject<VolFieldType>(fieldNames_[0]))
    {
        return store
        (
            resultName_,
            fo.template calcFieldType<VolFieldType>()
        );
    }
    else if (foundObject<SurfaceFieldType>(fieldNames_[0]))
    {
        return store
        (
            resultName_,
            fo.template calcFieldType<SurfaceFieldType>()
        );
    }
    else if (foundObject<SurfFieldType>(fieldNames_[0]))
    {
        return store
        (
            resultName_,
            fo.template calcFieldType<SurfFieldType>()
        );
    }
    else
    {
        return false;
    }
}


template<class Type, class FOType>
bool Foam::functionObjects::fieldsExpression::calcType(FOType& fo)
{
    return calcFieldTypes<Type>(fo);
}


template<class FOType>
bool Foam::functionObjects::fieldsExpression::calcAllTypes(FOType& fo)
{
    bool processed = false;

    processed = processed || fo.template calcType<scalar>(fo);
    processed = processed || fo.template calcType<vector>(fo);
    processed = processed || fo.template calcType<sphericalTensor>(fo);
    processed = processed || fo.template calcType<symmTensor>(fo);
    processed = processed || fo.template calcType<tensor>(fo);

    return processed;
}


// ************************************************************************* //
