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
    (c) 2016 Engys Ltd.

\*--------------------------------------------------------------------*/

#include "fields/volFields/volFields.H"
#include "meshes/pointMesh/pointMesh.H"
#include "fields/GeometricFields/pointFields/pointFieldsFwd.H"

// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * //

template<class Type>
bool Foam::fieldInitializations::fixedValue::setValue()
{
    typedef VolField<Type> GeoField;

    if (localDb().foundObject<GeoField>(name()))
    {
        GeoField& f
            = const_cast<GeoField&>(localDb().lookupObject<GeoField>(name()));

        f.primitiveFieldRef() = Field<Type>
        (
            word("value"), initDict(), f.size()
        );
        return true;
    }
    else
    {
        return false;
    }
}

template<class Type>
bool Foam::fieldInitializations::fixedValue::setPointValue()
{
    typedef PointField<Type> GeoField;

    if (mesh().foundObject<GeoField>(name(),true))
    {
        GeoField& f
            = const_cast<GeoField&>(mesh().lookupObject<GeoField>(name(),true));
        f.primitiveFieldRef() = Field<Type>
        (
            word("value"), initDict(), f.size()
        );
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
