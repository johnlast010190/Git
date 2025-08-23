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
    (c) 2015-2016 OpenCFD Ltd
    (c) 2019 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "surfaceReconstruct/surfaceReconstruct.H"
#include "surfFields/surfFields/surfFields.H"
#include "interpolation/surfaceInterpolation/schemes/linear/linear.H"
#include "finiteVolume/fvc/fvcAverage.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::functionObjects::surfaceReconstruct::reconstructFields()
{
    typedef SurfaceField<Type> SurfaceFieldType;

    // Convert field to map
    HashTable<word> fieldMap(2*fieldSet_.size());
    forAll(fieldSet_, i)
    {
        fieldMap.insert(fieldSet_[i].first(), fieldSet_[i].second());
    }


    HashTable<const SurfaceFieldType*> flds(obr_.lookupClass<SurfaceFieldType>());

    forAllConstIter(typename HashTable<const SurfaceFieldType*>, flds, iter)
    {
        const SurfaceFieldType& fld = *iter();

        if (fieldMap.found(fld.name()))
        {
            // const word sName = "reconstruct(" + fld.name() + ')';
            word& sName = fieldMap[fld.name()];

            if (obr_.found(sName))
            {
                Log << "        updating field " << sName << endl;
            }
            else
            {
                Log << "        reconstructing " << fld.name() << " to create "
                    << sName << endl;
            }

            store(sName, fvc::average(fld));
        }
    }
}


// ************************************************************************* //
