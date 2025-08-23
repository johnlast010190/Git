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
    (c) 2015-2019 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "parLagrangianRedistributor.H"
#include "db/Time/Time.H"
#include "db/IOobjectList/IOobjectList.H"
#include "meshes/polyMesh/polyDistributionMap/polyDistributionMap.H"
#include "fields/cloud/cloud.H"
#include "db/IOobjects/CompactIOField/CompactIOField.H"
#include "passiveParticle/passiveParticleCloud.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Container>
Foam::wordList Foam::parLagrangianRedistributor::filterObjects
(
    const IOobjectList& objects,
    const HashSet<word>& selectedFields
)
{
    const word fieldClassName(Container::typeName);

    // Parallel synchronise
    wordList fieldNames(objects.names(fieldClassName));
    Pstream::combineGather(fieldNames, ListUniqueEqOp<word>());
    Pstream::combineScatter(fieldNames);

    if (!selectedFields.empty())
    {
        DynamicList<word> selectedNames(fieldNames.size());
        forAll(fieldNames, i)
        {
            if (selectedFields.found(fieldNames[i]))
            {
                selectedNames.append(fieldNames[i]);
            }
        }
        fieldNames.transfer(selectedNames);
    }
    return fieldNames;
}


template<class Type>
void Foam::parLagrangianRedistributor::redistributeLagrangianFields
(
    const distributionMapBase& map,
    const word& cloudName,
    const IOobjectList& objects,
    const HashSet<word>& selectedFields
) const
{
    const wordList objectNames
    (
        filterObjects<IOField<Type>>
        (
            objects,
            selectedFields
        )
    );

    if (objectNames.size())
    {
        const word fieldClassName(IOField<Type>::typeName);

        Info<< "    Redistributing lagrangian "
            << fieldClassName << "s\n" << endl;

        forAll(objectNames, i)
        {
            Info<< "        " <<  objectNames[i] << endl;

            // Read if present
            IOField<Type> field
            (
                IOobject
                (
                    objectNames[i],
                    srcMesh_.time().timeName(),
                    cloud::prefix/cloudName,
                    srcMesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE,
                    false
                ),
                label(0)
            );

            map.distribute(field);


            if (field.size())
            {
                IOField<Type>
                (
                    IOobject
                    (
                        objectNames[i],
                        tgtMesh_.time().timeName(),
                        cloud::prefix/cloudName,
                        tgtMesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    std::move(field)
                ).write();
            }
        }

        Info<< endl;
    }
}


template<class Type>
void Foam::parLagrangianRedistributor::redistributeLagrangianFieldFields
(
    const distributionMapBase& map,
    const word& cloudName,
    const IOobjectList& objects,
    const HashSet<word>& selectedFields
) const
{
    wordList objectNames
    (
        filterObjects<CompactIOField<Field<Type>, Type>>
        (
            objects,
            selectedFields
        )
    );

    // Append IOField names
    {
        const wordList ioFieldNames
        (
            filterObjects<IOField<Field<Type>>>
            (
                objects,
                selectedFields
            )
        );
        objectNames.append(ioFieldNames);
    }


    if (objectNames.size())
    {
        const word fieldClassName(CompactIOField<Field<Type>, Type>::typeName);

        Info<< "    Redistributing lagrangian "
            << fieldClassName << "s\n" << endl;

        forAll(objectNames, i)
        {
            Info<< "        " <<  objectNames[i] << endl;

            // Read if present
            CompactIOField<Field<Type>, Type > field
            (
                IOobject
                (
                    objectNames[i],
                    srcMesh_.time().timeName(),
                    cloud::prefix/cloudName,
                    srcMesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE,
                    false
                ),
                label(0)
            );

            // Distribute
            map.distribute(field);

            // Write
            if (field.size())
            {
                CompactIOField<Field<Type>, Type>
                (
                    IOobject
                    (
                        objectNames[i],
                        tgtMesh_.time().timeName(),
                        cloud::prefix/cloudName,
                        tgtMesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    std::move(field)
                ).write();
            }
        }
    }
}


template<class Container>
void Foam::parLagrangianRedistributor::readLagrangianFields
(
    const passiveParticleCloud& cloud,
    const IOobjectList& objects,
    const HashSet<word>& selectedFields
)
{
    const wordList objectNames
    (
        filterObjects<Container>
        (
            objects,
            selectedFields
        )
    );

    if (objectNames.size())
    {
        const word fieldClassName(Container::typeName);

        Info<< "    Reading lagrangian "
            << fieldClassName << "s\n" << endl;

        forAll(objectNames, i)
        {
            Info<< "        " <<  objectNames[i] << endl;

            // Read if present
            Container* fieldPtr = new Container
            (
                IOobject
                (
                    objectNames[i],
                    cloud.time().timeName(),
                    cloud,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                label(0)
            );

            fieldPtr->store();
        }
    }
}


template<class Container>
void Foam::parLagrangianRedistributor::redistributeStoredLagrangianFields
(
    const distributionMapBase& map,
    passiveParticleCloud& cloud
) const
{
    HashTable<Container*> fields
    (
        cloud.lookupClass<Container >()
    );

    if (fields.size())
    {
        const word fieldClassName(Container::typeName);

        Info<< "    Redistributing lagrangian "
            << fieldClassName << "s\n" << endl;

        forAllIter(typename HashTable<Container*>, fields, iter)
        {
            Container& field = *iter();

            Info<< "        " <<  field.name() << endl;

            map.distribute(field);

            if (field.size())
            {
                Container
                (
                    IOobject
                    (
                        field.name(),
                        tgtMesh_.time().timeName(),
                        cloud::prefix/cloud.name(),
                        tgtMesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    std::move(field)
                ).write();
            }
        }
    }
}


// ************************************************************************* //
