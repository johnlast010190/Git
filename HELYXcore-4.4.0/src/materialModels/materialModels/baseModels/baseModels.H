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
    (c) 2021 Engys Ltd.

Class
    Foam::baseModels

Description
    Interface class.

SourceFiles
    baseModels.C

\*---------------------------------------------------------------------------*/

#ifndef baseModels_H
#define baseModels_H

#include "fields/volFields/volFields.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
class baseBaseModels;
struct depSubList
{
    baseBaseModels* model = nullptr;
    UPtrList<baseBaseModels> dependencies;
};
typedef List<depSubList> depList;

/*---------------------------------------------------------------------------*\
                           Class baseModels Declaration
\*---------------------------------------------------------------------------*/

// Class to use for circular dependency check to have same pointer type
class baseBaseModels
{
public:

    //- Destructor
    virtual ~baseBaseModels(){};

    //- Type independent ptr dependency list
    virtual const depList& dep() const = 0;

    //- Type of the template
    virtual const word modelType() const = 0;

    //- Type (scalar/vector/tensor)
    virtual const word funcType() const = 0;

    //- Type defined from TypeName macros
    virtual const word& type() const = 0;
};

template
<
    class Type,
    template<class> class PatchField=fvPatchField,
    class GeoMesh=volMesh
>
class baseModels
:
    public baseBaseModels
{
public:

    //- Runtime type information
    TypeName("baseModels");

    //- Destructor
    virtual ~baseModels()
    {};

    //- Convenience hash table type definition
    typedef HashTable<const dimensionSet*> dimsHashTable;

    //- Pointer to the dimensions hash table
    static dimsHashTable* dimsHashTablePtr_;

    //- Convenience hash table type definition
    typedef HashTable<const bool*> isMixtureHashTable;

    //- Pointer to the mixture identification
    static isMixtureHashTable* isMixtureHashTablePtr_;


    // Constructors

        //- Construcing dimension hash table
        static void constructTables();

        //- Destroy dimension has table
        static void destroyTables();


    // Member functions

        //- Class adding the dimension pointers to the table
        template<class matModType>
        class addToTable
        {
        public:
            addToTable()
            {
                constructTables();
                if
                (
                    !dimsHashTablePtr_->insert
                    (
                        matModType::typeName,
                        &matModType::modelDims
                    )
                )
                {
                    FatalErrorInFunction
                        << "Duplicate entry" << exit(FatalError);
                }
                if
                (
                   !isMixtureHashTablePtr_->insert
                    (
                        matModType::typeName,
                        &matModType::isMixture
                    )
                )
                {
                    FatalErrorInFunction
                        << "Duplicate entry" << exit(FatalError);
                }
            }
            ~addToTable() {destroyTables();}
        };

        //- Internal field
        virtual tmp<Field<Type>> primitiveField() const = 0;

        //- Class for boundary field access boundaryField()[patchi]
        class Boundary
        {
        public:
            virtual ~Boundary(){};
            virtual tmp<Field<Type>> operator[](const label patchi) const = 0;
        };

        //- Boundary field
        virtual Boundary& boundaryField() const = 0;

        //- Dimensions
        virtual const dimensionSet& dimensions() const = 0;

        //- Return dimensions of model from function name
        static const dimensionSet& dimensions(const word& funcName)
        {
            return *dimsHashTablePtr_->operator[](funcName);
        }

        //- Return dimensions of model from function name
        static bool isMixture(const word& funcName)
        {
            if (isMixtureHashTablePtr_->found(funcName))
            {
                return *isMixtureHashTablePtr_->operator[](funcName);
            }
            return false;
        }

        //- Mixture?
        virtual bool mixture() const = 0;

        //- Return whether the returned field is constant
        virtual bool isConst() const = 0;

        //- Single cell access
        virtual Type operator[](const label celli) const = 0;

        //- Calculate with specified values for pressure and temperature
        virtual Type value(scalar p, scalar T) const = 0;

        //- Geometric field access
        virtual tmp<GeometricField<Type, PatchField, GeoMesh>>
        operator()() const = 0;

        //- Update table
        virtual void updateTable(const word&) = 0;

        //- Pointer dependence list
        virtual const depList& dep() const = 0;

        //- Template type (scalar/vector/tensor)
        virtual const word modelType() const = 0;

        //- Property name (rho, Cp, D etc.)
        virtual const word funcType() const = 0;

        //- Is the equation of state is incompressible i.e. rho != f(p)
        virtual bool incompressible() const
        {
            return false;
        }

        //- Is the equation of state is isochoric i.e. rho = const
        virtual bool isochoric() const
        {
            return false;
        }

        //- Return true if thermal conductivity is isotropic
        virtual bool isotropic() const
        {
            return true;
        }

        //- Read function
        virtual bool read() = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


template<class Type, template<class> class PatchField, class GeoMesh>
typename baseModels<Type, PatchField, GeoMesh>::dimsHashTable*
baseModels<Type, PatchField, GeoMesh>::dimsHashTablePtr_ = nullptr;

template<class Type, template<class> class PatchField, class GeoMesh>
typename baseModels<Type, PatchField, GeoMesh>::isMixtureHashTable*
baseModels<Type, PatchField, GeoMesh>::isMixtureHashTablePtr_ = nullptr;


template<class Type,template<class> class PatchField, class GeoMesh>
void baseModels<Type, PatchField, GeoMesh>::constructTables()
{
    static bool constructed = false;
    if (!constructed)
    {
        constructed = true;
        dimsHashTablePtr_ = new dimsHashTable;
        isMixtureHashTablePtr_ = new isMixtureHashTable;
    }
};


template<class Type,template<class> class PatchField, class GeoMesh>
void baseModels<Type, PatchField, GeoMesh>::destroyTables()
{
    if (dimsHashTablePtr_)
    {
        delete dimsHashTablePtr_;
        delete isMixtureHashTablePtr_;
        dimsHashTablePtr_ = nullptr;
        isMixtureHashTablePtr_ = nullptr;
    }
};

template<class Type, template<class> class PatchField, class GeoMesh>
const word baseModels<Type, PatchField, GeoMesh>::typeName("baseModels");

} // End namespace Foam



#endif

// ************************************************************************* //
