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
    (c) 2011-2019 OpenFOAM Foundation
    (c) 2019-2024 Engys Ltd.

Class
    Foam::fv::gradScheme

Description
    Abstract base class for gradient schemes.

SourceFiles
    gradScheme.C

\*---------------------------------------------------------------------------*/

#ifndef gradScheme_H
#define gradScheme_H

#include "memory/tmp/tmp.H"
#include "fields/volFields/volFieldsFwd.H"
#include "fields/surfaceFields/surfaceFieldsFwd.H"
#include "db/typeInfo/typeInfo.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "blockLduSystem/BlockLduSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                           Class gradScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class gradScheme
:
    public tmp<gradScheme<Type>>::refCount
{
    // Private data

        const fvMesh& mesh_;


public:

    //- Runtime type information
    virtual const word& type() const = 0;


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            tmp,
            gradScheme,
            Istream,
            (
                const fvMesh& mesh,
                const objectRegistry& db,
                Istream& schemeData
            ),
            (mesh, db, schemeData)
        );


    // Constructors

        //- Construct from mesh
        gradScheme(const fvMesh& mesh)
        :
            mesh_(mesh)
        {}

        //- Disallow default bitwise copy construction
        gradScheme(const gradScheme&) = delete;


    // Selectors

        //- Return a pointer to a new gradScheme created on freestore
        static tmp<gradScheme<Type>> New
        (
            const fvMesh& mesh,
            const objectRegistry& db,
            Istream& schemeData
        );


    //- Destructor
    virtual ~gradScheme();


    // Member Functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        //- Calculate and return the grad of the given field.
        //  Used by grad either to recalculate the cached gradient when it is
        //  out of date with respect to the field or when it is not cached.
        virtual tmp<VolField<typename outerProduct<vector, Type>::type>>
        calcGrad
        (
            const VolField<Type>&,
            const word& name
        ) const = 0;

        //- Calculate and return the grad of the given field
        //  which may have been cached
        tmp<VolField<typename outerProduct<vector, Type>::type>> grad
        (
            const VolField<Type>&,
            const word& name
        ) const;

        //- Calculate and return the grad of the given field
        //  with the default name
        //  which may have been cached
        tmp<VolField<typename outerProduct<vector, Type>::type>> grad
        (
            const VolField<Type>&
        ) const;

        //- Calculate and return the grad of the given field
        //  with the default name
        //  which may have been cached
        tmp<VolField<typename outerProduct<vector, Type>::type>> grad
        (
            const tmp<VolField<Type>>&
        ) const;

        //- Return BlockLduSystem based on grad operator
        virtual tmp
        <
            BlockLduSystem<vector, typename outerProduct<vector, Type>::type>
        > fvmGrad
        (
            const VolField<Type>&
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const gradScheme&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables

#define makeFvGradTypeScheme(SS, Type)                                         \
    defineNamedTemplateTypeNameAndDebug(Foam::fv::SS<Foam::Type>, 0);          \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        namespace fv                                                           \
        {                                                                      \
            gradScheme<Type>::addIstreamConstructorToTable<SS<Type>>           \
                add##SS##Type##IstreamConstructorToTable_;                     \
        }                                                                      \
    }


#define makeFvGradScheme(SS)                                                   \
                                                                               \
makeFvGradTypeScheme(SS, scalar)                                               \
makeFvGradTypeScheme(SS, vector)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "finiteVolume/gradSchemes/gradScheme/gradScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
