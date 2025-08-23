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
    Foam::fv::divScheme

Group
    grpFvDivSchemes

Description
    Abstract base class for div schemes.

SourceFiles
    divScheme.C

\*---------------------------------------------------------------------------*/

#ifndef divScheme_H
#define divScheme_H

#include "memory/tmp/tmp.H"
#include "fields/volFields/volFieldsFwd.H"
#include "fields/surfaceFields/surfaceFieldsFwd.H"
#include "interpolation/surfaceInterpolation/schemes/linear/linear.H"
#include "db/typeInfo/typeInfo.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "blockLduSystem/BlockLduSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
class fvMatrix;

class fvMesh;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                           Class divScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class divScheme
:
    public tmp<divScheme<Type>>::refCount
{
protected:

    // Protected data

        const fvMesh& mesh_;
        tmp<surfaceInterpolationScheme<Type>> tinterpScheme_;


public:

    //- Runtime type information
    virtual const word& type() const = 0;


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            tmp,
            divScheme,
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
        divScheme(const fvMesh& mesh)
        :
            mesh_(mesh),
            tinterpScheme_(new linear<Type>(mesh))
        {}

        //- Construct from mesh and Istream
        divScheme(const fvMesh& mesh, const objectRegistry& db, Istream& is)
        :
            mesh_(mesh),
            tinterpScheme_(surfaceInterpolationScheme<Type>::New(mesh, db, is))
        {}

        //- Disallow default bitwise copy construction
        divScheme(const divScheme&) = delete;


    // Selectors

        //- Return a pointer to a new divScheme created on freestore
        static tmp<divScheme<Type>> New
        (
            const fvMesh& mesh,
            const objectRegistry& db,
            Istream& schemeData
        );


    //- Destructor
    virtual ~divScheme();


    // Member Functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        virtual tmp<VolField<typename innerProduct<vector, Type>::type>> fvcDiv
        (
            const VolField<Type>&
        ) = 0;

        //- Return BlockLduSystem based on div operator
        virtual tmp
        <
            BlockLduSystem<vector, typename innerProduct<vector, Type>::type>
        > fvmDiv
        (
            const VolField<Type>&
        ) const;

        //- Return BlockLduSystem based on div operator
        virtual tmp
        <
            BlockLduSystem<vector, typename innerProduct<vector, Type>::type>
        > fvmBCDiv
        (
            const SurfaceField<scalar>&,
            const VolField<Type>&
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const divScheme&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables

#define makeFvDivTypeScheme(SS, Type)                                          \
    defineNamedTemplateTypeNameAndDebug(Foam::fv::SS<Foam::Type>, 0);          \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        namespace fv                                                           \
        {                                                                      \
            divScheme<Type>::addIstreamConstructorToTable<SS<Type>>            \
                add##SS##Type##IstreamConstructorToTable_;                     \
        }                                                                      \
    }

#define makeFvDivScheme(SS)                                                    \
                                                                               \
makeFvDivTypeScheme(SS, vector)                                                \
makeFvDivTypeScheme(SS, sphericalTensor)                                       \
makeFvDivTypeScheme(SS, symmTensor)                                            \
makeFvDivTypeScheme(SS, tensor)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "finiteVolume/divSchemes/divScheme/divScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
