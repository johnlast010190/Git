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
    (c) 1991-2005 OpenCFD Ltd.
    (c) 2016-2024 Engys Ltd.

Class
    lnGradScheme

Description
    Abstract base class for lnGrad schemes.

SourceFiles
    lnGradScheme.C

\*---------------------------------------------------------------------------*/

#ifndef lnGradScheme_H
#define lnGradScheme_H

#include "memory/tmp/tmp.H"
#include "fields/areaFields/areaFieldsFwd.H"
#include "fields/edgeFields/edgeFieldsFwd.H"
#include "db/typeInfo/typeInfo.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class faMesh;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

/*---------------------------------------------------------------------------*\
                 Class lnGradScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class lnGradScheme
:
    public refCount
{
    // Private data

        //- Hold reference to mesh
        const faMesh& mesh_;


public:

    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            tmp,
            lnGradScheme,
            Mesh,
            (
                const faMesh& mesh,
                Istream& schemeData,
                const word& forceGradSchemeName
            ),
            (mesh, schemeData, forceGradSchemeName)
        );


    // Constructors

        //- Construct from mesh
        lnGradScheme(const faMesh& mesh)
        :
            mesh_(mesh)
        {}


    // Selectors

        //- Return new tmp interpolation scheme
        static tmp<lnGradScheme<Type>> New
        (
            const faMesh& mesh,
            Istream& schemeData,
            const word& forceGradSchemeName=word::null
        );


    // Destructor

        virtual ~lnGradScheme();


    // Member Functions

        //- Return mesh reference
        const faMesh& mesh() const
        {
            return mesh_;
        }

        //- Return the lnGrad of the given cell field
        //  with the given weigting factors
        static tmp<GeometricField<Type, faePatchField, faEdgeMesh>>
        lnGrad
        (
            const GeometricField<Type, faPatchField, areaMesh>&,
            const tmp<edgeScalarField>&,
            const word& lnGradName = "lnGrad"
        );

        //- Return the interpolation weighting factors for the given field
        virtual tmp<edgeScalarField> deltaCoeffs
        (
            const GeometricField<Type, faPatchField, areaMesh>&
        ) const = 0;

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return false;
        }

        //- Return the explicit correction to the lnGrad
        //  for the given field
        virtual tmp<GeometricField<Type, faePatchField, faEdgeMesh>>
        correction(const GeometricField<Type, faPatchField, areaMesh>&) const
        {
            return tmp<GeometricField<Type, faePatchField, faEdgeMesh>>(nullptr);
        }

        //- Return the lnGrad of the given cell field
        //  with explicit correction
        virtual tmp<GeometricField<Type, faePatchField, faEdgeMesh>>
        lnGrad(const GeometricField<Type, faPatchField, areaMesh>&) const;

        //- Return the lnGrad of the given tmp cell field
        //  with explicit correction
        tmp<GeometricField<Type, faePatchField, faEdgeMesh>>
        lnGrad
        (
            const tmp<GeometricField<Type, faPatchField, areaMesh>>&
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const lnGradScheme&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables
#define makeLnGradTypeScheme(SS, Type)                                         \
    defineNamedTemplateTypeNameAndDebug(Foam::fa::SS<Foam::Type>, 0);          \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        namespace fa                                                           \
        {                                                                      \
            lnGradScheme<Type>::addMeshConstructorToTable<SS<Type>>           \
                add##SS##Type##MeshConstructorToTable_;                        \
        }                                                                      \
    }

#define makeLnGradScheme(SS)                                                   \
                                                                               \
makeLnGradTypeScheme(SS, scalar)                                               \
makeLnGradTypeScheme(SS, vector)                                               \
makeLnGradTypeScheme(SS, tensor)                                               \
makeLnGradTypeScheme(SS, symmTensor)                                           \
makeLnGradTypeScheme(SS, sphericalTensor)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "finiteArea/lnGradSchemes/lnGradScheme/lnGradScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
