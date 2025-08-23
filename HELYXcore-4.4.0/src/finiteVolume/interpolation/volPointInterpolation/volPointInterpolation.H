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
    (c) 2016 OpenCFD Ltd.
    (c) 2011-2022 OpenFOAM Foundation
    (c) 2010-2025 Engys Ltd.

Class
    Foam::volPointInterpolation

Description
    Interpolate from cell centres to points (vertices) using inverse distance
    weighting

SourceFiles
    volPointInterpolation.C
    volPointInterpolationTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef volPointInterpolation_H
#define volPointInterpolation_H

#include "meshes/MeshObject/MeshObject.H"
#include "fields/volFields/volFields.H"
#include "fields/GeometricFields/pointFields/pointFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class volPointInterpolation Declaration
\*---------------------------------------------------------------------------*/

class volPointInterpolation
:
    public MeshObject<fvMesh, UpdateableMeshObject, volPointInterpolation>
{
    // Private data

        //- Interpolation scheme weighting factor array
        scalarListList pointWeights_;

        //- Boundary addressing
        autoPtr<primitivePatch> boundaryPtr_;

        //- Per boundary point the weights per pointFaces
        List<scalarListList> boundaryPointWeights_;

        //- Per boundary point the weights per pointFaces.
        List<scalarListList> boundaryPointNbrWeights_;


    // Private Member Functions

        //- Construct all point weighting factors
        void makeWeights();


public:

    // Declare name of the class and its debug switch
    ClassName("volPointInterpolation");


    // Constructors

        //- Constructor given fvMesh
        explicit volPointInterpolation(const fvMesh&);

        //- Disallow default bitwise copy construction
        volPointInterpolation(const volPointInterpolation&) = delete;


    // Member functions

        //- Correct weighting factors for moving mesh
        virtual bool movePoints();

        //- Update mesh topology using the morph engine
        virtual void topoChange(const polyTopoChangeMap&);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);

        //- Redistribute or update using the given distribution map
        virtual void distribute(const polyDistributionMap&);


    // Interpolation functions

        //- Interpolate volField using inverse distance weighting
        //  returning pointField
        template<class Type>
        tmp<PointField<Type>> interpolate
        (
            const VolField<Type>&
        ) const;

        //- Interpolate tmp<volField> using inverse distance weighting
        //  returning pointField
        template<class Type>
        tmp<PointField<Type>> interpolate
        (
            const tmp<VolField<Type>>&
        ) const;


        // Low level

            //- Interpolate from volField to pointField
            //  using inverse distance weighting
            template<class Type>
            void interpolateUnconstrained
            (
                const VolField<Type>&,
                PointField<Type>&
            ) const;

            //- Interpolate from volField to pointField
            //  using inverse distance weighting
            template<class Type>
            void interpolate
            (
                const VolField<Type>&,
                PointField<Type>&
            ) const;

            //- Interpolate volField using inverse distance weighting
            //  returning pointField with name. Optionally caches.
            template<class Type>
            tmp<PointField<Type>> interpolate
            (
                const VolField<Type>&,
                const word& name,
                const bool cache
            ) const;


        // Interpolation for displacement (applies 2D correction)

            //- Interpolate from volField to pointField
            //  using inverse distance weighting
            void interpolateDisplacement
            (
                const volVectorField&,
                pointVectorField&
            ) const;


        // Interpolation for dimensionedFields (needed by foamVtk)

            //- Interpolate dimensioned internal field from cells to points
            //  using inverse distance weighting
            template<class Type>
            void interpolateDimensionedInternalField
            (
                const DimensionedField<Type, volMesh>& vf,
                DimensionedField<Type, pointMesh>& pf
            ) const;

            //- Interpolate dimensionedField using inverse distance weighting
            //  returning pointField
            template<class Type>
            tmp<DimensionedField<Type, pointMesh>> interpolate
            (
                const DimensionedField<Type, volMesh>&
            ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const volPointInterpolation&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "interpolation/volPointInterpolation/volPointInterpolationTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
