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
    (c) 2010-2024 Engys Ltd.

Class
    Foam::pointVolInterpolation

Description
    Interpolate from cell centres to points (vertices) using inverse distance
    weighting

SourceFiles
    pointVolInterpolation.C
    pointVolInterpolate.C

\*---------------------------------------------------------------------------*/

#ifndef pointVolInterpolation_H
#define pointVolInterpolation_H

#include "meshes/MeshObject/MeshObject.H"
#include "primitives/Scalar/lists/scalarList.H"
#include "fields/volFields/volFields.H"
#include "fields/GeometricFields/pointFields/pointFields.H"
#include "primitives/bools/Switch/Switch.H"
#include "interpolations/primitivePatchInterpolation/primitivePatchInterpolation.H"
#include "fields/FieldFields/FieldField/FieldFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;
class pointMesh;

/*---------------------------------------------------------------------------*\
                       Class pointVolInterpolation Declaration
\*---------------------------------------------------------------------------*/

class pointVolInterpolation
:
    public MeshObject<fvMesh, UpdateableMeshObject, pointVolInterpolation>
{
    // Private data

        const pointMesh& pointMesh_;
        const fvMesh& fvMesh_;

        //- Interpolation scheme weighting factor array.
        mutable FieldField<Field, scalar>* volWeightsPtr_;

        //- Primitive patch interpolators
        mutable PtrList<primitivePatchInterpolation>* patchInterpolatorsPtr_;


    // Private member functions

        //- Return patch interpolators
        const PtrList<primitivePatchInterpolation>& patchInterpolators() const;

        //- Construct point weighting factors
        void makeWeights() const;

        //- Clear addressing
        void clearAddressing() const;

        //- Clear geometry
        void clearGeom() const;


    protected:
        const pointMesh& pMesh() const
        {
            return pointMesh_;
        }
        const fvMesh& vMesh() const
        {
            return fvMesh_;
        }


    public:

        // Declare name of the class and it's debug switch
        ClassName("pointVolInterpolation");

    // Constructors
        //- Constructor given pointMesh and fvMesh.
        //pointVolInterpolation(const pointMesh&, const fvMesh&);
        pointVolInterpolation(const fvMesh&);

    // Destructor
        ~pointVolInterpolation();

    // Member functions

        // Access
            //- Return reference to weights arrays.
            //  This also constructs the weighting factors if necessary.
            const FieldField<Field, scalar>& volWeights() const;

        // Edit
            //- Correct weighting factors for moving mesh
            virtual bool movePoints();

            //- Update mesh topology using the morph engine
            virtual void topoChange(const polyTopoChangeMap&);

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap&);

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap&);


    // Interpolation functions

        //- Interpolate from pointField to volField
        //  using inverse distance weighting
        template<class Type>
        void interpolate
        (
            const GeometricField<Type, pointPatchField, pointMesh>&,
            GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        //- Interpolate pointField returning volField
        //  using inverse distance weighting
        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh>> interpolate
        (
            const GeometricField<Type, pointPatchField, pointMesh>&
        ) const;

        //- Interpolate tmp<pointField> returning volField
        //  using inverse distance weighting
        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh>> interpolate
        (
            const tmp<GeometricField<Type, pointPatchField, pointMesh>>&
        ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "pointVolInterpolate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
