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
    (c) 2017-1017 Engys Ltd.

Class
    Foam::edgeClassification

Description
    Classify feature edge type for every primitive patch edge

SourceFiles
    edgeClassification.C

\*---------------------------------------------------------------------------*/

#ifndef edgeClassification_H
#define edgeClassification_H

#include "meshes/primitiveMesh/primitivePatch/indirectPrimitivePatch.H"
#include "cfdTools/general/include/fvCFD.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class edgeClassification Declaration
\*---------------------------------------------------------------------------*/

template<typename Op>
class packedEltwiseAggOp {
    template<size_t... Is, typename... Args>
    __attribute__((always_inline))
    void apply(
        std::index_sequence<Is...>,
        std::tuple<Args...>& x,
        const std::tuple<Args...>& y
    ) const {
        ((Op{}(std::get<Is>(x), std::get<Is>(y))), ...);
    }

public:
    template<typename T, typename... Args>
    void operator()(T&& x, const std::tuple<Args...>& y) const {
        apply(std::index_sequence_for<Args...>{}, std::forward<T>(x), y);
    }
};

class edgeClassification
{

public:

    // Public data types


        // To combineReduce a labelList. Filters out duplicates.
        class uniqueEqOp
        {
            public:
           void operator()(labelList& x, const labelList& y) const
           {
               if (x.size() == 0)
               {
                   if (y.size() > 0)
                   {
                       x = y;
                   }
               }
               else
               {
                   forAll(y, yi)
                   {
                       if (findIndex(x, y[yi]) == -1)
                       {
                           label sz = x.size();
                           x.setSize(sz+1);
                           x[sz] = y[yi];
                       }
                   }
               }
           }
        };

        class sumEqOp
        {
            public:
            void operator()
            ( pointField& x, const pointField& y ) const
            {
                if (x.size() == 0)
                {
                    x = y;
                }
                else
                {
                    label sz = x.size();
                    label newSz = sz + y.size();
                    x.setSize(newSz);
                    label j = 0;
                    for (label i = sz; i < newSz; i++)
                    {
                        x[i] = y[j];
                        j++;
                    }
                }
            }
        };

        //- Extrusion controls
        enum edgeType
        {
            CONCAVE,
            CONVEX,
            BAFFLE,
            MANIFOLD,
            NONMANIFOLD,
            BOUNDARY
        };

private:

    // Private data

        //- Mesh
        const polyMesh& mesh_;

        const pointField& pts_;

        const indirectPrimitivePatch& pp_;

        const labelList& ppMeshEdges_;

        //- optional flipMap for handling internal faces
        const autoPtr<boolList> flipMap_;

        //- primitive patch feature edge types and angle
        const List<Tuple2<edgeType,scalar>> eType_;

    // Private Member Functions

        List<Tuple2<edgeType,scalar>> setEdgeType
        (
            const scalar convexFeatureAngle,
            const scalar concaveFeatureAngle,
            const scalar baffleFeatureAngle
        );

public:

    //- Runtime type information
    TypeName("edgeClassification");

    // Constructors

        //- Construct from components
        edgeClassification
        (
            const polyMesh& mesh,
            const pointField& pts,
            const indirectPrimitivePatch& pp,
            const labelList& ppMeshEdges,
            const scalar convexFeatureAngle,
            const scalar concaveFeatureAngle,
            const scalar baffleFeatureAngle = -GREAT,
            const autoPtr<boolList> flipMap = nullptr
        );

    // Destructor

    virtual ~edgeClassification(){};

    // Member Functions

    // Calculate point normals (with optional smoothing)
    vectorField calculatePointNormals
    (
        const boolList& excludedFaces,
        const label nSmoothIter,
        const bool correctBoundaryNormals
    );

    //For each edge the patches attached
    List<labelList> uniqueEdgePatches();

    // Return for every primitive patch edge the type and angle
    const List<Tuple2<edgeType,scalar>>& edgeTypes() const
    {
        return eType_;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
