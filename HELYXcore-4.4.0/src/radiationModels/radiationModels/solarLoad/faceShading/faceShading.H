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
    (c) 2024 Engys Ltd.

Class
    Foam::faceShading

Description

    faceShading uses the transmissivity value in the boundaryRadiationProperties
    in order to evaluate which faces are "hit" by the "direction" vector.
    NOTE: Only transmissivity values of zero are considered for opaque walls.

SourceFiles
    faceShading.C

\*---------------------------------------------------------------------------*/

#ifndef faceShading_H
#define faceShading_H

#include "fvMesh/fvMesh.H"
#include "db/Time/Time.H"
#include "meshTools/meshTools.H"
#include "fields/Fields/DynamicField/DynamicField.H"
#include "primitives/ints/lists/labelIOList.H"
#include "meshes/polyMesh/polyPatches/derived/wall/wallPolyPatch.H"
#include "triSurface/triSurfaceTools/triSurfaceTools.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class faceShading Declaration
\*---------------------------------------------------------------------------*/

class faceShading
{
    // Private data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Direction
        vector direction_;

        //- Faces directly hit by vector direction
        labelList rayStartFaces_;

        //- Absolute temperature offset
        const scalar TRef_;


    // Private members

        //- Calculate ray start faces
        void calculate();


        //- Construct a triSurface from patches and faces on global local index
        triSurface triangulate
        (
            const labelHashSet& includePatches,
            const List<labelHashSet>& includeAllFaces
        );

        //- Write rays
        void writeRays
        (
            const fileName& fName,
            const DynamicField<point>& endCf,
            const pointField& myFc
        );


public:

    // Declare name of the class and its debug switch
    ClassName("faceShading");


    // Constructors

        //- Construct from components
        faceShading
        (
            const fvMesh& mesh,
            const vector dir,
            const labelList& hitFaceList,
            const scalar& TRef
        );

        //- Construct from mesh and vector
        faceShading
        (
            const fvMesh& mesh,
            const vector,
            const scalar& TRef
        );

        //- Disallow default bitwise copy construction
        faceShading(const faceShading&) = delete;


    //- Destructor
    ~faceShading();


    // Member Functions

        //- Const acess to direction
        const vector direction() const
        {
            return direction_;
        }

        //- Non-const access to direction
        vector& direction()
        {
            return direction_;
        }

        //- Access to rayStartFaces
        const labelList& rayStartFaces() const
        {
            return rayStartFaces_;
        }

        //- Recalculate rayStartFaces using direction vector
        void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const faceShading&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
