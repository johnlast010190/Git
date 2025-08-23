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
    (c) 2011-2012 OpenFOAM Foundation
    (c) 2022-2024 Engys Ltd.

Class
    Foam::pointPatch

Description
    Basic pointPatch represents a set of points from the mesh.

SourceFiles
    pointPatch.C

\*---------------------------------------------------------------------------*/

#ifndef pointPatch_H
#define pointPatch_H

#include "primitives/ints/lists/labelList.H"
#include "fields/Fields/vectorField/vectorField.H"
#include "meshes/primitiveShapes/point/pointField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

class pointBoundaryMesh;
class pointConstraint;
class PstreamBuffers;

/*---------------------------------------------------------------------------*\
                      Class pointPatch Declaration
\*---------------------------------------------------------------------------*/

class pointPatch
{
    // Private data

        //- Reference to boundary mesh
        const pointBoundaryMesh& boundaryMesh_;


protected:

    // Protected Member Functions

        // The pointPatch geometry initialisation is called by pointBoundaryMesh
        friend class pointBoundaryMesh;

        //- Initialise the calculation of the patch geometry
        virtual void initCalcGeometry(PstreamBuffers&)
        {}

        //- Calculate the patch geometry
        virtual void calcGeometry(PstreamBuffers&)
        {}

        //- Initialise the patches for moving points
        virtual void initMovePoints(PstreamBuffers&, const pointField&)
        {}

        //- Correct patches after moving points
        virtual void movePoints(PstreamBuffers&, const pointField&)
        {}

        //- Initialise the update of the patch topology
        virtual void initTopoChange(PstreamBuffers&)
        {}

        //- Update of the patch topology
        virtual void topoChange(PstreamBuffers&)
        {}


public:

    //- Runtime type information
    TypeName("basePatch");


    // Constructor

        pointPatch(const pointBoundaryMesh& bm)
        :
            boundaryMesh_(bm)
        {}

        //- Disallow default bitwise copy construction
        pointPatch(const pointPatch&) = delete;


    //- Destructor
    virtual ~pointPatch()
    {}


    // Member Functions

        //- Return name
        virtual const word& name() const = 0;

        //- Return size
        virtual label size() const = 0;

        //- Return the index of this patch in the pointBoundaryMesh
        virtual label index() const = 0;

        //- Return boundaryMesh reference
        const pointBoundaryMesh& boundaryMesh() const
        {
            return boundaryMesh_;
        }

        //- Return true if this patch field is coupled
        virtual bool coupled() const
        {
            return false;
        }

        //- Return mesh points
        virtual const labelList& meshPoints() const = 0;

        //- Return mesh points
        virtual const vectorField& localPoints() const = 0;

        //- Return  point normals
        virtual const vectorField& pointNormals() const = 0;

        //- Return the constraint type this pointPatch implements.
        virtual const word& constraintType() const
        {
            return word::null;
        }

        //- Accumulate the effect of constraint direction of this patch
        virtual void applyConstraint
        (
            const label pointi,
            pointConstraint&
        ) const
        {}


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const pointPatch&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
