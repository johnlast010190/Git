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
    Foam::faAreaMapper

Description
    FA area mapper.

SourceFiles
    faAreaMapper.C

\*---------------------------------------------------------------------------*/

#ifndef faAreaMapper_H
#define faAreaMapper_H

#include "meshes/polyMesh/polyTopoChangeMap/morphFieldMapper.H"
#include "faMesh/faMesh.H"
#include "meshes/polyMesh/polyTopoChangeMap/faceMapper/faceMapper.H"
#include "containers/HashTables/HashSet/HashSet.H"
#include "meshes/polyMesh/polyTopoChangeMap/polyTopoChangeMap.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class faAreaMapper Declaration
\*---------------------------------------------------------------------------*/

class faAreaMapper
:
    public morphFieldMapper
{
    // Private data

        //- Reference to mesh mapper
        const faMesh& mesh_;

        //- Reference to polyTopoChangeMap
        const polyTopoChangeMap& map_;

        //- Is the mapping direct
        bool direct_;


    // Demand-driven private data

        //- Old mesh size
        label sizeBeforeMapping_;

        //- New face labels after mapping
        mutable labelList* newFaceLabelsPtr_;

        //- New face labels after mapping
        mutable labelList* newFaceLabelsMapPtr_;

        //- Direct addressing (only one form of addressing is used)
        mutable labelList* directAddrPtr_;

        //- Interpolated addressing (only one form of addressing is used)
        mutable labelListList* interpolationAddrPtr_;

        //- Interpolation weights
        mutable scalarListList* weightsPtr_;

        //- Inserted faces
        mutable labelList* insertedObjectLabelsPtr_;



    // Private Member Functions

        //- Calculate addressing
        void calcAddressing() const;

        //- Clear out local storage
        void clearOut();


public:

    // Constructors

        //- Construct from components
        faAreaMapper
        (
            const faMesh& mesh,
            const polyTopoChangeMap& map
        );

        //- Disallow default bitwise copy construction
        faAreaMapper(const faAreaMapper&) = delete;


    //- Destructor
    virtual ~faAreaMapper();


    // Member Functions

        //- Return new face labels
        const labelList& newFaceLabels() const;

        //- Return new face labels map
        //  For new faces return old face index if it exists
        //  If the face has been added, index will be -1
        const labelList& newFaceLabelsMap() const;

        //- Return size
        virtual label size() const
        {
            return newFaceLabels().size();
        }

        //- Return size of field before mapping
        virtual label sizeBeforeMapping() const
        {
            return sizeBeforeMapping_;
        }

        //- Is the mapping direct
        virtual bool direct() const
        {
            return direct_;
        }

        //- Has unmapped elements
        virtual bool hasUnmapped() const
        {
            return insertedObjects();
        }

        //- Return direct addressing
        virtual const labelUList& directAddressing() const;

        //- Return interpolated addressing
        virtual const labelListList& addressing() const;

        //- Return interpolaion weights
        virtual const scalarListList& weights() const;

        //- Are there any inserted faces
        virtual bool insertedObjects() const
        {
            return !insertedObjectLabels().empty();
        }

        //- Return list of inserted faces
        virtual const labelList& insertedObjectLabels() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const faAreaMapper&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
