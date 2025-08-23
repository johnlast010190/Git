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
    (c) 2024 Engys Ltd.

Class
    Foam::faEdgeMapper

Description
    FA edge mapper.  Currently, edge-based finite area data is not mapped,
    but only resized, since edge-based mapping data is not available

SourceFiles
    faEdgeMapper.C

\*---------------------------------------------------------------------------*/

#ifndef faEdgeMapper_H
#define faEdgeMapper_H

#include "meshes/polyMesh/polyTopoChangeMap/morphFieldMapper.H"
#include "faMesh/faMesh.H"
#include "meshes/polyMesh/polyTopoChangeMap/faceMapper/faceMapper.H"
#include "containers/HashTables/HashSet/HashSet.H"
#include "meshes/polyMesh/polyTopoChangeMap/polyTopoChangeMap.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class faEdgeMapper Declaration
\*---------------------------------------------------------------------------*/

class faEdgeMapper
:
    public morphFieldMapper
{
    // Private data

        //- Reference to mesh
        const faMesh& mesh_;

        //- Old mesh size
        label sizeBeforeMapping_;


    // Demand-driven private data

        //- Direct addressing
        mutable labelList* directAddrPtr_;


    // Private Member Functions

        //- Calculate addressing
        void calcAddressing() const;

        //- Clear out local storage
        void clearOut();


public:

    // Constructors

        //- Construct from components
        faEdgeMapper
        (
            const faMesh& mesh,
            const polyTopoChangeMap& map
        );

        //- Disallow default bitwise copy construction
        faEdgeMapper(const faEdgeMapper&) = delete;


    //- Destructor
    virtual ~faEdgeMapper();


    // Member Functions

        //- Return size
        virtual label size() const
        {
            return mesh_.nInternalEdges();
        }

        //- Return size of field before mapping
        virtual label sizeBeforeMapping() const
        {
            return sizeBeforeMapping_;
        }

        //- Is the mapping direct
        virtual bool direct() const
        {
            return true;
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
            return false;
        }

        //- Return list of inserted faces
        virtual const labelList& insertedObjectLabels() const
        {
            return labelList::null();
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const faEdgeMapper&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
