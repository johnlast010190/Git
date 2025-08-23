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
    (c) 2011-2013 OpenFOAM Foundation
    (c) 2022-2024 Engys Ltd.

Class
    Foam::cellMapper

Description
    This object provides mapping and fill-in information for cell data
    between the two meshes after the topological change.  It is
    constructed from polyTopoChangeMap.

SourceFiles
    cellMapper.C

\*---------------------------------------------------------------------------*/

#ifndef cellMapper_H
#define cellMapper_H

#include "meshes/polyMesh/polyTopoChangeMap/morphFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class polyTopoChangeMap;

/*---------------------------------------------------------------------------*\
                           Class cellMapper Declaration
\*---------------------------------------------------------------------------*/

class cellMapper
:
    public morphFieldMapper
{
    // Private data

        //- Reference to polyMesh
        const polyMesh& mesh_;

        //- Reference to polyTopoChangeMap
        const polyTopoChangeMap& map_;

        //- Are there any inserted (unmapped) cells
        bool insertedCells_;

        //- Is the mapping direct
        bool direct_;


    // Demand-driven private data

        //- Direct addressing (only one for of addressing is used)
        mutable labelList* directAddrPtr_;

        //- Interpolated addressing (only one for of addressing is used)
        mutable labelListList* interpolationAddrPtr_;

        //- Interpolation weights
        mutable scalarListList* weightsPtr_;

        //- Inserted cells
        mutable labelList* insertedCellLabelsPtr_;


    // Private Member Functions

        //- Calculate addressing for mapping with inserted cells
        void calcAddressing() const;

        //- Clear out local storage
        void clearOut();


public:

    // Constructors

        //- Construct from polyTopoChangeMap
        cellMapper(const polyTopoChangeMap& map);

        //- Disallow default bitwise copy construction
        cellMapper(const cellMapper&) = delete;


    //- Destructor
    virtual ~cellMapper();


    // Member Functions

        //- Return size before mapping
        virtual label sizeBeforeMapping() const override;

        //- Is the mapping direct
        virtual bool direct() const override
        {
            return direct_;
        }

        //- Are there unmapped values? I.e. do all size() elements get value?
        virtual bool hasUnmapped() const override
        {
            return insertedObjects();
        }

        //- Return the direct map addressing
        virtual const labelUList& directAddressing() const override;

        //- Return the interpolated map addressing
        virtual const labelListList& addressing() const override;

        //- Return interpolation weights
        virtual const scalarListList& weights() const override;

        //- Are there any inserted cells
        virtual bool insertedObjects() const override
        {
            return insertedCells_;
        }

        //- Return list of inserted cells
        const virtual labelList& insertedObjectLabels() const override;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const cellMapper&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
