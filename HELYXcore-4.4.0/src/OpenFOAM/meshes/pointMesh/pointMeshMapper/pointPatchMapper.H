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
    Foam::pointPatchMapper

Description
    Mapping class for a pointPatchField.

SourceFiles
    pointPatchMapper.C

\*---------------------------------------------------------------------------*/

#ifndef pointPatchMapper_H
#define pointPatchMapper_H

#include "meshes/pointMesh/pointMeshMapper/pointMapper.H"
#include "fields/pointPatchFields/pointPatchField/pointPatchFieldMapper.H"
#include "meshes/pointMesh/pointPatches/pointPatch/pointPatch.H"
#include "fields/Fields/primitiveFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class pointPatch;
class polyTopoChangeMap;

/*---------------------------------------------------------------------------*\
                       Class pointPatchMapper Declaration
\*---------------------------------------------------------------------------*/

class pointPatchMapper
:
    public generalPointPatchFieldMapper
{
    // Private data

        //- Reference to patch
        const pointPatch& patch_;

        //- Reference to point mapper for internal field
        const morphFieldMapper& pointMapper_;

        //- Reference to mapping data
        const polyTopoChangeMap& map_;

        //- Size before mapping
        const label sizeBeforeMapping_;


    // Demand-driven private data

        mutable bool hasUnmapped_;

        //- Direct addressing (only one for of addressing is used)
        mutable labelList* directAddrPtr_;

        //- Interpolated addressing (only one for of addressing is used)
        mutable labelListList* interpolationAddrPtr_;

        //- Interpolation weights
        mutable scalarListList* weightsPtr_;


    // Private Member Functions

        //- Calculate addressing for mapping with inserted cells
        void calcAddressing() const;

        //- Clear out local storage
        void clearOut();


public:

    // Constructors

        //- Construct from mappers
        pointPatchMapper
        (
            const pointPatch& patch,
            const pointMapper& pointMap,
            const polyTopoChangeMap& map
        );

        //- Disallow default bitwise copy construction
        pointPatchMapper(const pointPatchMapper&) = delete;


    //- Destructor
    virtual ~pointPatchMapper();


    // Member Functions

        //- Return size of field before mapping
        virtual label sizeBeforeMapping() const
        {
            return sizeBeforeMapping_;
        }

        //- Is the mapping direct
        virtual bool direct() const override
        {
            return pointMapper_.direct();
        }

        //- No indirect mapping
        virtual bool indirect() const override
        {
            return false;
        }

        //- Are there unmapped values? I.e. do all size() elements get value?
        virtual bool hasUnmapped() const override
        {
            return hasUnmapped_;
        }

        //- Return the direct map addressing
        virtual const labelUList& directAddressing() const override;

        //- Return the interpolated map addressing
        virtual const labelListList& addressing() const override;

        //- Return interpolation weights
        virtual const scalarListList& weights() const override;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const pointPatchMapper&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
