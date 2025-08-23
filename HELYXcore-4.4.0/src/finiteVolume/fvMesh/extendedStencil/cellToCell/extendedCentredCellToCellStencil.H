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
    (c) 2013-2016 OpenFOAM Foundation

Class
    Foam::extendedCentredCellToCellStencil

Description

SourceFiles
    extendedCentredCellToCellStencil.C

\*---------------------------------------------------------------------------*/

#ifndef extendedCentredCellToCellStencil_H
#define extendedCentredCellToCellStencil_H

#include "fvMesh/extendedStencil/cellToCell/extendedCellToCellStencil.H"
#include "fvMesh/extendedStencil/cellToFace/extendedCellToFaceStencil.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class cellToCellStencil;

/*---------------------------------------------------------------------------*\
              Class extendedCentredCellToCellStencil Declaration
\*---------------------------------------------------------------------------*/

class extendedCentredCellToCellStencil
:
    public extendedCellToCellStencil
{
    // Private data

        //- Swap map for getting neigbouring data
        autoPtr<distributionMap> mapPtr_;

        //- Per cell the stencil.
        labelListList stencil_;


public:

    // Constructors

        //- Construct from uncompacted cell stencil
        explicit extendedCentredCellToCellStencil(const cellToCellStencil&);

        //- Disallow default bitwise copy construction
        extendedCentredCellToCellStencil
        (
            const extendedCentredCellToCellStencil&
        ) = delete;


    // Member Functions

        //- Return reference to the parallel distribution map
        const distributionMap& map() const
        {
            return mapPtr_();
        }

        //- Return reference to the stencil
        const labelListList& stencil() const
        {
            return stencil_;
        }

        //- After removing elements from the stencil adapt the schedule (map).
        void compact();

        //- Use map to get the data into stencil order
        template<class Type>
        void collectData
        (
            const VolField<Type>& fld,
            List<List<Type>>& stencilFld
        ) const
        {
            extendedCellToFaceStencil::collectData
            (
                map(),
                stencil(),
                fld,
                stencilFld
            );
        }

        //- Sum vol field contributions to create cell values
        template<class Type, class WeightType>
        tmp
        <
            GeometricField
            <
                typename outerProduct<WeightType, Type>::type,
                fvPatchField,
                volMesh
            >
        > weightedSum
        (
            const VolField<Type>& fld,
            const List<List<WeightType>>& stencilWeights
        ) const
        {
            return extendedCellToCellStencil::weightedSum
            (
                map(),
                stencil(),
                fld,
                stencilWeights
            );
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const extendedCentredCellToCellStencil&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
