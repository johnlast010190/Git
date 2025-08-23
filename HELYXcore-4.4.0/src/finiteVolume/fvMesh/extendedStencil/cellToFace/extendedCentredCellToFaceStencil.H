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
    (c) 2024 Engys Ltd.

Class
    Foam::extendedCentredCellToFaceStencil

Description

SourceFiles
    extendedCentredCellToFaceStencil.C

\*---------------------------------------------------------------------------*/

#ifndef extendedCentredCellToFaceStencil_H
#define extendedCentredCellToFaceStencil_H

#include "fvMesh/extendedStencil/cellToFace/extendedCellToFaceStencil.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class cellToFaceStencil;

/*---------------------------------------------------------------------------*\
              Class extendedCentredCellToFaceStencil Declaration
\*---------------------------------------------------------------------------*/

class extendedCentredCellToFaceStencil
:
    public extendedCellToFaceStencil
{
    // Private data

        //- Swap map for getting neigbouring data
        autoPtr<distributionMap> mapPtr_;

        //- Per face the stencil.
        labelListList stencil_;


public:

    // Constructors

        //- Construct from uncompacted face stencil
        explicit extendedCentredCellToFaceStencil(const cellToFaceStencil&);

        //- Disallow default bitwise copy construction
        extendedCentredCellToFaceStencil
        (
            const extendedCentredCellToFaceStencil&
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
        template<class T>
        void collectData
        (
            const VolField<T>& fld,
            List<List<T>>& stencilFld
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

        //- Sum vol field contributions to create face values
        template<class Type>
        tmp<SurfaceField<Type>> weightedSum
        (
            const VolField<Type>& fld,
            const List<List<scalar>>& stencilWeights
        ) const
        {
            return extendedCellToFaceStencil::weightedSum
            (
                map(),
                stencil(),
                fld,
                stencilWeights
            );
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const extendedCentredCellToFaceStencil&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
