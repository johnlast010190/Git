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

Class
    Foam::extendedCellToFaceStencil

Description
    Calculates/constains the extended cell-to-face stencil.

    The stencil is a list of indices into either cells or boundary faces
    in a compact way. (element 0 is owner, 1 is neighbour). The index numbering
    is
    - cells first
    - then all (non-empty patch) boundary faces

    When used in evaluation is a two stage process:
    - collect the data (cell data and non-empty boundaries) into a
    single field
    - (parallel) distribute the field
    - sum the weights*field.

SourceFiles
    extendedCellToFaceStencil.C
    extendedCellToFaceStencilTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef extendedCellToFaceStencil_H
#define extendedCellToFaceStencil_H

#include "meshes/polyMesh/polyDistributionMap/distributionMap.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class globalIndex;

/*---------------------------------------------------------------------------*\
                  Class extendedCellToFaceStencil Declaration
\*---------------------------------------------------------------------------*/

class extendedCellToFaceStencil
{
protected:

    // Protected data

        const polyMesh& mesh_;


    // Protected member functions

        //- Write some statistics about stencil
        static void writeStencilStats
        (
            Ostream& os,
            const labelListList& stencil,
            const distributionMap& map
        );


public:

    // Declare name of the class and its debug switch
    ClassName("extendedCellToFaceStencil");


    // Constructors

        //- Construct from mesh
        explicit extendedCellToFaceStencil(const polyMesh&);


    // Member Functions

        //- Use map to get the data into stencil order
        template<class T>
        static void collectData
        (
            const distributionMap& map,
            const labelListList& stencil,
            const VolField<T>& fld,
            List<List<T>>& stencilFld
        );

        //- Sum vol field contributions to create face values
        template<class Type>
        static tmp<SurfaceField<Type>>
        weightedSum
        (
            const distributionMap& map,
            const labelListList& stencil,
            const VolField<Type>& fld,
            const List<List<scalar>>& stencilWeights
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvMesh/extendedStencil/cellToFace/extendedCellToFaceStencilTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
