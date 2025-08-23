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
    Foam::extendedFaceToCellStencil

Description
    Note: transformations on coupled patches not supported. Problem is the
    positions of cells reachable through these patches.

SourceFiles
    extendedFaceToCellStencil.C
    extendedFaceToCellStencilTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef extendedFaceToCellStencil_H
#define extendedFaceToCellStencil_H

#include "meshes/polyMesh/polyDistributionMap/distributionMap.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class globalIndex;

/*---------------------------------------------------------------------------*\
                  Class extendedFaceToCellStencil Declaration
\*---------------------------------------------------------------------------*/

class extendedFaceToCellStencil
{
protected:

    // Protected data

        const polyMesh& mesh_;

public:

    // Constructors

        //- Construct from mesh
        explicit extendedFaceToCellStencil(const polyMesh&);


    // Member Functions

        //- Use map to get the data into stencil order
        template<class T>
        static void collectData
        (
            const distributionMap& map,
            const labelListList& stencil,
            const SurfaceField<T>& fld,
            List<List<T>>& stencilFld
        );

        //- Sum surface field contributions to create cell values
        template<class Type>
        static tmp<VolField<Type>> weightedSum
        (
            const distributionMap& map,
            const labelListList& stencil,
            const SurfaceField<Type>& fld,
            const List<List<scalar>>& stencilWeights
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvMesh/extendedStencil/faceToCell/extendedFaceToCellStencilTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
