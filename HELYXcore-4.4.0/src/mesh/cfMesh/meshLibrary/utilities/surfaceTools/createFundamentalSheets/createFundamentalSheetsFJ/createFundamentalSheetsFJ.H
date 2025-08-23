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
    (c) Creative Fields, Ltd.
    (c) 2024 Engys Ltd.

Authors
    Franjo Juretic (franjo.juretic@c-fields.com)

Class
    createFundamentalSheetsFJ

Description
    Inserts sheets at the boundary of the mesh to capture all feature edges.
    An initial o-layer is inserted first, and followed by the insertion of
    additional sheets which exit the domain at feature edges. The additional
    sheets are inserted such that they exit at other patches next the patch
    they are extruded at.

SourceFiles
    createFundamentalSheetsFJ.C

\*---------------------------------------------------------------------------*/

#ifndef createFundamentalSheetsFJ_H
#define createFundamentalSheetsFJ_H

#include "utilities/meshes/polyMeshGenModifier/polyMeshGenModifier.H"
#include "utilities/surfaceTools/createFundamentalSheets/createFundamentalSheets/createFundamentalSheets.H"
#include "utilities/containers/Lists/labelLongList.H"
#include "primitives/bools/lists/boolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class meshSurfaceEngine;

/*---------------------------------------------------------------------------*\
                Class createFundamentalSheetsFJ Declaration
\*---------------------------------------------------------------------------*/

class createFundamentalSheetsFJ
:
    public createFundamentalSheets
{
    // Private member functions

        //- Create inital sheet from all boundary faces of the surface mesh
        void createInitialSheet();

        //- Create fundamental sheets for all feature edges
        void createSheetsAtFeatureEdges();


public:

    // Runtime type information
    TypeName("modifiedPartitions");


    // Constructors

        //- Construct from mesh data
        createFundamentalSheetsFJ
        (
            polyMeshGen& mesh,
            const bool createWrapperSheet = true
        );

        //- Disallow default construct
        createFundamentalSheetsFJ() = delete;

        //- Disallow default bitwise copy construction
        createFundamentalSheetsFJ(const createFundamentalSheetsFJ&) = delete;


    //- Destructor
    ~createFundamentalSheetsFJ();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const createFundamentalSheetsFJ&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
