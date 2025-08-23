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
    createFundamentalSheets

Description
    A base class for various method to generate fundamental sheets
    necessary to capture feature edges

SourceFiles
    createFundamentalSheets.C

\*---------------------------------------------------------------------------*/

#ifndef createFundamentalSheets_H
#define createFundamentalSheets_H

#include "utilities/meshes/polyMeshGenModifier/polyMeshGenModifier.H"
#include "utilities/containers/Lists/labelLongList.H"
#include "primitives/bools/lists/boolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class meshOctree;
class meshSurfaceEngine;

/*---------------------------------------------------------------------------*\
                Class createFundamentalSheets Declaration
\*---------------------------------------------------------------------------*/

class createFundamentalSheets
{
protected:

    // Protected data

        //- type of sheet creation
        word name_;

        //- reference to mesh
        polyMeshGen& mesh_;

        //- shall the procedure create the intial wrapper sheet
        const bool createWrapperSheet_;


public:

    //- Runtime type information
    TypeName("createFundamentalSheets");


    // Constructors

        //- Construct from mesh data
        createFundamentalSheets
        (
            polyMeshGen& mesh,
            const bool createWrapperSheet = true
        );

        //- Disallow default bitwise copy construction
        createFundamentalSheets(const createFundamentalSheets&) = delete;


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        createFundamentalSheets,
        polyMeshGen,
        (
            polyMeshGen& mesh,
            const bool createWrapperSheet
        ),
        (mesh, createWrapperSheet)
    );


    // Selectors

        //- Select constructed from dictionary
        static autoPtr<createFundamentalSheets> New
        (
            polyMeshGen& mesh, const bool createWrapperSheet = true
        );

    //- Destructor
    virtual ~createFundamentalSheets();


    // Member Functions

        //- Return the method name
        const word& name() const
        {
            return name_;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const createFundamentalSheets&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
