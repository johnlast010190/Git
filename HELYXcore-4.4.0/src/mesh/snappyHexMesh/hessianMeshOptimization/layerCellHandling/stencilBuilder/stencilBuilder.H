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
    (c) 2017-2024 Engys Ltd.

Class
    Foam::stencilBuilder

Description
    Abstract class for defining the desirable anisotropy in specific in
    mesh cells

\*---------------------------------------------------------------------------*/

#ifndef stencilBuilder_H
#define stencilBuilder_H

#include "fvMesh/fvMesh.H"
#include "cfdTools/general/include/fvCFD.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

namespace Foam
{

class stencilBuilder
{
public:

    //- Runtime type information
    TypeName("stencilBuilder");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        stencilBuilder,
        dictionary,
        (
            const dictionary& dict,
            const fvMesh& mesh
        ),
        (dict,mesh)
    );


    // Constructors

        //- Construct from mesh
        stencilBuilder(const fvMesh& mesh);

        //- Disallow default bitwise copy construction
        stencilBuilder(const stencilBuilder&) = delete;


    //- Destructor
    virtual ~stencilBuilder();


    // Selectors

        //- Select null constructed
        static autoPtr<stencilBuilder> New
        (
            const dictionary& dict,
            const fvMesh& mesh
        );

    // Public Member Functions

        virtual const vectorField& anisotropicDirection() const = 0;

        virtual const scalarField& anisotropy() const = 0;
};

}//end namespace Foam

#endif
