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
    (c) 2011-2015 OpenFOAM Foundation

Class
    Foam::motionDiffusivity

Description
    Abstract base class for cell-centre mesh motion diffusivity.

SourceFiles
    motionDiffusivity.C

\*---------------------------------------------------------------------------*/

#ifndef motionDiffusivity_H
#define motionDiffusivity_H

#include "fields/surfaceFields/surfaceFieldsFwd.H"
#include "fvMesh/fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class motionDiffusivity Declaration
\*---------------------------------------------------------------------------*/

class motionDiffusivity
{
    // Private data

        //- Mesh reference
        const fvMesh& mesh_;

public:

    //- Runtime type information
    TypeName("motionDiffusivity");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            motionDiffusivity,
            Istream,
            (
                const fvMesh& mesh,
                Istream& mdData
            ),
            (mesh, mdData)
        );


    // Selectors

        //- Select null constructed
        static autoPtr<motionDiffusivity> New
        (
            const fvMesh& mesh,
            Istream& mdData
        );


    // Constructors

        //- Construct for the given fvMesh
        motionDiffusivity(const fvMesh& mesh);


    //- Destructor
    virtual ~motionDiffusivity();


    // Member Functions

        //- Return reference to the mesh
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        //- Return diffusivity field
        virtual tmp<surfaceScalarField> operator()() const = 0;

        //- Correct the motion diffusivity
        virtual void correct() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
