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
    (c) 2016 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::blockVertices::projectVertex

Description
    Projects the vertex onto the selected surfacees of the
    geometry provided as a searchableSurfaces object.

SourceFiles
    projectVertex.C

\*---------------------------------------------------------------------------*/

#ifndef blockVertices_projectVertex_H
#define blockVertices_projectVertex_H

#include "blockVertices/pointVertex/pointVertex.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace blockVertices
{

/*---------------------------------------------------------------------------*\
                        Class projectVertex Declaration
\*---------------------------------------------------------------------------*/

class projectVertex
:
    public pointVertex
{
    // Private data

        const searchableSurfaces& geometry_;

        //- The indices of surfaces onto which the points are projected
        labelList surfaces_;


public:

    //- Runtime type information
    TypeName("project");


    // Constructors

        //- Construct from Istream setting pointsList
        projectVertex
        (
            const dictionary&,
            const label index,
            const searchableSurfaces& geometry,
            Istream&
        );

        //- Disallow default bitwise copy construction
        projectVertex(const projectVertex&) = delete;


    //- Destructor
    virtual ~projectVertex()
    {}


    // Member Functions

        //- Project the given points onto the surface
        virtual operator point() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const projectVertex&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace blockVertices
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
