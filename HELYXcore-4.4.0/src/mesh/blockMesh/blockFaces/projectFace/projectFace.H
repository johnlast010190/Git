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
    Foam::blockFaces::projectFace

Description
    Projects the given set of face points onto the selected surface of the
    geometry provided as a searchableSurfaces object.

SourceFiles
    projectFace.C

\*---------------------------------------------------------------------------*/

#ifndef blockFaces_projectFace_H
#define blockFaces_projectFace_H

#include "blockFaces/blockFace/blockFace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace blockFaces
{

/*---------------------------------------------------------------------------*\
                           Class projectFace Declaration
\*---------------------------------------------------------------------------*/

class projectFace
:
    public blockFace
{
    // Private data

        //- The surface onto which the points are projected
        const searchableSurface& surface_;


    // Private Member Functions

        const searchableSurface& lookupSurface
        (
            const searchableSurfaces& geometry,
            Istream& is
        ) const;

        //- Convert i,j to single index
        label index
        (
            const labelPair& n,
            const labelPair& coord
        ) const;

        //- Calulate lambdas (but unnormalised)
        void calcLambdas
        (
            const labelPair& n,
            const pointField& points,
            scalarField& lambdaI,
            scalarField& lambdaJ
        ) const;


public:

    //- Runtime type information
    TypeName("project");


    // Constructors

        //- Construct from Istream setting pointsList
        projectFace
        (
            const dictionary& dict,
            const label index,
            const searchableSurfaces& geometry,
            Istream&
        );

        //- Disallow default bitwise copy construction
        projectFace(const projectFace&) = delete;


    //- Destructor
    virtual ~projectFace()
    {}


    // Member Functions

        //- Project the given points onto the surface
        virtual void project
        (
            const blockDescriptor&,
            const label blockFacei,
            pointField& points
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const projectFace&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace blockFaces
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
