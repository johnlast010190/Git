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
    (c) held by original author

Class
    Foam::relaxationShapes::relaxationShapeCylindricalMoving

Description
    An object used to define CYLINDRICAL relaxation zones that move with the mesh.
    See e.g.

    @article { jacobsenFuhrmanFredsoe2011,
        Author = {Jacobsen, N G and Fuhrman, D R and Freds\o{}e, J},
        title = {{A Wave Generation Toolbox for the Open-Source CFD Library: OpenFoam\textregistered{}}},
        Journal = {{Int. J. for Numer. Meth. Fluids}},
        Year = {2012},
        Volume = {70},
        Number = {9},
        Pages = {1073-1088},
        DOI = {{10.1002/fld.2726}},
    }
SourceFiles
    relaxationShapeCylindricalMoving.C

Author
    Derived from: Niels Gjøl Jacobsen, Technical University of Denmark.  All rights reserved.
    Created by: Kevin Silva, CSHL, The University of Michigan

\*---------------------------------------------------------------------------*/

#ifndef relaxationShapeCylindricalMoving_H
#define relaxationShapeCylindricalMoving_H

#include "relaxationZone/relaxationShape/relaxationShape.H"
#include "referenceFrames/coordinateFrame.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relaxationShapes
{

/*---------------------------------------------------------------------------*\
                           Class relaxationShapeCylindricalMoving Declaration
\*---------------------------------------------------------------------------*/

class relaxationShapeCylindricalMoving
:
    public relaxationShape
{

protected:

    // Protected data

        //- Coordinate frame
        coordinateFrame* coorFramePtr_;

        // Initial center of relaxation zone
        vector centre0_;

        // Center of relaxation zone at the current time step
        vector centre_;

        // Inner radius of cylindrical relaxatiion zone
        scalar rInner_;

        // Outer radius of cylindrical relaxatiion zone
        scalar rOuter_;

        // Width of the cylindrical relaxatiion zone
        scalar width_;

        // Initial center of domain
        vector domainCentre0_;

    // Protected member functions
        bool insideZone
        (
            const label&
        ) const;

        virtual void computeSigmaCoordinate();

        virtual void findComputationalCells();

public:

    //- Runtime type information
    TypeName("relaxationShapeCylindricalMoving");

    // Constructors

        //- from components
        relaxationShapeCylindricalMoving
        (
            const word& ,
            const fvMesh& mesh_
        );


    // Destructor

        ~relaxationShapeCylindricalMoving()
        {}


    // Member Functions
        virtual const pointField& pointSet();

        virtual scalar interpolation(const scalarField&, const point&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace relaxationShapes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
