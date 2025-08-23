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
    Foam::relaxationShapes::relaxationShapeFrozen

Description
    This frozen-type relaxation zone is a wrapper relaxationZone for all of the
    other relaxation zones, but it does not update the relaxationZone coeffi-
    cients on moving meshes.
    For certain applications this is desired over the updating, and rather than
    adding an if-sentence to be evaluated at rather low level, this wrapper was
    developed.

SourceFiles
    relaxationShapeFrozen.C

Author
    Niels Gjøl Jacobsen, Deltares.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef relaxationShapeFrozen_H
#define relaxationShapeFrozen_H

#include "relaxationZone/relaxationShape/relaxationShape.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relaxationShapes
{

/*---------------------------------------------------------------------------*\
                           Class relaxationShapeFrozen Declaration
\*---------------------------------------------------------------------------*/

class relaxationShapeFrozen
:
    public relaxationShape
{

protected:

    // Protected data

    // Protected member functions
        virtual void computeSigmaCoordinate();

        virtual void findComputationalCells();

public:

    //- Runtime type information
    TypeName("relaxationShapeFrozen");

    // Constructors

        //- from components
        relaxationShapeFrozen
        (
            const word& ,
            const fvMesh& mesh_
        );


    // Destructor

        ~relaxationShapeFrozen()
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
