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
    (c) 2019-2021 OpenCFD Ltd

Class
    Foam::blockEdges::lineEdge

Description
    A straight edge between the start point and the end point.

SourceFiles
    lineEdge.C

\*---------------------------------------------------------------------------*/

#ifndef blockEdges_lineEdge_H
#define blockEdges_lineEdge_H

#include "blockEdges/blockEdge/blockEdge.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace blockEdges
{

/*---------------------------------------------------------------------------*\
                          Class lineEdge Declaration
\*---------------------------------------------------------------------------*/

class lineEdge
:
    public blockEdge
{
public:

    //- Runtime type information
    TypeName("line");


    // Constructors

        //- Construct from components
        lineEdge
        (
            const pointField& points,   //!< Referenced point field
            const edge& fromTo          //!< Start/end in point field
        );

        //- Construct from components
        lineEdge
        (
            const pointField& points,   //!< Referenced point field
            const label from,           //!< Start point in point field
            const label to              //!< End point in point field
        );

        //- Construct from Istream and point field.
        lineEdge
        (
            const dictionary& dict,
            const label index,
            const searchableSurfaces&  /*unused*/,
            const pointField& points,  //!< Referenced point field
            Istream& is
        );


    //- Destructor
    virtual ~lineEdge() = default;


    // Member Functions

        //- The point position corresponding to the curve parameter
        //  0 <= lambda <= 1
        point position(const scalar) const;

        //- The length of the curve
        scalar length() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace blockEdges
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
