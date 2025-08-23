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

Authors
    Franjo Juretic (franjo.juretic@c-fields.com)

Class
    cellIOGraph

Description
    A graph of cells which support automated output

SourceFiles
    cellIOGraph.C

\*---------------------------------------------------------------------------*/

#ifndef cellIOGraph_H
#define cellIOGraph_H

#include "utilities/containers/VRWGraph/VRWGraph.H"
#include "db/regIOobject/regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class cellIOGraph Declaration
\*---------------------------------------------------------------------------*/

class cellIOGraph
:
    public regIOobject,
    public VRWGraph
{

public:

    //- Runtime type information
        TypeName("cellIOGraph");

    // Constructors

        //- Construct from IOobject
        cellIOGraph(const IOobject&);

        //- Construct from IOobject and size of cellIOGraph
        cellIOGraph(const IOobject&, const label);

        //- Construct from IOobject and a VRWGraph
        cellIOGraph(const IOobject&, const VRWGraph&);


    // Member functions

        bool writeData(Ostream&) const;


    // Member operators

        void operator=(const cellIOGraph&);
        void operator=(const VRWGraph&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
