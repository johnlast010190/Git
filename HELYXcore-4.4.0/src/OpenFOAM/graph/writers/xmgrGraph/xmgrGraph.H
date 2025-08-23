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
    (c) 2011 OpenFOAM Foundation

Class
    Foam::xmgrGraph

Description
    Output and \b agr file for \em xmgrace
    (http://plasma-gate.weizmann.ac.il/Grace/)

SourceFiles
    xmgrGraph.C

\*---------------------------------------------------------------------------*/

#ifndef xmgrGraph_H
#define xmgrGraph_H

#include "graph/graph.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class xmgrGraph Declaration
\*---------------------------------------------------------------------------*/

class xmgrGraph
:
    public graph::writer
{

public:

    //- Runtime type information
    TypeName("xmgr");

    //- FileName extension  for this graph format
    static const word ext_;


    // Constructors

        //- Construct null
        xmgrGraph()
        {}


    //- Destructor
    ~xmgrGraph()
    {}


    // Member Functions

        // Access

            //- Return the appropriate fileName extension
            //  for this graph format
            const word& ext() const
            {
                return ext_;
            }


        // Write

            void write(const graph&, Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
