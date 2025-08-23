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
    faceIOGraph

Description
    A graph of faces which support automated output

SourceFiles
    faceIOGraph.C

\*---------------------------------------------------------------------------*/

#ifndef faceIOGraph_H
#define faceIOGraph_H

#include "utilities/containers/VRWGraph/VRWGraph.H"
#include "db/regIOobject/regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class faceIOGraph Declaration
\*---------------------------------------------------------------------------*/

class faceIOGraph
:
    public regIOobject,
    public VRWGraph
{

public:

    //- Runtime type information
        TypeName("faceIOGraph");

    // Constructors

        //- Construct from IOobject
        faceIOGraph(const IOobject&);

        //- Construct from IOobject and size of faceIOGraph
        faceIOGraph(const IOobject&, const label);

        //- Construct from IOobject and a VRWGraph
        faceIOGraph(const IOobject&, const VRWGraph&);


    // Member functions

        bool writeData(Ostream&) const;


    // Member operators

        void operator=(const faceIOGraph&);
        void operator=(const VRWGraph&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
