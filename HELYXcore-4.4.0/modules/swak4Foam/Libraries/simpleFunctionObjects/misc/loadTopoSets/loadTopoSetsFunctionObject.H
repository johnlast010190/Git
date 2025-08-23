/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : dev
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
    (c) ICE Stroemungsfoschungs GmbH
    (c) 2024 Engys Ltd.

Class
    loadTopoSetsFunctionObject

Description
    Load topoSets like cellSets and faceSets in the beginning and register them
    with the mesh

SourceFiles
    loadTopoSetsFunctionObject.C

Contributors/Copyright:
    2016 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#ifndef loadTopoSetsFunctionObject_H
#define loadTopoSetsFunctionObject_H

#include "simpleFunctionObject/simpleFunctionObject.H"
#include "containers/HashTables/HashSet/HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class loadTopoSetsFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class loadTopoSetsFunctionObject
:
    public simpleFunctionObject
{
    // Private Data

        //- Load sets even if they are already in memory
        bool forceLoading_;

        //- Load cellSets
        bool loadCellSets_;

        //- Load faceSets
        bool loadFaceSets_;

        //- Load pointSets
        bool loadPointSets_;

        //- Should we write?
        bool writeSets_;

        HashSet<word> cellSetNames_;
        HashSet<word> faceSetNames_;
        HashSet<word> pointSetNames_;


    // Private Member Functions

        //- load the sets
        void loadSets();

        //- load sets of a specific type
        template<class TopoSetType>
        void loadAllSets(HashSet<word>& names);

        //- write sets of a specific type
        template<class TopoSetType>
        void writeAllSets(HashSet<word>& names);


public:

    //- Runtime type information
    TypeName("loadTopoSets");


    // Constructors

        //- Construct from components
        loadTopoSetsFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        loadTopoSetsFunctionObject(const loadTopoSetsFunctionObject&) = delete;


    // Member Functions

        bool start();

        void writeSimple();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const loadTopoSetsFunctionObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
