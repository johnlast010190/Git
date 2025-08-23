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
    (c) 2022 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

Class
    Foam::fvMeshTopoChangers::list

Description
    Sequence of mesh topology changes applied in order.

SourceFiles
    fvMeshTopoChangersList.C

\*---------------------------------------------------------------------------*/

#ifndef fvMeshTopoChangersList_H
#define fvMeshTopoChangersList_H

#include "fvMesh/fvMeshTopoChangers/fvMeshTopoChanger/fvMeshTopoChanger.H"
#include "containers/Dictionaries/PtrDictionary/PtrDictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshTopoChangers
{

/*---------------------------------------------------------------------------*\
                   Class fvMeshTopoChangers::list Declaration
\*---------------------------------------------------------------------------*/

class list
:
    public fvMeshTopoChanger
{
    // Private Data

        PtrDictionary<fvMeshTopoChanger> list_;


public:

    //- Runtime type information
    TypeName("list");


    // Constructors

        //- Construct from fvMesh and dictionary
        list(fvMesh& mesh, const dictionary& dict);

        //- Disallow default bitwise copy construction
        list(const list&) = delete;


    //- Destructor
    virtual ~list();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();

        //- Update corresponding to the given map
        virtual void topoChange(const polyTopoChangeMap&);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);

        //- Update corresponding to the given distribution map
        virtual void distribute(const polyDistributionMap&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const list&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvMeshTopoChangers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
