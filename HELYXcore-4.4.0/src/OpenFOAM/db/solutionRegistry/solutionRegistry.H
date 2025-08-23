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
    (c) 2020-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#ifndef solutionRegistry_H
#define solutionRegistry_H

#include "meshes/polyMesh/polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class solutionRegistry
:
    public objectRegistry
{
protected:

    word regionName_;

    word meshName_;

    polyMesh& mesh_;

    bool default_;


public:

    // Constructors

        //- Construct from components
        solutionRegistry
        (
            const Time& runTime,
            polyMesh&,
            const word& regionName,
            bool defaultBehavior = false,
            bool adjointMultipoint = false
        );

        //- Disallow default bitwise copy construction
        solutionRegistry(const solutionRegistry&) = delete;

    //- Destructor
    ~solutionRegistry() override;


    // Member Functions

        //- Mesh name
        inline const word& meshName() const
        {
            return meshName_;
        }

        //- Region name
        inline const word& regionName() const
        {
            return regionName_;
        }

        //- Object registry
        const objectRegistry& registry() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const solutionRegistry&) = delete;
};

}

#endif
