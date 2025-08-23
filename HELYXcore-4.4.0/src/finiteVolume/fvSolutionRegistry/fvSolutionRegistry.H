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
    (c) 2023-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#ifndef fvSolutionRegistry_H
#define fvSolutionRegistry_H

#include "db/solutionRegistry/solutionRegistry.H"
#include "fvMesh/fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class fvSolutionRegistry Declaration
\*---------------------------------------------------------------------------*/

class fvSolutionRegistry
:
    public solutionRegistry
{
protected:

    fvMesh& mesh_;


public:

    // Constructors

        //- Construct from components
        fvSolutionRegistry
        (
            const Time& runTime,
            fvMesh&,
            const word& regionName,
            bool defaultBehavior = false,
            bool adjointMultipoint = false
        );

        //- Disallow default bitwise copy construction
        fvSolutionRegistry(const fvSolutionRegistry&) = delete;


    // Member Functions

        //- Access to mesh
        inline fvMesh& mesh() const
        {
            return mesh_;
        }

        //- Get refference to the mesh based on if the object registry
        //  is a fvSolutionRegistry
        static const fvMesh& getMesh(const objectRegistry& obr)
        {
            if (isA<fvSolutionRegistry>(obr))
            {
                return dynamic_cast<const fvSolutionRegistry&>(obr).mesh();
            }
            return dynamic_cast<const fvMesh&>(obr);
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fvSolutionRegistry&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
