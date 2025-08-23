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
    (c) 1991-2005 OpenCFD Ltd.
    (c) 2024 Engys Ltd.

Class
    Foam::faBoundaryMeshMapper

Description
    Foam::faBoundaryMeshMapper

\*---------------------------------------------------------------------------*/

#ifndef faBoundaryMeshMapper_H
#define faBoundaryMeshMapper_H

#include "containers/Lists/PtrList/PtrList.H"
#include "faMesh/faMeshMapper/faPatchMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class faBoundaryMeshMapper Declaration
\*---------------------------------------------------------------------------*/

class faBoundaryMeshMapper
:
    public PtrList<faPatchMapper>
{
public:

    // Constructors

        //- Construct from components
        faBoundaryMeshMapper
        (
            const faMesh& mesh,
            const polyTopoChangeMap& map
        )
        :
            PtrList<faPatchMapper>(mesh.boundary().size())
        {
            const faBoundaryMesh& patches = mesh.boundary();

            forAll(patches, patchI)
            {
                set
                (
                    patchI,
                    new faPatchMapper
                    (
                        patches[patchI],
                        map
                    )
                );
            }
        }

        //- Disallow default bitwise copy construction
        faBoundaryMeshMapper(const faBoundaryMeshMapper&) = delete;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const faBoundaryMeshMapper&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
