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
    (c) 2024 Engys Ltd.

Class
    Foam::fvBoundaryMeshMapper

Description
    Foam::fvBoundaryMeshMapper

\*---------------------------------------------------------------------------*/

#ifndef fvBoundaryMeshMapper_H
#define fvBoundaryMeshMapper_H

#include "containers/Lists/PtrList/PtrList.H"
#include "fvMesh/fvMeshMapper/fvPatchMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class fvBoundaryMeshMapper Declaration
\*---------------------------------------------------------------------------*/

class fvBoundaryMeshMapper
:
    public PtrList<fvPatchMapper>
{
public:

    // Constructors

        //- Construct from components
        fvBoundaryMeshMapper
        (
            const fvMesh& mesh,
            const faceMapper& faceMap
        )
        :
            PtrList<fvPatchMapper>(mesh.boundary().size())
        {
            const fvBoundaryMesh& patches = mesh.boundary();

            forAll(patches, patchi)
            {
                set
                (
                    patchi,
                    new fvPatchMapper
                    (
                        patches[patchi],
                        faceMap
                    )
                );
            }
        }

        //- Disallow default bitwise copy construction
        fvBoundaryMeshMapper(const fvBoundaryMeshMapper&) = delete;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fvBoundaryMeshMapper&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
