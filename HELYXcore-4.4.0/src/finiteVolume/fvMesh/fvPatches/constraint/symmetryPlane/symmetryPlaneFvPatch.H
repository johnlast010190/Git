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
    (c) 2013 OpenFOAM Foundation

Class
    Foam::symmetryPlaneFvPatch

Description
    Symmetry-plane patch.

SourceFiles
    symmetryPlaneFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef symmetryPlaneFvPatch_H
#define symmetryPlaneFvPatch_H

#include "fvMesh/fvPatches/fvPatch/fvPatch.H"
#include "meshes/polyMesh/polyPatches/constraint/symmetryPlane/symmetryPlanePolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class symmetryPlaneFvPatch Declaration
\*---------------------------------------------------------------------------*/

class symmetryPlaneFvPatch
:
    public fvPatch
{
    // Private data

        //- Local reference cast into the symmetryPlane patch
        const symmetryPlanePolyPatch& symmetryPlanePolyPatch_;


public:

    //- Runtime type information
    TypeName(symmetryPlanePolyPatch::typeName_());


    // Constructors

        //- Construct from polyPatch
        symmetryPlaneFvPatch(const polyPatch& patch, const fvBoundaryMesh& bm)
        :
            fvPatch(patch, bm),
            symmetryPlanePolyPatch_
            (
                refCast<const symmetryPlanePolyPatch>(patch)
            )
        {}


    // Member Functions

        //- Return symmetry plane normal
        const vector& n() const
        {
            return symmetryPlanePolyPatch_.n();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
