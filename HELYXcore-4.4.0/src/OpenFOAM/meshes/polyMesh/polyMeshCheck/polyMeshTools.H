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
    (c) 2012-2013 OpenFOAM Foundation

Namespace
    Foam::polyMeshTools

Description
    Collection of static functions operating on polyMesh (mainly checks) so
    that need access to patch information.

SourceFiles
    polyMeshTools.C

\*---------------------------------------------------------------------------*/
#ifndef polyMeshTools_H
#define polyMeshTools_H

#include "meshes/polyMesh/polyMesh.H"
#include "meshes/primitiveMesh/primitiveMeshCheck/primitiveMeshTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Namespace polyMeshTools Declaration
\*---------------------------------------------------------------------------*/

class polyMeshTools
:
    public primitiveMeshTools
{

public:

    //- Generate orthogonality field. (1 for fully orthogonal, < 1 for
    //  non-orthogonal)
    static tmp<scalarField> faceOrthogonality
    (
        const polyMesh& mesh,
        const vectorField& fAreas,
        const vectorField& cellCtrs
    );

    //- Generate skewness field
    static tmp<scalarField> faceSkewness
    (
        const polyMesh& mesh,
        const pointField& points,
        const vectorField& fCtrs,
        const vectorField& fAreas,
        const vectorField& cellCtrs
    );

    //- Generate interpolation factors field
    static tmp<scalarField> faceWeights
    (
        const polyMesh& mesh,
        const vectorField& fCtrs,
        const vectorField& fAreas,
        const vectorField& cellCtrs
    );

    //- Generate volume ratio field
    static tmp<scalarField> volRatio
    (
        const polyMesh& mesh,
        const scalarField& vol
    );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
