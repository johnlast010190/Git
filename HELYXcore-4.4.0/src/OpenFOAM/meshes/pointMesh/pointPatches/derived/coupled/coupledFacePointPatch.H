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
    Foam::coupledFacePointPatch

Description
    coupled patch for post-processing.  Used as the base class for processor
    and cyclic pointPatches

SourceFiles
    coupledFacePointPatch.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFacePointPatch_H
#define coupledFacePointPatch_H

#include "meshes/pointMesh/pointPatches/basic/coupled/coupledPointPatch.H"
#include "meshes/pointMesh/pointPatches/facePointPatch/facePointPatch.H"
#include "meshes/polyMesh/polyPatches/basic/coupled/coupledPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class pointBoundaryMesh;

/*---------------------------------------------------------------------------*\
                      Class coupledFacePointPatch Declaration
\*---------------------------------------------------------------------------*/

class coupledFacePointPatch
:
    public facePointPatch,
    public coupledPointPatch
{
    // Private data

        const coupledPolyPatch& coupledPolyPatch_;


    // Private Member Functions

        //- Disallow default construct as copy
        coupledFacePointPatch(const coupledFacePointPatch&);

        //- Disallow default assignment
        void operator=(const coupledFacePointPatch&);


protected:

        // Construction of demand-driven data

            //- Calculate mesh points
            virtual void calcGeometry(PstreamBuffers&) = 0;


public:

    //- Runtime type information
    TypeName(coupledPolyPatch::typeName_());


    // Constructors

        //- Construct from components
        coupledFacePointPatch
        (
            const polyPatch& patch,
            const pointBoundaryMesh& bm
        );


    //- Destructor
    virtual ~coupledFacePointPatch();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
