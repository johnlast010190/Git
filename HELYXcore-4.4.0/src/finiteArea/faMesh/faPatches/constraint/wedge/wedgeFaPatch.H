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
    (c) 2016 Engys Ltd.


Class
    wedgeFaPatch

Description
    Wedge front and back plane patch.

SourceFiles
    wedgeFaPatch.C

\*---------------------------------------------------------------------------*/

#ifndef wedgeFaPatch_H
#define wedgeFaPatch_H

#include "faMesh/faPatches/faPatch/faPatch.H"
#include "meshes/polyMesh/polyPatches/constraint/wedge/wedgePolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class wedgeFaPatch Declaration
\*---------------------------------------------------------------------------*/

class wedgeFaPatch
:
    public faPatch
{
    // Private data

        const wedgePolyPatch* wedgePolyPatchPtr_;

        //- Axis point label
        mutable label axisPoint_;

        //- Is it axis point looked for?
        mutable bool axisPointChecked_;

        //- Finde axis point
        void findAxisPoint() const;

public:

    //- Runtime type information
    TypeName("wedge");


    // Constructors

        //- Construct from dictionary
        wedgeFaPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const faBoundaryMesh& bm
        );


    // Destructor

        virtual ~wedgeFaPatch()
        {}


    // Member functions

        // Access

            //- Return axis of the wedge
            const vector& axis() const
            {
                return wedgePolyPatchPtr_->axis();
            }

            //- Return plane normal between the wedge boundaries
            const vector& centreNormal() const
            {
                return wedgePolyPatchPtr_->centreNormal();
            }

            //- Return face transformation tensor
            const tensor& edgeT() const
            {
                return wedgePolyPatchPtr_->faceT();
            }

            //- Return neighbour-cell transformation tensor
            const tensor& faceT() const
            {
                return wedgePolyPatchPtr_->cellT();
            }

            //- Return axis point label
            label axisPoint() const
            {
                if (axisPoint_ == -1 && !axisPointChecked_)
                {
                    findAxisPoint();
                }

                return axisPoint_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
