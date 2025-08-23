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
    (c) 2023-2024 Engys Ltd.

Class
    Foam::frameSources

Description

Usage

\*---------------------------------------------------------------------------*/

#ifndef frameSources_H
#define frameSources_H

#include "referenceFrames/coordinateFrame.H"
#include "meshes/MeshObject/MeshObject.H"
#include "cellSetOption/cellSetOption.H"
#include "referenceFrames/frameSourceFaces/frameSourceFaces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fv
{

/*---------------------------------------------------------------------------*\
                        Class frameSources Declaration
\*---------------------------------------------------------------------------*/

class frameSources
:
    public cellSetOption
{

protected:

    // Protected Member Data

        //- Dimensionality check
        void dimensionCheck
        (
            const dimensionSet& dim1,
            const dimensionSet& dim2
        ) const;


    // Protected Member Data

        //- Reference frame
        mutable coordinateFrame* coorFramePtr_;

        //- Frame source faces
        frameSourceFaces frameSourceFaces_;

        //- Should frame use acceleration?
        const Switch acceleration_;

        //- Should frame use conservative formulation?
        const Switch conservative_;

        //- Should frame use ddtPhiCor?
        const Switch ddtPhiCorr_;


public:

    //- Constructors

        //- Mesh constructor
        frameSources
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const objectRegistry& obr
        );

    //- Destructor
        virtual ~frameSources();


    // Public Member Functions

        //- Access coordinate frame of GRF
        coordinateFrame& coorFrame() const
        {
            if (!coorFramePtr_)
            {
                FatalErrorInFunction
                    << "Coordinate frame not allocated"
                    << nl << exit(FatalError);
            }
            return *coorFramePtr_;
        }

        //- Access frame source faces
        const frameSourceFaces& getFrameSourceFaces() const;

        //- Should it use acceleration?
        inline Switch acceleration() const
        {
            return acceleration_;
        }

        //- Should it use conservative formulation?
        inline Switch conservative() const
        {
            return conservative_;
        }

        //- Should it use ddtPhiCor?
        inline Switch ddtPhiCorr() const
        {
            return ddtPhiCorr_;
        }

        //- Check if corriois force is active
        bool accelerationActive() const
        {
            return (acceleration() && active_);
        }

        //- Frame acceleration
        virtual tmp<vectorField> frameAcceleration
        (
            const vectorField& U,
            bool addCoriolis = true
        ) const;

        //- Add coriolis contribution to ddtU therm
        void addAcceleration
        (
            const volVectorField& U,
            volVectorField& ddtU
        ) const override;

        //- Add the Coriolis force contribution to the momentum equation
        //  Adds to the lhs of the equation; optionally add to rhs
        void addAcceleration
        (
            fvVectorMatrix& UEqn,
            const bool rhs = false
        ) const override;

        //- Add the Coriolis force contribution to the momentum equation
        //  Adds to the lhs of the equation; optionally add to rhs
        void addAcceleration
        (
            const volScalarField& rho,
            fvVectorMatrix& UEqn,
            const bool rhs = false
        ) const override;

        //- Add the Coriolis force contribution to the block momentum equation
        void addAcceleration
        (
            fvBlockMatrix<vector>& UEqn
        ) const override;

        //- Add the Coriolis force contribution to the block momentum equation
        void addAcceleration
        (
            const volScalarField& rho,
            fvBlockMatrix<vector>& UEqn
        ) const override;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}//end namespace fv

}//end namespace FOAM

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif
