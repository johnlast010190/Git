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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2017 OpenCFD Ltd
    (c) 2020-2024 Engys Ltd.

Class
    Foam::helyxCoupledControl

Description

\*---------------------------------------------------------------------------*/

#ifndef helyxCoupledControl_H
#define helyxCoupledControl_H

#include "cfdTools/general/solutionControl/solutionControl/solutionControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class helyxCoupledControl Declaration
\*---------------------------------------------------------------------------*/

class helyxCoupledControl
:
    public solutionControl
{
protected:

    // Protected data

        // Solution controls

            //- Flag to indicate whether to solve for the flow
            bool solveFlow_;

            //- Maximum number of PIMPLE correctors
            label nCorr_;

            //- Current PISO corrector
            label corr_;

            //- Flag to indicate whether to only solve turbulence on final iter
            bool turbOnFinalIterOnly_;

            //- Linearization type flag (0 for 1)
            label stabilityMode_;

            //- Temporal term momentum interpolation correction
            Switch choiCorrection_;

            //- Converged flag
            bool converged_;


    // Protected Member Functions

        //- Set the firstIteration flag on the mesh data dictionary
        virtual void setFirstIterFlag
        (
            const bool check = true,
            const bool force = false
        ) override;


public:

    // Static Data Members

        //- Run-time type information
        TypeName("blockSolver");


    // Constructors

        //- Construct from mesh and the name of control sub-dictionary
        helyxCoupledControl
        (
            fvMesh&,
            const word& dictName="blockSolver"
        );

        //- Construct from mesh and the name of control sub-dictionary
        helyxCoupledControl
        (
            fvMesh&,
            const objectRegistry&,
            const word& dictName="blockSolver"
        );

        //- Disallow default bitwise copy construction
        helyxCoupledControl(const helyxCoupledControl&) = delete;


    //- Destructor
    virtual ~helyxCoupledControl();


    // Member Functions

        //- Read controls from fvSolution dictionary
        virtual void readControls() override;

        // Access

            //- Maximum number of outer correctors
            virtual label nOuterCorr() const override
            {
                return nCorr();
            }

            //- Maximum number of correctors
            inline label nCorr() const;

            //- Current PISO corrector index
            inline label corr() const;

            //- Current PISO corrector index
            inline label& corr();

            //- Stability Mode
            inline label stabilityMode() const;

            //- Choi correction
            inline Switch choiCorrection() const;

            inline bool turbOnFinalIterOnly() const
            {
                return turbOnFinalIterOnly_;
            }


        // Solution control

            //- Coupled loop
            virtual bool loop() override;

            //- Pressure corrector loop control
            inline bool correct();

            //- Return true to store the intial residuals
            inline bool storeInitialResiduals() const;

            //- Return true for first outer iteration
            inline bool firstIter() const;

            //- Return true fore final outer iteration
            inline bool finalIter() const;

            //- Return true to solve for flow
            inline bool solveFlow() const;

            //- Return true to solve for turbulence
            inline bool turbCorr();

            //- Return true if all convergence checks are satisfied
            virtual bool criteriaSatisfied() override;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const helyxCoupledControl&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cfdTools/general/solutionControl/helyxCoupledControl/helyxCoupledControlI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
