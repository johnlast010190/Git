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
    (c) 2011-2014 OpenFOAM Foundation

Class
    Foam::twoPhaseMixture

Description
    A two-phase mixture model

SourceFiles
    twoPhaseMixture.C

\*---------------------------------------------------------------------------*/

#ifndef twoPhaseMixture_H
#define twoPhaseMixture_H

#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class twoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class twoPhaseMixture
{
protected:

    // Protected data

        word phase1Name_;
        word phase2Name_;

        volScalarField alpha1_;
        volScalarField alpha2_;


public:

    // Constructors

        //- Construct from components
        twoPhaseMixture
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    ~twoPhaseMixture()
    {}


    // Member Functions

        const word& phase1Name() const
        {
            return phase1Name_;
        }

        const word& phase2Name() const
        {
            return phase2Name_;
        }

        //- Return the phase-fraction of phase 1
        const volScalarField& alpha1() const
        {
            return alpha1_;
        }

        //- Return the phase-fraction of phase 1
        volScalarField& alpha1()
        {
            return alpha1_;
        }

        //- Return the phase-fraction of phase 2
        const volScalarField& alpha2() const
        {
            return alpha2_;
        }

        //- Return the phase-fraction of phase 2
        volScalarField& alpha2()
        {
            return alpha2_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
