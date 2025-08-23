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
    (c) held by original author
    (c) 2024 Engys Ltd.

Class
    Foam::userDefinedDistribution

Description
    The user can define the location of the individual wave gauges. Useful
    for comparison with experimental data sets.

SourceFiles
    userDefinedDistribution.C

Author
    Niels Gjøl Jacobsen, Technical University of Denmark.  All rights reserved.


\*---------------------------------------------------------------------------*/

#ifndef userDefinedDistribution_H
#define userDefinedDistribution_H

#include "preProcessing/probes/pointDistributions/pointDistributions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class waveTheory Declaration
\*---------------------------------------------------------------------------*/

class userDefinedDistribution
:
    public pointDistributions
{
public:

    //- Runtime type information
    TypeName("userDefinedDistribution");


    // Constructors

        //- Construct from components
        userDefinedDistribution
        (
            const fvMesh&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        userDefinedDistribution(const userDefinedDistribution&) = delete;


    //- Destructor
    virtual ~userDefinedDistribution();


    // Member Functions

        virtual pointField evaluate();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const userDefinedDistribution&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
