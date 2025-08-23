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
    Foam::interpolateOvertopping

Description
    Reads a stream of overtopping data.

    The data is interpolated onto an equidistant time axis.

SourceFiles
    interpolateOvertopping.C

Author
    Niels Gjøl Jacobsen, Technical University of Denmark.  All rights reserved.


\*---------------------------------------------------------------------------*/

#ifndef interpolateOvertopping_H
#define interpolateOvertopping_H

#include "postProcessing/postProcessingWaves/postProcessingWaves.H"
#include "postProcessing/postProcessingWaves/writeRawData/rawOvertopping/rawOvertopping.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class waveTheory Declaration
\*---------------------------------------------------------------------------*/

class interpolateOvertopping
:
    public postProcessingWaves
{
public:

    //- Runtime type information
    TypeName("interpolateOvertopping");


    // Constructors

        //- Construct from components
        interpolateOvertopping
        (
            const Time&,
            const dictionary&,
            const word&
        );

        //- Disallow default bitwise copy construction
        interpolateOvertopping(const interpolateOvertopping&) = delete;


    //- Destructor
    virtual ~interpolateOvertopping();


    // Member Functions

        virtual void evaluate();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const interpolateOvertopping&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
