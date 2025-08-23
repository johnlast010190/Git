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
    Foam::interpolateSurfaceElevation

Description
    Reads a stream of logged surface elevations and writes data for each gauge
    in an OF-formatted manner.

    The data is interpolated onto an equidistant time axis.

SourceFiles
    interpolateSurfaceElevation.C

Author
    Niels Gjøl Jacobsen, Technical University of Denmark.  All rights reserved.


\*---------------------------------------------------------------------------*/

#ifndef interpolateSurfaceElevation_H
#define interpolateSurfaceElevation_H

#include "postProcessing/postProcessingWaves/postProcessingWaves.H"
#include "postProcessing/postProcessingWaves/writeRawData/rawSurfaceElevation/rawSurfaceElevation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class waveTheory Declaration
\*---------------------------------------------------------------------------*/

class interpolateSurfaceElevation
:
    public postProcessingWaves
{
public:

    //- Runtime type information
    TypeName("interpolateSurfaceElevation");


    // Constructors

        //- Construct from components
        interpolateSurfaceElevation
        (
            const Time&,
            const dictionary&,
            const word&
        );

        //- Disallow default bitwise copy construction
        interpolateSurfaceElevation
        (
            const interpolateSurfaceElevation&
        ) = delete;


    //- Destructor
    virtual ~interpolateSurfaceElevation();


    // Member Functions

        virtual void evaluate();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const interpolateSurfaceElevation&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
