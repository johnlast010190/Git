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
    (c) 2014 OpenFOAM Foundation

Class
    Foam::immiscibleIncompressibleThreePhaseDriftMixture

Description
    An immiscible incompressible two-phase mixture transport model

SourceFiles
    immiscibleIncompressibleThreePhaseDriftMixture.C

\*---------------------------------------------------------------------------*/

#ifndef immiscibleIncompressibleThreePhaseDriftMixture_H
#define immiscibleIncompressibleThreePhaseDriftMixture_H

#include "incompressibleThreePhaseDriftMixture.H"
#include "threePhaseDriftInterfaceProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class immiscibleIncompressibleThreePhaseDriftMixture Declaration
\*---------------------------------------------------------------------------*/

class immiscibleIncompressibleThreePhaseDriftMixture
:
    public incompressibleThreePhaseDriftMixture,
    public threePhaseDriftInterfaceProperties
{

public:

    // Constructors

        //- Construct from components
        immiscibleIncompressibleThreePhaseDriftMixture
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~immiscibleIncompressibleThreePhaseDriftMixture()
    {}


    // Member Functions

        //- Correct the transport and interface properties
        virtual void correct()
        {
            incompressibleThreePhaseDriftMixture::correct();
            threePhaseDriftInterfaceProperties::correct();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
