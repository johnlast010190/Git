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
    (c) 2011-2022 OpenFOAM Foundation

Class
    Foam::laminarFlameSpeedModels::constant

Description
    Constant laminar flame speed model.

SourceFiles
    selectLaminarFlameSpeed.C

\*---------------------------------------------------------------------------*/

#ifndef constant_H
#define constant_H

#include "laminarFlameSpeed/laminarFlameSpeed.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarFlameSpeedModels
{

/*---------------------------------------------------------------------------*\
                           Class constant Declaration
\*---------------------------------------------------------------------------*/

class constant
:
    public laminarFlameSpeed
{

private:

    // Const laminar flame speed

        dimensionedScalar Su_;


    // Private Member Functions

        //- Construct as copy (not implemented)
        constant(const constant&);
        void operator=(const constant&);


public:

    //- Runtime type information
    TypeName("constant");


    // Constructors

        //- Construct from dictionary and psiuMulticomponentThermo
        constant
        (
            const dictionary&,
            const psiuMulticomponentThermo&
        );


    //- Destructor
    virtual ~constant();


    // Member functions

        //- Return the laminar flame speed [m/s]
        tmp<volScalarField> operator()() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End laminarFlameSpeedModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
