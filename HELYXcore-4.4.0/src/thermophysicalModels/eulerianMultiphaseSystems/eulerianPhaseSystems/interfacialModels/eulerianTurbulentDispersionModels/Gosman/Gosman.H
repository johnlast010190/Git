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
    (c) 2014-2015 OpenFOAM Foundation

Class
    Foam::eulerianTurbulentDispersionModels::Gosman

Description
    Turbulent dispersion model of Gosman et al.

    Reference:
    \verbatim
        "Multidimensional modeling of turbulent two-phase flows in stirred
         vessels"
        Gosman, A.D., Lekakou, C., Politis, S., Issa, R.I., and Looney, M.K.,
        AIChE Journal
        Volume 38, Issue 12, 1992, pp. 1946-1956
     \endverbatim

SourceFiles
    Gosman.C

\*---------------------------------------------------------------------------*/

#ifndef Gosman_H
#define Gosman_H

#include "../eulerianTurbulentDispersionModel/eulerianTurbulentDispersionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class eulerianPhasePair;

namespace eulerianTurbulentDispersionModels
{

/*---------------------------------------------------------------------------*\
                             Class Gosman Declaration
\*---------------------------------------------------------------------------*/

class Gosman
:
    public eulerianTurbulentDispersionModel
{
    // Private data

        //- Schmidt number
        const dimensionedScalar sigma_;


public:

    //- Runtime type information
    TypeName("Gosman");


    // Constructors

        //- Construct from a dictionary and a phase pair
        Gosman
        (
            const dictionary& dict,
            const eulerianPhasePair& pair
        );


    //- Destructor
    virtual ~Gosman();


    // Member Functions

        //- Turbulent diffusivity
        //  multiplying the gradient of the phase-fraction
        virtual tmp<volScalarField> D() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace eulerianTurbulentDispersionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
