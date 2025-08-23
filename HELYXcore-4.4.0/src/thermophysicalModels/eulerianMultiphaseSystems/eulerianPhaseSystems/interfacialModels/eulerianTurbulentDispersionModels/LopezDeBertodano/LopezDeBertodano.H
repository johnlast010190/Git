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
    Foam::eulerianTurbulentDispersionModels::LopezDeBertodano

Description
    Lopez de Bertodano (1992) turbulent dispersion model.

    \verbatim
        "Turbulent bubbly two-phase flow in a triangular
         duct"
        Lopez de Bertodano, M.
        Ph.D. Thesis, Rensselaer Polytechnic Institution, New York, USA, 1992.
    \endverbatim

    \verbatim
        "The Favre averaged drag model for turbulent dispersion in Eulerian
         multi-phase flows"
        Burns, A.D., Frank, T., Hamill, I., Shi, J.M.,
        5th international conference on multiphase flow
        Volume 4, Paper 392, May 2004
    \endverbatim

SourceFiles
    LopezDeBertodano.C

\*---------------------------------------------------------------------------*/

#ifndef LopezDeBertodano_H
#define LopezDeBertodano_H

#include "..//eulerianTurbulentDispersionModel/eulerianTurbulentDispersionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class eulerianPhasePair;

namespace eulerianTurbulentDispersionModels
{

/*---------------------------------------------------------------------------*\
              Class LopezDeBertodano Declaration
\*---------------------------------------------------------------------------*/

class LopezDeBertodano
:
    public eulerianTurbulentDispersionModel
{
    // Private data

        //- Constant turbulent dispersion coefficient
        const dimensionedScalar Ctd_;


public:

    //- Runtime type information
    TypeName("LopezDeBertodano");


    // Constructors

        //- Construct from a dictionary and a phase pair
        LopezDeBertodano
        (
            const dictionary& dict,
            const eulerianPhasePair& pair
        );


    //- Destructor
    virtual ~LopezDeBertodano();


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
