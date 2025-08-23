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
    Foam::eulerianLiftModels::Moraga

Description
    Lift model of Moraga et al.

    References:
    \verbatim
        "Implementation and Comparison of Correlations for interfacial Forces
         in a Gas-Liquid System within an Euler-Euler Framework"
        Otromke, M.,
        PhD Thesis, April 2013
    \endverbatim

    \verbatim
        "Lateral forces on spheres in turbulent uniform shear flow"
        Moraga, F.J., Bonetto, F.J., Lahey, R.T.,
        International Journal of Multiphase Flow
        Volume 25, Issues 6-7, September 1999, pp. 1321-1372
    \endverbatim

SourceFiles
    Moraga.C

\*---------------------------------------------------------------------------*/

#ifndef Moraga_H
#define Moraga_H

#include "../eulerianLiftModel/eulerianLiftModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class eulerianPhasePair;

namespace eulerianLiftModels
{

/*---------------------------------------------------------------------------*\
                           Class Moraga Declaration
\*---------------------------------------------------------------------------*/

class Moraga
:
    public eulerianLiftModel
{
public:

    //- Runtime type information
    TypeName("Moraga");


    // Constructors

        //- Construct from a dictionary and a phase pair
        Moraga
        (
            const dictionary& dict,
            const eulerianPhasePair& pair
        );


    //- Destructor
    virtual ~Moraga();


    // Member Functions

        //- Lift coefficient
        virtual tmp<volScalarField> Cl() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace eulerianLiftModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
