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
    Foam::eulerianAspectRatioModels::Wellek

Description
    Aspect ratio model of Wellek et al.

    References:
    \verbatim
        "Implementation and Comparison of Correlations for interfacial Forces
         in a Gas-Liquid System within an Euler-Euler Framework"
        Otromke, M.,
        PhD Thesis, April 2013
    \endverbatim

    \verbatim
        "Shape of liquid drops moving in liquid media"
        Wellek, R.M., Agrawal, A.K., Skelland, A.H.P.,
        International Journal of Multiphase Flow
        Volume 12, Issue 5, September 1966, pp. 854-862
    \endverbatim

SourceFiles
    Wellek.C

\*---------------------------------------------------------------------------*/

#ifndef Wellek_H
#define Wellek_H

#include "../eulerianAspectRatioModel/eulerianAspectRatioModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace eulerianAspectRatioModels
{

/*---------------------------------------------------------------------------*\
                       Class Wellek Declaration
\*---------------------------------------------------------------------------*/

class Wellek
:
    public eulerianAspectRatioModel
{
public:

    //- Runtime type information
    TypeName("Wellek");


    // Constructors

        //- Construct from a dictionary and an ordered phase pair
        Wellek
        (
            const dictionary& dict,
            const eulerianPhasePair& pair
        );


    //- Destructor
    virtual ~Wellek();


    // Member Functions

        //- Aspect ratio
        virtual tmp<volScalarField> E() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace eulerianAspectRatioModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
