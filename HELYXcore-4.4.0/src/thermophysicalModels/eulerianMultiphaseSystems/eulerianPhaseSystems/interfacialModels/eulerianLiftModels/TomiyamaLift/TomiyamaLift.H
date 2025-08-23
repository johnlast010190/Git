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
    Foam::eulerianLiftModels::TomiyamaLift

Description
    Lift model of Tomiyama et al.

    Reference:
    \verbatim
        "Transverse migration of single bubbles in simple shear flows"
        Tomiyama, A., Tamai, H., Zun, I., Hosokawa, S.,
        Chemical Engineering Science
        Volume 57, Issue 11, June 2002, pp. 1849-1858
    \endverbatim

    The coefficient for pow3(EoH) proposed by Tomiyama (2002) has been modified
    to make the model continuous at EoH = 10.7 while maintaining the
    lift coefficient proposed by Tomiyama (2002) when EoH >= 10.7.

SourceFiles
    TomiyamaLift.C

\*---------------------------------------------------------------------------*/

#ifndef TomiyamaLift_H
#define TomiyamaLift_H

#include "../eulerianLiftModel/eulerianLiftModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class eulerianPhasePair;

namespace eulerianLiftModels
{

/*---------------------------------------------------------------------------*\
                           Class TomiyamaLift Declaration
\*---------------------------------------------------------------------------*/

class TomiyamaLift
:
    public eulerianLiftModel
{
public:

    //- Runtime type information
    TypeName("Tomiyama");


    // Constructors

        //- Construct from a dictionary and a phase pair
        TomiyamaLift
        (
            const dictionary& dict,
            const eulerianPhasePair& pair
        );


    //- Destructor
    virtual ~TomiyamaLift();


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
