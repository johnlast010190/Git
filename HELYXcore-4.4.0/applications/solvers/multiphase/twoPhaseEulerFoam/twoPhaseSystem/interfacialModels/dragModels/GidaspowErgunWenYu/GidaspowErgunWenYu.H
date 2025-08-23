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
    (c) 2011-2014 OpenFOAM Foundation

Class
    Foam::dragModels::GidaspowErgunWenYu

Description
    Gidaspow, Ergun, Wen and Yu drag model

    Reference:
    \verbatim
        "Multiphase flow and fluidization",
        Gidaspow, D.,
        Academic Press, New York, 1994.
    \endverbatim

SourceFiles
    GidaspowErgunWenYu.C

\*---------------------------------------------------------------------------*/

#ifndef GidaspowErgunWenYu_H
#define GidaspowErgunWenYu_H

#include "interfacialModels/dragModels/dragModel/dragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace dragModels
{

class Ergun;
class WenYu;

/*---------------------------------------------------------------------------*\
                           Class GidaspowErgunWenYu Declaration
\*---------------------------------------------------------------------------*/

class GidaspowErgunWenYu
:
    public dragModel
{
    // Private data

        //- Ergun drag model
        autoPtr<Ergun> Ergun_;

        //- Wen Yu drag model
        autoPtr<WenYu> WenYu_;


public:

    //- Runtime type information
    TypeName("GidaspowErgunWenYu");


    // Constructors

        //- Construct from a dictionary and an ordered phase pair
        GidaspowErgunWenYu
        (
            const dictionary& interfaceDict,
            const phasePair& pair,
            const bool registerObject
        );


    //- Destructor
    virtual ~GidaspowErgunWenYu();


    // Member Functions

        //- Drag coefficient
        virtual tmp<volScalarField> CdRe() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dragModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
