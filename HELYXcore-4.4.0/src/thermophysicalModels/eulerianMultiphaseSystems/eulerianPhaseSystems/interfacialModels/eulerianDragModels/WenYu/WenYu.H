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
    (c) 2011-2015 OpenFOAM Foundation

Class
    Foam::eulerianDragModels::WenYu

Description
    Wen and Yu drag model

    Reference:
    \verbatim
        "Eulerian Two-Phase Flow Theory Applied to Fluidization"
        Enwald, H., Peirano, E., Almstedt, A-E.,
        Int. J. Multiphase Flow, Vol. 22, Suppl, 1996, pp. 21-66
        Eq. 86-87, p. 40
    \endverbatim

SourceFiles
    WenYu.C

\*---------------------------------------------------------------------------*/

#ifndef WenYu_H
#define WenYu_H

#include "../eulerianDragModel/eulerianDragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class eulerianPhasePair;

namespace eulerianDragModels
{

/*---------------------------------------------------------------------------*\
                           Class WenYu Declaration
\*---------------------------------------------------------------------------*/

class WenYu
:
    public eulerianDragModel
{
    // Private data

        //- Residual Reynolds Number
        const dimensionedScalar residualRe_;


public:

    //- Runtime type information
    TypeName("WenYu");


    // Constructors

        //- Construct from a dictionary and a phase pair
        WenYu
        (
            const dictionary& dict,
            const eulerianPhasePair& pair,
            const bool registerObject
        );


    //- Destructor
    virtual ~WenYu();


    // Member Functions

        //- Drag coefficient
        virtual tmp<volScalarField> CdRe() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace eulerianDragModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
