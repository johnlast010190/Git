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
    (c) 2015 OpenFOAM Foundation

Class
    Foam::eulerianDragModels::TomiyamaKataokaZunSakaguchi

Description
    Drag model for gas-liquid system of Tomiyama et al.

    Reference:
    \verbatim
        "Drag coefficients of single bubbles under normal and microgravity
         conditions"
        Tomiyama, A., Kataoka, I., Zun, I., Sakaguchi, T.
        JSME International Series B, Fluids and Thermal Engineering,
        Vol. 41, 1998, pp. 472-479
    \endverbatim

SourceFiles
    TomiyamaKataokaZunSakaguchi.C

\*---------------------------------------------------------------------------*/

#ifndef TomiyamaKataokaZunSakaguchi_H
#define TomiyamaKataokaZunSakaguchi_H

#include "../eulerianDragModel/eulerianDragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class eulerianPhasePair;

namespace eulerianDragModels
{

/*---------------------------------------------------------------------------*\
                           Class TomiyamaKataokaZunSakaguchi Declaration
\*---------------------------------------------------------------------------*/

class TomiyamaKataokaZunSakaguchi
:
    public eulerianDragModel
{
    // Private data

        //- Residual Reynolds Number
        const dimensionedScalar residualRe_;

        //- Residual Eotvos number
        const dimensionedScalar residualEo_;


public:

    //- Runtime type information
    TypeName("TomiyamaKataokaZunSakaguchi");


    // Constructors

        //- Construct from a dictionary and a phase pair
        TomiyamaKataokaZunSakaguchi
        (
            const dictionary& dict,
            const eulerianPhasePair& pair,
            const bool registerObject
        );


    //- Destructor
    virtual ~TomiyamaKataokaZunSakaguchi();


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
