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
    Foam::dragModels::TomiyamaAnalytic

Description
    Analytical drag model of Tomiyama et al.

    Reference:
    \verbatim
        "Drag Coefficients of Bubbles. 1st Report. Drag Coefficients of a
         Single Bubble in a Stagnant Liquid."
        Tomiyama, A., Kataoka, I., and Sakaguchi, T.,
        Nippon Kikai Gakkai Ronbunshu
        Volume 61, Issue 587, 1995, pp. 2357-2364
    \endverbatim

SourceFiles
    TomiyamaAnalytic.C

\*---------------------------------------------------------------------------*/

#ifndef TomiyamaAnalytic_H
#define TomiyamaAnalytic_H

#include "interfacialModels/dragModels/dragModel/dragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace dragModels
{

/*---------------------------------------------------------------------------*\
                           Class TomiyamaAnalytic Declaration
\*---------------------------------------------------------------------------*/

class TomiyamaAnalytic
:
    public dragModel
{
    // Private data

        //- Residual Reynolds Number
        const dimensionedScalar residualRe_;

        //- Residual Eotvos number
        const dimensionedScalar residualEo_;

        //- Residual aspect ratio
        const dimensionedScalar residualE_;


public:

    //- Runtime type information
    TypeName("TomiyamaAnalytic");


    // Constructors

        //- Construct from a dictionary and a phase pair
        TomiyamaAnalytic
        (
            const dictionary& dict,
            const phasePair& pair,
            const bool registerObject
        );


    //- Destructor
    virtual ~TomiyamaAnalytic();


    // Member Functions

        // Drag coefficient
        virtual tmp<volScalarField> CdRe() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dragModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
