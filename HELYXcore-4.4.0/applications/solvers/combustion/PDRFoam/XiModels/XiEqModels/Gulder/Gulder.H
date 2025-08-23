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
    (c) 2011-2012 OpenFOAM Foundation
    (c) 2024 Engys Ltd

Class
    Foam::XiEqModels::Gulder

Description
    Simple Gulder model for XiEq based on Gulders correlation
    with a linear correction function to give a plausible profile for XiEq.

SourceFiles
    Gulder.C

\*---------------------------------------------------------------------------*/

#ifndef Gulder_H
#define Gulder_H

#include "XiEqModel/XiEqModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace XiEqModels
{

/*---------------------------------------------------------------------------*\
                          Class Gulder Declaration
\*---------------------------------------------------------------------------*/

class Gulder
:
    public XiEqModel
{
    // Private data

        //- Model constant
        scalar XiEqCoef_;

        //- Minimum laminar burning velocity
        const dimensionedScalar SuMin_;

        //- Schelkin effect Model constant
        scalar uPrimeCoef_;

        //- Use sub-grid Schelkin effect
        bool subGridSchelkin_;


public:

    //- Runtime type information
    TypeName("Gulder");


    // Constructors

        //- Construct from components
        Gulder
        (
            const dictionary& XiEqProperties,
            const psiuMulticomponentThermo& thermo,
            const compressible::RASModel& turbulence,
            const volScalarField& Su
        );

        //- Disallow default bitwise copy construction
        Gulder(const Gulder&) = delete;


    //- Destructor
    virtual ~Gulder();


    // Member Functions

        //- Return the flame-wrinking XiEq
        virtual tmp<volScalarField> XiEq() const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& XiEqProperties);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const Gulder&) = delete;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace XiEqModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
