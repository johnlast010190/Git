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
    (c) 2011-2022 OpenFOAM Foundation
    (c) 2021-2025 Engys Ltd.

Class
    Foam::noChemistrySolver

Description
    Dummy chemistry solver for 'none' option

SourceFiles
    noChemistrySolver.H
    noChemistrySolver.C

\*---------------------------------------------------------------------------*/

#ifndef noChemistySolver_H
#define noChemistySolver_H

#include "chemistrySolver/chemistrySolver/chemistrySolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class noChemistrySolver Declaration
\*---------------------------------------------------------------------------*/

template<class ChemistryModel>
class noChemistrySolver
:
    public chemistrySolver<ChemistryModel>
{

public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from thermo
        noChemistrySolver(const fluidMulticomponentThermo& thermo);

        //- Construct from thermo
        noChemistrySolver(const solidMulticomponentThermo& thermo);


    //- Destructor
    virtual ~noChemistrySolver();


    // Member Functions

        //- Update the concentrations and return the chemical time
        virtual void solve
        (
            scalar& p,
            scalar& T,
            scalarField& c,
            const label li,
            scalar& deltaT,
            scalar& subDeltaT
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "chemistrySolver/noChemistrySolver/noChemistrySolver.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
