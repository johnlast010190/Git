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
    (c) 2021 Engys Ltd.

Class
    Foam::EulerImplicit

Description
    An Euler implicit solver for chemistry

    Euler implicit integration based on the reaction rate Jacobian is used to
    solver for the composition and Euler explicit integration for the
    temperature.  Separating the integration of temperature from composition
    in this manner is significantly more stable for exothermic systems

SourceFiles
    EulerImplicit.C

\*---------------------------------------------------------------------------*/

#ifndef EulerImplicit_H
#define EulerImplicit_H

#include "chemistrySolver/chemistrySolver/chemistrySolver.H"
#include "matrices/simpleMatrix/simpleMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class EulerImplicit Declaration
\*---------------------------------------------------------------------------*/

template<class ChemistryModel>
class EulerImplicit
:
    public chemistrySolver<ChemistryModel>
{
    // Private data

        //- Coefficients dictionary
        dictionary coeffsDict_;

        //- Chemistry timescale coefficient
        scalar cTauChem_;

        //- Field encapsulating the composition, temperature and pressure
        mutable scalarField cTp_;

        //- Reaction rate field
        mutable scalarField R_;

        //- Reaction Jacobian
        mutable scalarSquareMatrix J_;

        //- Euler implicit integration matrix for composition
        mutable simpleMatrix<scalar> E_;


public:

    //- Runtime type information
    TypeName("EulerImplicit");


    // Constructors

        //- Construct from thermo
        EulerImplicit(const fluidMulticomponentThermo& thermo);


    //- Destructor
    virtual ~EulerImplicit();


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
    #include "chemistrySolver/EulerImplicit/EulerImplicit.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
