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
    (c) 2024 Engys Ltd.

Class
    Foam::coefficientMixture

Description
    The base class for coefficient mixtues. At the moment this is required for
    for making chemistry computations viable because polynomial computations
    are usually employed and by normal mising multiple polynmials are required
    to be computed. This is computationally expensive. However, by precomputing
    coefficients and combining them to one polynomial, the computation can be
    done in one step. So far it is done for the delta Gibbs free energy.

\*---------------------------------------------------------------------------*/

#ifndef coefficientMixture_H
#define coefficientMixture_H

#include "materialModels/materialTables/materialTables.H"
#include "fvMesh/fvMesh.H"
#include "fvSolutionRegistry/fvSolutionRegistry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class coefficientMixture Declaration
\*---------------------------------------------------------------------------*/

class coefficientMixture
{
public:

    //- Runtime type information
    TypeName("coefficientMixture");

    // Declare runtime construction
    declareRunTimeSelectionTable
    (
        autoPtr,
        coefficientMixture,
        dictionary,
        (
            const UPtrList<baseModels<scalar>>& components,
            const scalarField& coeffs,
            const word& name
        ),
        (components, coeffs, name)
    );


    // Selectors

        //- Return reference to selected rho model
        static autoPtr<coefficientMixture> New
        (
            const objectRegistry& obr,
            const scalarField& stoichCoeffs,
            const wordList& speciesNames,
            const word& phaseName,
            const word& name
        );


    // Constructors

        //- Default constructor
        coefficientMixture()
        {}

        //- Construct from components
        coefficientMixture
        (
            const UPtrList<baseModels<scalar>>& components,
            const scalarField& stoichCoeffs,
            const wordList& speciesNames,
            const word& name
        );


    //- Destructor
    virtual ~coefficientMixture();


    // Public member functions

        //- Return the Gibbs free energy delta
        virtual scalar gStd(const scalar T) const = 0;

        //- Delta of entropy
        virtual scalar deltaS(const scalar T) const = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
