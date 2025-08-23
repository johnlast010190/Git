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
    (c) 2012-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#ifndef coalCloudList_H
#define coalCloudList_H

#include "coalCloud/coalCloud.H"
#include "fields/volFields/volFieldsFwd.H"
#include "fvMatrices/fvMatricesFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class coalCloudList Declaration
\*---------------------------------------------------------------------------*/

class coalCloudList
:
    public PtrList<coalCloud>
{
private:

    //- Reference to the mesh
    const fvMesh& mesh_;


public:

    // Constructor

        coalCloudList
        (
            const volScalarField& rho,
            const volVectorField& U,
            const dimensionedVector& g,
            const SLGThermo& slgThermo
        );


    // Member Functions

        // Evolution

            //- Evolve the cloud collection
            void evolve();


        // Source terms

            //- Return const reference to momentum source
            inline tmp<volVectorField::Internal> UTrans() const;

            //- Return tmp momentum source term
            inline tmp<fvVectorMatrix> SU(volVectorField& U) const;

            //- Sensible enthalpy transfer [J/kg]
            inline tmp<volScalarField::Internal> hsTrans() const;

            //- Return sensible enthalpy source term [J/kg/m3/s]
            inline tmp<fvScalarMatrix> Sh(volScalarField& hs) const;


            //- Return mass source term for specie i - specie eqn
            inline tmp<fvScalarMatrix> SYi
            (
                const label i,
                volScalarField& Yi
            ) const;

            //- Return total mass transfer [kg/m3]
            inline tmp<volScalarField::Internal> rhoTrans() const;

            //- Return tmp total mass source for carrier phase
            //  - fully explicit
            inline tmp<volScalarField::Internal> Srho() const;

            //- Return tmp total mass source for carrier phase specie i
            //  - fully explicit
            inline tmp<volScalarField::Internal> Srho
            (
                const label i
            ) const;

            //- Return total mass source term [kg/m3/s]
            inline tmp<fvScalarMatrix> Srho(volScalarField& rho) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "coalCloudList/coalCloudListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
