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
    (c) 2015-2017 OpenFOAM Foundation
    (c) 2021-2024 Engys Ltd

Class
    Foam::fv::mistThermalAndMassSource

Group
    grpFvOptionsSources

Description
    Calculates and applies the heat transfer between mist and humidity
    and heat transfer into the energy equation of the continuous phase.

Usage
    Example usage:
    mistThermalAndMassSource
    {
        type            mistThermalAndMassSource;
        active          true;
        mistThermalAndMassSourceCoeffs
        {
            mistProperties
            {
                // To be defined
                diam    0.001;         // droplets diameter
                terminalVelocity 0.0;

                // Fine tuning
                Dc    2e-05;           // if 0.0, it uses Dw
                Schmidt 0.7;
                liqDensity 1000;       // water density
                humSat 0.0;             //if 0.0, it is calculated in the code

                // Others
                mistName mist;
                humidityName w;
                pName p: // optional
                phiName phi; // optional
            }
            // Thermal equation field. This should be omitted for USF
            fields ( T );
        }
    }
    \endverbatim

SourceFiles
    mistThermalAndMassSource.C

\*---------------------------------------------------------------------------*/

#ifndef mistThermalAndMassSource_H
#define mistThermalAndMassSource_H

#include "solverObjects/solverObject/solverObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
               Class mistThermalAndMassSource Declaration
\*---------------------------------------------------------------------------*/

class mistThermalAndMassSource
:
    public solverObject
{
protected:

    // Protected data

        //- Mist phase name
        word mistName_;

        //- Humidity phase name
        word humidityName_;

        //- Energy equation field name
        word heTName_;

        //- Pressure field name
        word pName_;

        //- Phi name
        word phiName_;

        //- Molecular weight of vapour phase
        scalar Mvap_;

        //- Molecular weight of air
        scalar Mair_;

        dimensionedScalar Dc_;
        dimensionedScalar rhod_;
        dimensionedScalar diam_;
        dimensionedScalar dHevap_;
        dimensionedScalar Urel_;

        dimensionedScalar humSat_;
        dimensionedScalar Sc_;


public:

    //- Runtime type information
    TypeName("mistThermalAndMassSource");


    // Constructors

        //- Construct from explicit source name and mesh
        mistThermalAndMassSource
        (
            const word& sourceName,
            const objectRegistry& obr,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        mistThermalAndMassSource(const mistThermalAndMassSource&) = delete;

    //- Destructor
    virtual ~mistThermalAndMassSource()
    {}


    // Member Functions

        virtual void getSolveGraph
        (
            wordList& solveNames,
            HashTable<wordList>& requiredDependencies,
            HashTable<wordList>& optionalDependencies,
            HashTable<wordList>& correctorMembers
        ) override
        {}

        virtual void getSourceGraph
        (
            wordList& fieldNames,
            HashTable<wordList>& sourceDependencies
        ) override;

        // Evaluate

            //- Source term to mass and temperature equations
            virtual void addSup
            (
                fvMatrix<scalar>& eqn,
                const label fieldi
            ) override;

            //- Source term to energy equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<scalar>& eqn,
                const label fieldi
            ) override;

            //- Compute interfacial area density [1/m]
            tmp<volScalarField> iad() const;

            //- Compute the Sherwood number [-]
            tmp<volScalarField> Sh() const;

            //- Compute the driving concentration difference [kg vapour/kg air]
            tmp<volScalarField> deltaC() const;

            //- Compute mass transfer rate mDot [1/s]
            tmp<volScalarField> mDot() const;

            //- Return air density
            tmp<volScalarField> getRho() const;

            //- Return air Cp
            tmp<volScalarField> getCp() const;

            //- Return absolute pressure in Pa
            tmp<volScalarField> getPabs() const;

            //- Return saturation humidity
            tmp<volScalarField> humSat() const;

            //- Disperse phase Reynolds number
            tmp<volScalarField> Re() const;

            //- Return mass diffusion coeff
            tmp<volScalarField> massDCoeff() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const mistThermalAndMassSource&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
