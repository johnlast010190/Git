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
    (c) 2011-2019 OpenFOAM Foundation

Class
    Foam::regionModels::thermalBaffle::thermalBaffle

Description
    2D thermal baffle

SourceFiles
    thermalBaffle.C
    thermalBaffleI.H

\*---------------------------------------------------------------------------*/

#ifndef thermalBaffle_H
#define thermalBaffle_H

#include "fields/volFields/volFieldsFwd.H"
#include "regionModel/regionModel/regionModel.H"
#include "solidThermo/solidThermo.H"
#include "radiationModels/radiationModel/radiationModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{

/*---------------------------------------------------------------------------*\
                        Class thermalBaffle Declaration
\*---------------------------------------------------------------------------*/

class thermalBaffle
:
    public regionModel
{
    // Private Member Functions

        //- Calculate the baffle mesh thickness
        tmp<volScalarField::Internal> calcDelta() const;


protected:

    // Protected data

        // Geometric parameters

            //- Baffle mesh thickness
            volScalarField::Internal delta_;

            //- Baffle physical thickness
            volScalarField::Internal thickness_;


        // Solution parameters

            //- Number of non orthogonal correctors
            label nNonOrthCorr_;


        // Thermo properties

            //- Solid thermo
            autoPtr<solidThermo> thermo_;

            //- Enthalpy/internal energy
            volScalarField& he_;


        // Source term fields

            //- Surface energy source [J/m2/s]
            volScalarField Qs_;

            //- Volumetric energy source [J/m3/s]
            volScalarField Q_;


        // Sub models

            //- Pointer to radiation model
            autoPtr<radiationModel> radiation_;


    // Protected member functions

        // IO

            //- Read control parameters IO dictionary
            virtual bool read();

            //- Read control parameters from dictionary
            virtual bool read(const dictionary& dict);


        // Equations

            //- Solve energy equation
            void solveEnergy();


public:

    //- Runtime type information
    TypeName("thermalBaffle");


    // Constructors

        //- Construct from components and dict
        thermalBaffle
        (
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        thermalBaffle(const thermalBaffle&) = delete;


    //- Destructor
    virtual ~thermalBaffle();


    // Member Functions

        // Thermo properties

            //- Return const reference to the solidThermo
            const solidThermo& thermo() const;


            // Fields

                //- Return the film specific heat capacity [J/kg/K]
                const volScalarField& Cp() const;

                //- Return solid absorptivity [1/m]
                const tmp<volScalarField> kappaRad() const;

                //- Return temperature [K]
                const volScalarField& T() const;

                //- Return density [Kg/m^3]
                const volScalarField& rho() const;

                //- Return thermal conductivity [W/m/K]
                const tmp<volScalarField> kappa() const;


        // Evolution

            //- Pre-evolve  thermal baffle
            virtual void preEvolveRegion();

            //- Evolve the thermal baffle
            virtual void evolveRegion();


       // I-O

            //- Provide some feedback
            virtual void info();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const thermalBaffle&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //