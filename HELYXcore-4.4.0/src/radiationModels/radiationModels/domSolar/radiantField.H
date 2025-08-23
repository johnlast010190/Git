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
    (c) 2019-2024 Engys Ltd.

Class
    Foam::radiantField

Description
    A storage class for basic radiant field.

    Example of the object specification:
    @verbatim
    radiantField
    {
        radiantInputType userDefined; //solarCalculator

        intensity constant 10; //radiation intensity for userDefined
        direction constant (0 -1 0); //radiation direction for userDefined

        //for solarCalculator input see
        // radiation/solarCalculator/solarCalculator.H

    }
    @endverbatim

SourceFiles
    radiantField.C

\*---------------------------------------------------------------------------*/

#ifndef radiantField_H
#define radiantField_H

#include "fields/volFields/volFields.H"
#include "primitives/functions/Function1/Function1/Function1.H"
#include "solarCalculator/solarCalculator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class radiantField Declaration
\*---------------------------------------------------------------------------*/

class radiantField
:
    public volScalarField
{
public:

        enum radiantInputType
        {
            ritUserDefined,
            ritSolarCalculator
        };

        static const NamedEnum<radiantInputType, 2> radiantInputTypeNames_;


private:

        //- Direction/intensity input mode - valid entries are:
        //  userDefined, solarCalculator
        radiantInputType inputMode_;

        // For now just store both input options as autoPtrs
        // run time selection should be added if options expand

            //- Radiation intensity
            autoPtr<Function1<scalar>> radiationFlux_;

            //- Radiation direction vector
            autoPtr<Function1<vector>> radiationDirection_;

            //- Solar calculater
            autoPtr<solarCalculator> solarCalc_;

        //- Patch source coeffs
        const FieldField<Field, scalar>& patchSourceCoeffs_;

        //- Current solar flux
        scalar currentFlux_;

        //- Current solar direction
        vector currentDirection_;

        //- Update solar tracking index
        label updateTimeIndex_;

        //- Initialisation switch
        bool updateAfterRead_;

        //- Initialisation switch
        bool initializeField_;

        //- Initial tolerance for explicit wave solve
        scalar initTol_;

        //- Initialisation matrix solver switch
        label nSolverSolutions_;

        //- Matrix solver switch
        bool solveField_;

        //- Inlet outlet scaling
        scalar qScale_;


    // Private functions

        //- set fv dicts if matrix solvers to be used
        void setFvDicts();

        //- Provide boundary conditions for intensity field
        const dictionary radiantBoundaryConditions
        (
            const fvMesh& mesh,
            word
        ) const;

        //- create initial intensity field
        tmp<volScalarField> autoCreateRadiantField
        (
            const word& fieldName,
            const fvMesh& mesh
        );

        //- Calculate scaling for incedent radiation to match incoming
        //  radiation transmission
        scalar inletOutletScaling();

        //- Update radiant power flux and direction
        void updateRadiationProperties();

        //- Add diffuse radiation over patches
        //  based on ASHRAE hanbook
        void updateSkyDiffusiveRadiation(volScalarField&);

        //- Return vertical direction
        vector verticalDirection();

        //- Initialise radiant field and input functions
        void initialise(const dictionary&);


public:

    // Constructors

        //- Construct for given objectRegistry and dictionary.
        radiantField
        (
            const word& fieldName,
            const fvMesh&,
            const dictionary&,
            const FieldField<Field, scalar>& transmissivity
        );

        //- Disallow default bitwise copy construction
        radiantField(const radiantField&) = delete;


    //- Destructor
    virtual ~radiantField();


    // Member Functions

        //- Update solver
        void calculate(volScalarField&);

        scalar getZenithAngle() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const radiantField&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
