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
    (c) 2024-2025 Engys Ltd.

Class
    Foam::basicChemistryModel

Description
    Base class for chemistry models

SourceFiles
    basicChemistryModelI.H
    basicChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef basicChemistryModel_H
#define basicChemistryModel_H

#include "fluidMulticomponentThermo/fluidMulticomponentThermo.H"
#include "solidMulticomponentThermo/solidMulticomponentThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                     class basicChemistryModel Declaration
\*---------------------------------------------------------------------------*/

class basicChemistryModel
:
    public IOdictionary
{
public:

    // Public Enumerations

        //- Enumeration for the type of Jacobian to be calculated by the
        //  chemistry model
        enum class jacobianType
        {
            fast,
            exact
        };

        //- Jacobian type names
        static const NamedEnum<jacobianType, 2> jacobianTypeNames_;


protected:

    // Protected data

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Reference to the local registry
        const objectRegistry& obr_;

        //- Reference to the thermo
        const fluidMulticomponentThermo* thermo_;

        //- Reference to the thermo
        const solidMulticomponentThermo* solidThermo_;

        //- Chemistry activation switch
        Switch chemistry_;

        //- Initial chemical time step
        const scalar deltaTChemIni_;

        //- Maximum chemical time step
        const scalar deltaTChemMax_;

        //- Latest estimation of integration step
        volScalarField::Internal deltaTChem_;


public:

    //- Runtime type information
    TypeName("basicChemistryModel");

    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        basicChemistryModel,
        thermo,
        (const fluidMulticomponentThermo& thermo),
        (thermo)
    );


    // Constructors

        //- Construct from thermo
        basicChemistryModel(const fluidMulticomponentThermo& thermo);

        //- Construct from solid thermo
        basicChemistryModel(const solidMulticomponentThermo& thermo);

        //- Disallow default bitwise copy construction
        basicChemistryModel(const basicChemistryModel&) = delete;


    // Selectors

        //- Select based on fluid reaction thermo
        static autoPtr<basicChemistryModel> New
        (
            const fluidMulticomponentThermo& thermo
        );


    //- Destructor
    virtual ~basicChemistryModel();


    // Member Functions

        //- Return const access to the mesh database
        inline const fvMesh& mesh() const;

        //- Return object registry
        inline const objectRegistry& obr() const;

        //- Return const access to the thermo
        inline const fluidMulticomponentThermo& thermo() const;

        //- Return const access to the thermo
        inline const solidMulticomponentThermo& solidThermo() const;

        //- Chemistry activation switch
        inline Switch chemistry() const;

        //- The number of species
        virtual label nSpecie() const = 0;

        //- The number of reactions
        virtual label nReaction() const = 0;

        //- Return the latest estimation of integration step
        inline const volScalarField::Internal& deltaTChem() const;

        //- Return reaction rates of the species [kg/m^3/s]
        virtual const PtrList<volScalarField::Internal>& RR() const = 0;

        //- Return reaction rates of the species in reactioni [kg/m^3/s]
        virtual PtrList<volScalarField::Internal> reactionRR
        (
            const label reactioni
        ) const = 0;

        //- Return heat release of the reactioni
        virtual tmp<volScalarField::Internal> reactionQdot
        (
            const label reactioni
        ) const = 0;

        //- Calculates the reaction rates
        virtual void calculate() = 0;

        //- Solve the reaction system for the given time step
        //  and return the characteristic time
        virtual scalar solve(const scalar deltaT) = 0;

        //- Solve the reaction system for the given time step
        //  and return the characteristic time
        virtual scalar solve(const scalarField& deltaT) = 0;

        //- Return the chemical time scale
        virtual tmp<volScalarField> tc() const = 0;

        //- Return the heat release rate [kg/m/s^3]
        virtual tmp<volScalarField> Qdot() const = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const basicChemistryModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "basicChemistryModel/basicChemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
