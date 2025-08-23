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
    (c) 2013-2022 OpenFOAM Foundation
    (c) 2024-2025 Engys Ltd.

Class
    Foam::basicSolidChemistryModel

Description
    Chemistry model for solid thermodynamics

SourceFiles
    basicSolidChemistryModelI.H
    basicSolidChemistryModel.C
    newChemistrySolidModel.C

\*---------------------------------------------------------------------------*/

#ifndef basicSolidChemistryModel_H
#define basicSolidChemistryModel_H

#include "basicChemistryModel/basicChemistryModel.H"
#include "memory/autoPtr/autoPtr.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "solidMulticomponentThermo/solidMulticomponentThermo.H"
#include "speciesTable/speciesTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                     class basicSolidChemistryModel Declaration
\*---------------------------------------------------------------------------*/

class basicSolidChemistryModel
:
    public basicChemistryModel
{
protected:

    // Protected data

        //- Solid thermo
        const solidMulticomponentThermo& solidThermo_;


public:

    //- Runtime type information
    TypeName("basicSolidChemistryModel");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        basicSolidChemistryModel,
        thermo,
        (const solidMulticomponentThermo& thermo),
        (thermo)
    );


    // Constructors

        //- Construct from thermo
        basicSolidChemistryModel(const solidMulticomponentThermo& thermo);

        //- Disallow default bitwise copy construction
        basicSolidChemistryModel(const basicSolidChemistryModel&) = delete;


    //- Selector
    static autoPtr<basicSolidChemistryModel> New
    (
        const solidMulticomponentThermo& thermo
    );


    //- Destructor
    virtual ~basicSolidChemistryModel();


    // Member Functions

        //- Return const access to the solid thermo package
        inline const solidMulticomponentThermo& solidThermo() const
        {
            return solidThermo_;
        }

        //- Return total gases mass source term [kg/m3/s]
        virtual tmp<volScalarField::Internal> RRg() const = 0;

        //- Return total solids mass source term [kg/m3/s]
        virtual tmp<volScalarField::Internal> RRs() const = 0;

        //- Return chemical source terms for solids [kg/m3/s]
        virtual const volScalarField::Internal& RRs(const label i) const = 0;

        //- Return chemical source terms for gases [kg/m3/s]
        virtual const volScalarField::Internal& RRg(const label i) const = 0;

        //- Return net solid sensible enthalpy [J/Kg]
        virtual tmp<DimensionedField<scalar, volMesh>> RRsHs() const = 0;

        //- Set reacting status of cell, celli
        virtual void setCellReacting(const label celli, const bool active) = 0;

        //- Calculates the reaction rates
        virtual void calculate() = 0;

        //- Returns the reaction rate of the reactioni
        virtual PtrList<volScalarField::Internal> reactionRR
        (
            const label reactioni
        ) const
        {
            NotImplemented;
        }

        //- Returns the reaction rate of the reactioni
        virtual tmp<volScalarField::Internal> reactionQdot
        (
            const label reactioni
        ) const
        {
            NotImplemented;
        }

        //- Return reaction rates of the species [kg/m^3/s]
        virtual const PtrList<volScalarField::Internal>& RR() const
        {
            NotImplemented;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const basicSolidChemistryModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
