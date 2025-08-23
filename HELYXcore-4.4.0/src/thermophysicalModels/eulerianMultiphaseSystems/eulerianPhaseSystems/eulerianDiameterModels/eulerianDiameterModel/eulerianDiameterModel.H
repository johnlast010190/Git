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
    (c) 2011-2015 OpenFOAM Foundation

Class
    Foam::eulerianDiameterModel

Description
    A2stract base-class for dispersed-phase particle diameter models.

SourceFiles
    eulerianDiameterModel.C
    newDiameterModel.C

\*---------------------------------------------------------------------------*/

#ifndef diameterModel_H
#define diameterModel_H

#include "db/dictionary/dictionary.H"
#include "../../eulerianPhaseModel/eulerianPhaseModel/eulerianPhaseModel.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class eulerianDiameterModel Declaration
\*---------------------------------------------------------------------------*/

class eulerianDiameterModel
{

protected:

    // Protected data

        dictionary diameterProperties_;
        const eulerianPhaseModel& phase_;


public:

    //- Runtime type information
    TypeName("diameterModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            eulerianDiameterModel,
            dictionary,
            (
                const dictionary& diameterProperties,
                const eulerianPhaseModel& phase
            ),
            (diameterProperties, phase)
        );


    // Constructors

        eulerianDiameterModel
        (
            const dictionary& diameterProperties,
            const eulerianPhaseModel& phase
        );


    //- Destructor
    virtual ~eulerianDiameterModel();


    // Selectors

        static autoPtr<eulerianDiameterModel> New
        (
            const dictionary& diameterProperties,
            const eulerianPhaseModel& phase
        );


    // Member Functions

        //- Return the phase diameter properties dictionary
        const dictionary& diameterProperties() const
        {
            return diameterProperties_;
        }

        //- Return the phase
        const eulerianPhaseModel& phase() const
        {
            return phase_;
        }

        //- Return the phase mean diameter field
        virtual tmp<volScalarField> d() const = 0;

        //- Correct the diameter field
        virtual void correct();

        //- Read phaseProperties dictionary
        virtual bool read(const dictionary& phaseProperties) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
