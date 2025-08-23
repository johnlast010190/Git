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
    (c) 2024 Engys Ltd.

Class
    Foam::LESdelta

Description
    Abstract base class for LES deltas

SourceFiles
    LESdelta.C

\*---------------------------------------------------------------------------*/

#ifndef LESdelta_H
#define LESdelta_H

#include "turbulenceModel.H"
#include "fields/volFields/volFields.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class LESdelta Declaration
\*---------------------------------------------------------------------------*/

class LESdelta
{
protected:

    // Protected data

        const turbulenceModel& turbulenceModel_;

        volScalarField delta_;


public:

    //- Runtime type information
    TypeName("LESdelta");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            LESdelta,
            dictionary,
            (
                const word& name,
                const turbulenceModel& turbulence,
                const dictionary& dict
            ),
            (name, turbulence, dict)
        );


    // Constructors

        //- Construct from name, turbulenceModel and dictionary
        LESdelta
        (
            const word& name,
            const turbulenceModel& turbulence
        );

        //- Disallow default bitwise copy construction
        LESdelta(const LESdelta&) = delete;


    // Selectors

        //- Return a reference to the selected LES delta
        static autoPtr<LESdelta> New
        (
            const word& name,
            const turbulenceModel& turbulence,
            const dictionary& dict
        );

        //- Return a reference to the selected LES delta
        static autoPtr<LESdelta> New
        (
            const word& name,
            const turbulenceModel& turbulence,
            const dictionary& dict,
            const dictionaryConstructorTable&
        );


    //- Destructor
    virtual ~LESdelta()
    {}


    // Member Functions

        //- Return turbulenceModel reference
        const turbulenceModel& turbulence() const
        {
            return turbulenceModel_;
        }

        //- Read the LESdelta dictionary
        virtual void read(const dictionary&) = 0;

        // Correct values
        virtual void correct() = 0;


    // Member Operators

        virtual operator const volScalarField&() const
        {
            return delta_;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const LESdelta&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
