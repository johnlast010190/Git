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
    (c) 2011-2024 Engys Ltd.

Class
    Foam::blendingFunction

Description
    Abstract base class for run-time selectable blending functions

SourceFiles
    blendingFunction.C

\*---------------------------------------------------------------------------*/

#ifndef blendingFunction_H
#define blendingFunction_H

#include "db/typeInfo/typeInfo.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class blendingFunction Declaration
\*---------------------------------------------------------------------------*/

class blendingFunction
{
public:

    //- Runtime type information
    TypeName("blendingFunction");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        blendingFunction,
        dictionary,
        (
            const dictionary& blendingFunctionDict,
            bool dum
        ),
        (blendingFunctionDict, dum)
    );


    // Constructors

        //- Construct from components
        blendingFunction(const dictionary& data, bool dum = true);

        //- Disallow default bitwise copy construction
        blendingFunction(const blendingFunction&) = delete;


    // Selectors

        //- Return an autoPtr to the selected blendingFunction
        static autoPtr<blendingFunction> New
        (
            const dictionary& blendingFunctionDict
        );


    //- Destructor
    virtual ~blendingFunction()
    {}


    // Member Functions

        //- Update blending function
        virtual void update(scalar) = 0;


    // Member Operators

        //- convert value into blending factor
        virtual scalar operator()(scalar) const = 0;

        //- Disallow default bitwise assignment
        void operator=(const blendingFunction&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
