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
    (c) 2019 OpenFOAM Foundation

Class
    Foam::functionObjects::Qdot

Description
    Calculates and outputs the heat release rate for the current combustion
    model.

SourceFiles
    Qdot.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_Qdot_H
#define functionObjects_Qdot_H

#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class Qdot Declaration
\*---------------------------------------------------------------------------*/

class Qdot
:
    public fvMeshFunctionObject
{
private:

    // Private data

        //- The name of the phase
        word phaseName_;

        //- Heat realease rate
        autoPtr<volScalarField> Qdot_;


public:

    //- Runtime type information
    TypeName("Qdot");


    // Constructors

        //- Construct from Time and dictionary
        Qdot
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        Qdot(const Qdot&) = delete;


    //- Destructor
    virtual ~Qdot();


    // Member Functions

        //- Read the data
        virtual bool read(const dictionary&);

        //- Calculate the Qdot field
        virtual bool execute();

        //- Do nothing
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const Qdot&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //