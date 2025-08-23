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

Authors
    2019. Nikolaos Magoulas (Engys Ltd.). All rights reserved.

Class
    Foam::functionObjects::objectiveFunctionObject::
    pumpEfficiencyObjectiveFunctionObject

Description
    Calculates pump efficiency.

SourceFiles
    pumpEfficiencyObjectiveFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef pumpEfficiencyObjectiveFunctionObject_H
#define pumpEfficiencyObjectiveFunctionObject_H

#include "objectiveFunctionObject/objectiveFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
               Class pumpEfficiencyObjectiveFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class pumpEfficiencyObjectiveFunctionObject
:
    public objectiveFunctionObject
{
protected:

    // Protected data

        //- Switch on adiabatic behaviour for compressible flow (default false)
        Switch adiabatic_;

        autoPtr<scalar> pressureGain_;
        autoPtr<scalar> TOmega_;


    // Protected Member Functions

        void costDp() const;
        void costTw() const;


public:

    //- Runtime type information
    TypeName("pumpEfficiency");


    // Constructors

        //- Construct from Time and dictionary and adjoint file-format
        pumpEfficiencyObjectiveFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&,
            bool useAdjointFileFormat
        );

        //- Construct from Time and dictionary
        //  This constructor is needed so that the objective is added on the
        //  selection table of function objects.
        pumpEfficiencyObjectiveFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        pumpEfficiencyObjectiveFunctionObject
        (
            const pumpEfficiencyObjectiveFunctionObject&
        ) = delete;


    //- Destructor
    virtual ~pumpEfficiencyObjectiveFunctionObject();


    // Member Functions

        //- Read and update from dictionary
        virtual bool read(const dictionary&);

        //- Calculate the value of the cost function and write it to file
        virtual bool execute();

        //- Track and write
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const pumpEfficiencyObjectiveFunctionObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
