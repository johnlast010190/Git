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
    volumeUniformityObjectiveFunctionObject

Description
    Calculates the flow uniformity inside a volume.

SourceFiles
    volumeUniformityObjectiveFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef volumeUniformityObjectiveFunctionObject_H
#define volumeUniformityObjectiveFunctionObject_H

#include "objectiveFunctionObject/objectiveFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
               Class volumeUniformityObjectiveFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class volumeUniformityObjectiveFunctionObject
:
    public objectiveFunctionObject
{
protected:

    // Protected data

        //- Uniformity in the secified direction
        autoPtr<vector> desiredFlowDirection_;

        //- Specified target velocity
        autoPtr<vector> targetVelocity_;

        //- Objective multiplier
        scalar sourceCoeff_;

        //- Zone volume
        scalar volume_;


    // Protected Member Functions

        scalar zoneVolume() const;

        vector targetVelocity() const;


public:

    //- Runtime type information
    TypeName("volumeUniformity");


    // Constructors

        //- Construct from Time and dictionary and adjoint file-format
        volumeUniformityObjectiveFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&,
            bool useAdjointFileFormat
        );

        //- Construct from Time and dictionary
        //  This constructor is needed so that the objective is added on the
        //  selection table of function objects.
        volumeUniformityObjectiveFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        volumeUniformityObjectiveFunctionObject
        (
            const volumeUniformityObjectiveFunctionObject&
        ) = delete;


    //- Destructor
    virtual ~volumeUniformityObjectiveFunctionObject();


    // Member Functions

        //- Read and update from dictionary
        virtual bool read(const dictionary&);

        //- Calculate the value of the cost function and write it to file
        virtual bool execute();

        //- Track and write
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=
        (
            const volumeUniformityObjectiveFunctionObject&
        ) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
