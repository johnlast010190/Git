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
    rearLiftObjectiveFunctionObject

Description
    Calculates the value of rear lift.

SourceFiles
    rearLiftObjectiveFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef rearLiftObjectiveFunctionObject_H
#define rearLiftObjectiveFunctionObject_H

#include "objectiveFunctionObject/objectiveFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
               Class rearLiftObjectiveFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class rearLiftObjectiveFunctionObject
:
    public objectiveFunctionObject
{
protected:

    // Protected data

        //- Direction of lift
        vector forceDirection_;

        //- Direction of moment contributing to the objective
        vector momentDirection_;

        //- Rotation centre
        //  Usually placed on the ground, in the middle of the front and
        //  rear axis.
        vector rotationCentre_;

        //- Distance between the front and rear axis
        scalar axisLength_;


public:

    //- Runtime type information
    TypeName("rearLift");


    // Constructors

        //- Construct from Time and dictionary and adjoint file-format
        rearLiftObjectiveFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&,
            bool useAdjointFileFormat
        );

        //- Construct from Time and dictionary
        //  This constructor is needed so that the objective is added on the
        //  selection table of function objects.
        rearLiftObjectiveFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        rearLiftObjectiveFunctionObject
        (
            const rearLiftObjectiveFunctionObject&
        ) = delete;


    //- Destructor
    virtual ~rearLiftObjectiveFunctionObject();


    // Member Functions

        //- Read and update from dictionary
        virtual bool read(const dictionary&);

        //- Calculate the value of the cost function and write it to file
        virtual bool execute();

        //- Track and write
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const rearLiftObjectiveFunctionObject&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
