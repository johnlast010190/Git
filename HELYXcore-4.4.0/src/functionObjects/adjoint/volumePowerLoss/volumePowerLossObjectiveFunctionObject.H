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
    volumePowerLossObjectiveFunctionObject

Description
    Calculate pressure loss defined on the volume.

SourceFiles
    volumePowerLossObjectiveFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef volumePowerLossObjectiveFunctionObject_H
#define volumePowerLossObjectiveFunctionObject_H

#include "objectiveFunctionObject/objectiveFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
               Class volumePowerLossObjectiveFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class volumePowerLossObjectiveFunctionObject
:
    public objectiveFunctionObject
{
protected:

    // Protected data

        //- total objective patch inlet area
        scalar AreaObjPatchInlet_;

        //- total objective patch outlet area
        scalar AreaObjPatchOutlet_;

        //- Switch on adiabatic reporting for compressible flow (default false)
        Switch adiabatic_;

        //- Switch on tangential boundary stress energy transport
        // (default false)
        Switch stressWork_;

        Switch printLossesSeparately_;


public:

    //- Runtime type information
    TypeName("volumePowerLoss");


    // Constructors

        //- Construct from time, dictionary, adjoint solver and use old-style
        //  adjoint file-format
        volumePowerLossObjectiveFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&,
            bool useAdjointFileFormat
        );

        //- Construct from time, dictionary and adjoint solver
        volumePowerLossObjectiveFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        volumePowerLossObjectiveFunctionObject
        (
            const volumePowerLossObjectiveFunctionObject&
        ) = delete;


    //- Destructor
    virtual ~volumePowerLossObjectiveFunctionObject();


    // Member Functions

        //- read and update from dictionary
        virtual bool read(const dictionary&);

        //- Calculate cost function value
        virtual bool execute();

        //- Track and write
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const volumePowerLossObjectiveFunctionObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
