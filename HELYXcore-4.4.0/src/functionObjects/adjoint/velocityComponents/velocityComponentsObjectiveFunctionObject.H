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
    (c) 2020-2024 Engys Ltd.

Class
    Foam::functionObjects::objectiveFunctionObject::velocityComponentsFunctionObject

Description


SourceFiles
    velocityComponentsObjectiveFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef velocityComponentsObjectiveFunctionObject_H
#define velocityComponentsObjectiveFunctionObject_H

#include "objectiveFunctionObject/objectiveFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

class velocityComponentsObjectiveFunctionObject
:
    public objectiveFunctionObject
{

protected:

    vector axis_;

    vector origin_;

    vector swirlDirection_;

    bool volumeBased_;

    scalar objectiveVolume_;
    scalar objectiveArea_;

    //-Actual values
    class velocityData{
    public:
        scalar target_;
        scalar weights_;
        scalar deviation_;
        bool active_;

        velocityData()
        {
            target_ = 0.;
            weights_ = 0.;
            deviation_ = -1;
            active_ = false;
        }
    };

    velocityData axialVel_, peripheralVel_, radialVel_, peripheralOverAxialVel_;

    //- Write the header of the output file in postProcessing
    virtual void writeFileHeader(Ostream&);

    bool isAnActivePatch(const label& pI);

    void setVariables(const dictionary& objectiveDict);
    void calculateObjectiveVolume();
    void calculateObjectiveArea();

    scalar getVolumeAxialDeviation();
    scalar getVolumePeripheralDeviation();
    scalar getVolumeRadialDeviation();
    scalar getVolumePeripheralOverAxialDeviation();

    scalar getSurfaceAxialDeviation();
    scalar getSurfacePeripheralDeviation();
    scalar getSurfaceRadialDeviation();
    scalar getSurfacePeripheralOverAxialDeviation();


public:

    TypeName("velocityComponents");

    // Constructors

        velocityComponentsObjectiveFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&,
            bool useAdjointFileFormat
        );

        velocityComponentsObjectiveFunctionObject
        (
            const word&,
            const Time&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        velocityComponentsObjectiveFunctionObject
        (
            const velocityComponentsObjectiveFunctionObject&
        ) = delete;


    //- Destructor
    virtual ~velocityComponentsObjectiveFunctionObject();


    // Member Functions

        //- Read and update from dictionary
        virtual bool read(const dictionary&);

        //- Calculate cost function value
        virtual bool execute();

        //- Track and write
        virtual bool write();


    // Member Operators

    //- Disallow default bitwise assignment
    void operator=(const velocityComponentsObjectiveFunctionObject&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
