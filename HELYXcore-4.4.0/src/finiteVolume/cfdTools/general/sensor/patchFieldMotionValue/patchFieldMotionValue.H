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
    (c) 2017 Engys Ltd.

Class
    Foam::sensorTypes::patchFieldMotionValue

Description
    Sensor type based on patch field (surface<Type>Field) at a given patch.

    \verbatim
        patchName    <patch name>
        fieldName    <field name>
    \endverbatim

SourceFiles
    patchFieldMotionValue.C

\*---------------------------------------------------------------------------*/

#ifndef patchFieldMotionValue_H
#define patchFieldMotionValue_H

#include "cfdTools/general/sensor/sensor/sensor.H"
#include "sets/topoSetSource/topoSetSource.H"
#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace sensorTypes
{

/*---------------------------------------------------------------------------*\
                           Class patchFieldMotionValue Declaration
\*---------------------------------------------------------------------------*/

class patchFieldMotionValue
:
    public sensor<scalar>
{
    // Private data

        //- patch name
        word patchName_;

        //- reference point
        vector point0_;

        //- direction
        vector direction_;

        //- Linear velocity
        autoPtr<Function1<vector>> velocity_;


public:

    // Runtime type information
    TypeName("patchFieldMotionValue");


    // Constructors

        //- Construct from components
        patchFieldMotionValue
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Copy constructor
        patchFieldMotionValue(const patchFieldMotionValue& cv);


    //- Destructor
    virtual ~patchFieldMotionValue();


    // Member Functions

        //- Return patchFieldMotionValue value
        virtual tmp<Field<scalar>> valueField() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace sensorTypes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
