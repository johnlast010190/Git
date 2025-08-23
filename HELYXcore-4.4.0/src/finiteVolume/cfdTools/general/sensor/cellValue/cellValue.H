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
    Foam::sensorTypes::cellValue

Description
    Sensor type based on cell value of vol<Type>Field in a given cell.

    \verbatim
        cellID    <cell index>
    \endverbatim

SourceFiles
    cellValue.C

\*---------------------------------------------------------------------------*/

#ifndef cellValue_H
#define cellValue_H

#include "cfdTools/general/sensor/sensor/sensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace sensorTypes
{

/*---------------------------------------------------------------------------*\
                           Class cellValue Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class cellValue
:
    public sensor<Type>
{
    // Private data

        //- cell ID
        label cellID_;


public:

    // Runtime type information
    TypeName("cellValue");


    // Constructors

        //- Construct from components
        cellValue
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Copy constructor
        cellValue(const cellValue<Type>& cv);


    //- Destructor
    virtual ~cellValue();


    // Member Functions

        //- Return cellValue value
        Type value(const VolField<Type>& volField) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace sensorTypes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "cfdTools/general/sensor/cellValue/cellValue.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
