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
    Foam::sensorTypes::maxNonOrthogonality

Description
    Sensor based on the maximum non-orthogonality of the faces of a cell

SourceFiles
    maxNonOrthogonality.C

\*---------------------------------------------------------------------------*/

#ifndef maxNonOrthogonality_H
#define maxNonOrthogonality_H

#include "cfdTools/general/sensor/sensor/sensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace sensorTypes
{

/*---------------------------------------------------------------------------*\
                           Class maxNonOrthogonality Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class maxNonOrthogonality
:
    public sensor<Type>
{

public:

    // Runtime type information
    TypeName("maxNonOrthogonality");


    // Constructors

        //- Construct from components
        maxNonOrthogonality
        (
            const fvMesh& mesh,
            const dictionary& sensorDict
        );


    //- Destructor
    virtual ~maxNonOrthogonality();


    // Member Functions

        tmp<Field<scalar>> valueField() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace sensorTypes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "cfdTools/general/sensor/maxNonOrthogonality/maxNonOrthogonality.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
