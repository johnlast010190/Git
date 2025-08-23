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
    (c) 2011-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "cellValue/cellValue.H"
#include "set/set.H"
#include "cellDeterminant/cellDeterminant.H"
#include "maxFaceWeight/maxFaceWeight.H"
#include "maxSkewness/maxSkewness.H"
#include "maxNonOrthogonality/maxNonOrthogonality.H"
#include "magGradU/magGradU.H"
#include "maxSingularValue/maxSingularValue.H"
#include "maxSingularValueOverU/maxSingularValueOverU.H"

#include "fields/Fields/fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeSensors(Type)                                                      \
    makeSensor(Type);                                                          \
    makeSensorType(cellValue, Type);                                           \
    makeSensorType(set, Type);                                                 \
    makeSensorType(cellDeterminant, Type);                                     \
    makeSensorType(maxFaceWeight, Type);                                       \
    makeSensorType(maxSkewness, Type);                                         \
    makeSensorType(maxNonOrthogonality, Type);                                 \
    makeSensorType(magGradU, Type);                                            \
    makeSensorType(maxSingularValue, Type);                                    \
    makeSensorType(maxSingularValueOverU, Type);


namespace Foam
{
    makeSensors(scalar);
    makeSensors(vector);
    makeSensors(sphericalTensor);
    makeSensors(symmTensor);
    makeSensors(tensor);
}


// ************************************************************************* //
