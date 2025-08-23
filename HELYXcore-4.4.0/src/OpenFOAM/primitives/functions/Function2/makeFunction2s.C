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
    (c) 2020-2021 OpenFOAM Foundation
    (c) 2022 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "primitives/functions/Function2/None/None2.H"
#include "primitives/functions/Function2/Constant/Constant2.H"
#include "primitives/functions/Function2/Zero/ZeroConstant2.H"
#include "primitives/functions/Function2/One/OneConstant2.H"
#include "primitives/functions/Function2/Scale/Scale2.H"
#include "primitives/functions/Function2/UniformTable2/UniformTable2.H"
#include "primitives/functions/Function2/Table2/Table2.H"
#include "primitives/functions/Function2/APIdiffCoef/APIdiffCoef.H"

#include "fields/Fields/fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeFunction2s(Type)                                                   \
    makeFunction2(Type);                                                       \
    makeFunction2Type(None, Type);                                             \
    makeFunction2Type(Constant, Type);                                         \
    makeFunction2Type(ZeroConstant, Type);                                     \
    makeFunction2Type(OneConstant, Type);                                      \
    makeFunction2Type(Scale, Type);                                            \
    makeFunction2Type(UniformTable, Type);                                     \
    makeFunction2Type(Table, Type);


namespace Foam
{
    makeFunction2(label);
    makeFunction2Type(None, label);
    makeFunction2Type(Constant, label);

    makeFunction2s(scalar);
    makeFunction2s(vector);
    makeFunction2s(sphericalTensor);
    makeFunction2s(symmTensor);
    makeFunction2s(tensor);
}


// ************************************************************************* //
