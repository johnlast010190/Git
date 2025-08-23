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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2018-2023 ENGYS

\*---------------------------------------------------------------------------*/

#include "primitives/functions/Function1/Constant/Constant.H"
#include "primitives/functions/Function1/Zero/ZeroConstant.H"
#include "primitives/functions/Function1/One/OneConstant.H"
#include "primitives/functions/Function1/PolynomialEntry/PolynomialEntry.H"
#include "primitives/functions/Function1/NSRDS0/NSRDS0.H"
#include "primitives/functions/Function1/NSRDS1/NSRDS1.H"
#include "primitives/functions/Function1/NSRDS2/NSRDS2.H"
#include "primitives/functions/Function1/NSRDS3/NSRDS3.H"
#include "primitives/functions/Function1/NSRDS4/NSRDS4.H"
#include "primitives/functions/Function1/NSRDS5/NSRDS5.H"
#include "primitives/functions/Function1/NSRDS6/NSRDS6.H"
#include "primitives/functions/Function1/NSRDS7/NSRDS7.H"
#include "primitives/functions/Function1/NSRDS14/NSRDS14.H"
#include "primitives/functions/Function1/Sine/Sine.H"
#include "primitives/functions/Function1/Square/Square.H"
#include "primitives/functions/Function1/CSV/CSV.H"
#include "primitives/functions/Function1/Table/Table.H"
#include "primitives/functions/Function1/TableFile/TableFile.H"
#include "primitives/functions/Function1/Scale/Scale.H"
#include "primitives/functions/Function1/SelectEntry/SelectEntry.H"
#include "primitives/functions/Function1/xyzPolynomial/xyzPolynomial.H"
#include "primitives/functions/Function1/NonUniformTable/NonUniformTable1.H"
#include "primitives/functions/Function1/uniformTable1/uniformTable1.H"

#include "fields/Fields/fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeFunction1s(Type)                                                   \
    makeFunction1(Type);                                                       \
    makeFunction1Type(Constant, Type);                                         \
    makeFunction1Type(ZeroConstant, Type);                                     \
    makeFunction1Type(OneConstant, Type);                                      \
    makeFunction1Type(Polynomial, Type);                                       \
    makeFunction1Type(xyzPolynomial, Type);                                    \
    makeFunction1Type(Sine, Type);                                             \
    makeFunction1Type(NonUniformTable, Type);                                  \
    makeFunction1Type(Square, Type);                                           \
    makeFunction1Type(CSV, Type);                                              \
    makeFunction1Type(Table, Type);                                            \
    makeFunction1Type(TableFile, Type);                                        \
    makeFunction1Type(Scale, Type);                                            \
    makeFunction1Type(Select, Type);                                           \
    makeFunction1TypeAndLookup(NSRDS0, Type, NSRDSfunc0);                      \
    makeFunction1TypeAndLookup(NSRDS1, Type, NSRDSfunc1);                      \
    makeFunction1TypeAndLookup(NSRDS2, Type, NSRDSfunc2);                      \
    makeFunction1TypeAndLookup(NSRDS3, Type, NSRDSfunc3);                      \
    makeFunction1TypeAndLookup(NSRDS4, Type, NSRDSfunc4);                      \
    makeFunction1TypeAndLookup(NSRDS5, Type, NSRDSfunc5);                      \
    makeFunction1TypeAndLookup(NSRDS6, Type, NSRDSfunc6);                      \
    makeFunction1TypeAndLookup(NSRDS7, Type, NSRDSfunc7);                      \
    makeFunction1TypeAndLookup(NSRDS14, Type, NSRDSfunc14);                    \
    makeFunction1Type(UniformTable, Type);

#define makeFieldFunction1s(Type)                                              \
    makeFunction1(Type);                                                       \
    makeFunction1Type(Constant, Type);                                         \
    makeFunction1Type(Table, Type);                                            \
    makeFunction1Type(TableFile, Type);

namespace Foam
{
    makeFunction1(label);
    makeFunction1Type(Constant, label);
    // Polynomial functions and interpolation do evaluate to label
    // Instead evaluate a scalar and convert to label as appropriate

    makeFunction1s(scalar);
    makeFunction1s(vector);
    makeFunction1s(sphericalTensor);
    makeFunction1s(symmTensor);
    makeFunction1s(tensor);

    makeFunction1(bool);
    makeFunction1Type(Table, bool);
    makeFunction1Type(Constant, bool);

    makeFieldFunction1s(scalarField);
}


// ************************************************************************* //
