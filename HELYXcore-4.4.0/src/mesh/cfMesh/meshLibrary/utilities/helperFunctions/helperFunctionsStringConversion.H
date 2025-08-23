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
    (c) Creative Fields, Ltd.

Authors
    Franjo Juretic (franjo.juretic@c-fields.com)

InNamespace
    Foam::help

Description
    String conversion functions

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef helperFunctionsStringConversion_H
#define helperFunctionsStringConversion_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "primitives/strings/word/word.H"
#include "primitives/ints/label/label.H"
#include "primitives/Scalar/scalar/scalar.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Namespace help functions Declaration
\*---------------------------------------------------------------------------*/

namespace help
{

    //- convert the text to scalar
    scalar textToScalar(const word& w);

    //- convert the text to label
    label textToLabel(const word& w);

    //- convert the scalar value into text
    word scalarToText(const scalar s);

    //- convert the integer value into text
    word labelToText(const label l);

} // End namespace help

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
