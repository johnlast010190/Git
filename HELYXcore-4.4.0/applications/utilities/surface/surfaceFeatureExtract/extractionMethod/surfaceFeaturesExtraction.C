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
    (c) 2017 OpenCFD Ltd.
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "surfaceFeaturesExtraction.H"
#include "db/dictionary/dictionary.H"
#include "containers/Lists/ListOps/ListOps.H"
#include "db/error/error.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceFeaturesExtraction
{
    defineTypeName(method);
    defineRunTimeSelectionTable
    (
        method,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceFeaturesExtraction::method::method()
:
    includedAngle_(0),
    geometricTestOnly_(Switch::NO)
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::surfaceFeaturesExtraction::method::~method()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::surfaceFeaturesExtraction::method>
Foam::surfaceFeaturesExtraction::method::New
(
    const dictionary& dict
)
{
    const word methodName = dict.lookup("extractionMethod");

    const auto ctor =
        ctorTableLookup
        (
            "extractionMethod",
            dictionaryConstructorTable_(),
            methodName
        );
    return autoPtr<method>(ctor(dict));
}


// ************************************************************************* //
