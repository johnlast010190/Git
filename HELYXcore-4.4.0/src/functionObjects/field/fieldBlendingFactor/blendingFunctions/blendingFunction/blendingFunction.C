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
    (c) 2011-2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "db/dictionary/dictionary.H"
#include "fieldBlendingFactor/blendingFunctions/blendingFunction/blendingFunction.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{

defineTypeNameAndDebug(blendingFunction, 0);
defineRunTimeSelectionTable(blendingFunction, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

blendingFunction::blendingFunction(const dictionary& dict, bool dum)
{}

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<blendingFunction> blendingFunction::New
(
    const dictionary& dict
)
{
    word functionType = dict.lookup("functionType");

    const auto ctor =
        ctorTableLookup
        (
            "function type",
            dictionaryConstructorTable_(),
            functionType
        );

    // additional dummy variable needed for macros to work
    bool dum = true;
    return autoPtr<blendingFunction>(ctor(dict, dum));
}


} // End namespace Foam

// ************************************************************************* //
