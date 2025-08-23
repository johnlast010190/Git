/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : dev
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
    (c) ICE Stroemungsfoschungs GmbH

Contributors/Copyright:
    2011-2014, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>
    2013 Bruno Santos <wyldckat@gmail.com>

 SWAK Revision: $Id:  $
\*---------------------------------------------------------------------------*/

#include "removeGlobalVariable.H"

#include "repositories/GlobalVariablesRepository.H"

namespace Foam {
    defineTypeNameAndDebug(removeGlobalVariable,0);
}

Foam::removeGlobalVariable::removeGlobalVariable
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
    :
    obr_(obr)
{
    read(dict);
    execute();
}

Foam::removeGlobalVariable::~removeGlobalVariable()
{}

void Foam::removeGlobalVariable::timeSet()
{
    // Do nothing
}

void Foam::removeGlobalVariable::read(const dictionary& dict)
{
    names_=wordList(dict.lookup("globalVariables"));
    scope_=word(dict.lookup("globalScope"));
}

void Foam::removeGlobalVariable::execute()
{
   forAll(names_,i) {
        const word &name=names_[i];

        bool removed=GlobalVariablesRepository::getGlobalVariables(
            obr_
        ).removeValue(
            name,
            scope_
        );
        if (!removed) {
            WarningIn("Foam::removeGlobalVariable::read(const dictionary& dict)")
                << "Variable " << name << " in scope " << scope_
                    << " not removed" << endl;
        }
   }
}


void Foam::removeGlobalVariable::end()
{
}

#ifdef FOAM_IOFILTER_WRITE_NEEDS_BOOL
bool
#else
void
#endif
Foam::removeGlobalVariable::write()
{
#ifdef FOAM_IOFILTER_WRITE_NEEDS_BOOL
    return true;
#endif
}

void Foam::removeGlobalVariable::clearData()
{
}

// ************************************************************************* //
