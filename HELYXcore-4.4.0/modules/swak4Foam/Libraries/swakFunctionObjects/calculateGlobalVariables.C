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

#include "calculateGlobalVariables.H"

#include "repositories/GlobalVariablesRepository.H"

namespace Foam {
    defineTypeNameAndDebug(calculateGlobalVariables,0);
}

Foam::calculateGlobalVariables::calculateGlobalVariables
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
    :
    obr_(obr),
    driver_(
        CommonValueExpressionDriver::New(
            dict,
            refCast<const fvMesh>(obr)
        )
    ),
    toGlobalNamespace_(dict.lookup("toGlobalNamespace")),
    toGlobalVariables_(dict.lookup("toGlobalVariables")),
    noReset_(dict.lookupOrDefault<bool>("noReset",false))
{
    if (debug) {
        Info<< "calculateGlobalVariables " << name << " created" << endl;
    }

    if (!dict.found("noReset")) {
        WarningIn("calculateGlobalVariables::calculateGlobalVariables")
            << "No entry 'noReset' in " << dict.name()
                << ". Assumig 'false'"<< endl;

    }

    driver_->createWriterAndRead(name+"_"+type());

    executeAndWriteToGlobal();
}

Foam::calculateGlobalVariables::~calculateGlobalVariables()
{}

void Foam::calculateGlobalVariables::executeAndWriteToGlobal()
{
    // this also sets the variables
    driver_->clearVariables();

    forAll(toGlobalVariables_,i) {
        const word &name=toGlobalVariables_[i];
        if (debug) {
            Info<< "Getting variable " << name << endl;
        }

        ExpressionResult &res=GlobalVariablesRepository::getGlobalVariables(
            obr_
        ).addValue(
            name,
            toGlobalNamespace_,
            const_cast<const CommonValueExpressionDriver&>(
                driver_()
            ).variable(name)
        );

        if (noReset_) {
            res.noReset();
        }

        if (debug) {
            Pout<< "Has value "
                << const_cast<const CommonValueExpressionDriver&>(
                    driver_()
                ).variable(name) << endl;
        }
    }
}

void Foam::calculateGlobalVariables::timeSet()
{
    // Do nothing
}

void Foam::calculateGlobalVariables::read(const dictionary& dict)
{
    WarningIn("Foam::calculateGlobalVariables::read(const dictionary& dict)")
        << "This function object does not honor changes during runtime"
            << endl;
}

#ifdef FOAM_IOFILTER_WRITE_NEEDS_BOOL
bool
#else
void
#endif
Foam::calculateGlobalVariables::write()
{
    executeAndWriteToGlobal();

    // make sure that the stored Variables are consistently written
    driver_->tryWrite();

#ifdef FOAM_IOFILTER_WRITE_NEEDS_BOOL
    return true;
#endif
}


void Foam::calculateGlobalVariables::end()
{
}

void Foam::calculateGlobalVariables::execute()
{
}

void Foam::calculateGlobalVariables::clearData()
{
}

// ************************************************************************* //
