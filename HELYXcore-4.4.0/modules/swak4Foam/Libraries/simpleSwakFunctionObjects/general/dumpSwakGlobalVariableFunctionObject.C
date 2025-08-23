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
    2011-2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id:  $
\*---------------------------------------------------------------------------*/

#include "dumpSwakGlobalVariableFunctionObject.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dumpSwakGlobalVariableFunctionObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        dumpSwakGlobalVariableFunctionObject,
        dictionary
    );



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

dumpSwakGlobalVariableFunctionObject::dumpSwakGlobalVariableFunctionObject
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    timelineFunctionObject(name,t,dict),
    globalScope_(dict.lookup("globalScope")),
    globalName_(dict.lookup("globalName"))
{
    const string warnSwitch="IKnowThatThisFunctionObjectMayWriteExcessiveAmountsOfData";
    if (!dict.lookupOrDefault<bool>(warnSwitch,false)) {
        WarningIn("dumpSwakGlobalVariableFunctionObject::dumpSwakGlobalVariableFunctionObject")
            << "This functionObject may write huge amounts of data. "
                << "If you understand the risks set the switch " << warnSwitch
                << " to 'true' to get rid of this warning"
                << endl;
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

word dumpSwakGlobalVariableFunctionObject::dirName()
{
    return typeName;
}

wordList dumpSwakGlobalVariableFunctionObject::fileNames()
{
    return wordList(1,name());
}

stringList dumpSwakGlobalVariableFunctionObject::columnNames()
{
    return stringList(1,"No way to know how much data will follow");
}

void dumpSwakGlobalVariableFunctionObject::writeSimple()
{

    if (verbose()) {
        Info<< "Global " << name() << " : ";
    }

    ExpressionResult value(
        GlobalVariablesRepository::getGlobalVariables(
                obr_
        ).get(
            globalName_,
            wordList(1,globalScope_)
        )
    );

    word rType(value.valueType());

    if (rType==pTraits<scalar>::typeName) {
        writeTheData<scalar>(value);
    } else if (rType==pTraits<vector>::typeName) {
        writeTheData<vector>(value);
    } else if (rType==pTraits<tensor>::typeName) {
        writeTheData<tensor>(value);
    } else if (rType==pTraits<symmTensor>::typeName) {
        writeTheData<symmTensor>(value);
    } else if (rType==pTraits<sphericalTensor>::typeName) {
        writeTheData<sphericalTensor>(value);
    } else {
        WarningIn("dumpSwakGlobalVariableFunctionObject::writeSimple()")
            << "Don't know how to handle type " << rType
                << endl;
    }

    if (verbose()) {
        Info<< endl;
    }
}

} // namespace Foam

// ************************************************************************* //
