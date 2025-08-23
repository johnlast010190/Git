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
    (c) 2024 Engys Ltd.

Contributors/Copyright:
    2012-2013, 2015-2016 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "dynamicFunctionObjectListProxy.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

#include "meshes/polyMesh/polyMesh.H"
#include "db/IOstreams/IOstreams/IOmanip.H"
#include "include/swakTime.H"

#include "include/swak.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dynamicFunctionObjectListProxy, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        dynamicFunctionObjectListProxy,
        dictionary
    );

    typedef dynamicFunctionObjectListProxy::dynamicDictionaryProvider
        dynProvider;

    defineTypeNameAndDebug(dynProvider, 0);
    defineRunTimeSelectionTable
    (
        dynamicFunctionObjectListProxy::dynamicDictionaryProvider,
        dictionary
    );


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

dynamicFunctionObjectListProxy::dynamicFunctionObjectListProxy
(
    const word& name,
    const Time& t,
    const dictionary& dict,
    const char* providerNameStr
)
:
    functionObjectListProxy(name, t, dict, false)
{
    word providerName(providerNameStr);
    if (providerName.size() == 0)
    {
        providerName=word(dict.lookup("dictionaryProvider"));
    }
    provider_ =
        dynamicDictionaryProvider::New(providerName, dict, (*this));

    if (dict.lookup<bool>("readDuringConstruction"))
    {
        if (writeDebug())
        {
            Info<< this->name()
                << " list initialized during construction" << endl;
        }
        read(dict);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dynamicFunctionObjectListProxy::initFunctions()
{
    string text(provider_->getDictionaryText());
    if (Pstream::parRun())
    {
        string localText = text;
        Pstream::scatter(text);
        if (text != localText)
        {
            Pout<< "WARNING: In dynamicFunctionObjectListProxy::initFunctions() "
                << "for " << name()
                << " the text of the dictionary is different from the master"
                << endl
                << " Overwritten local version with master";
        }
    }
    {
        fileName fName=obr_.time().path()/word(this->name()+".dictionaryText");
        OFstream o(fName);
        o << text.c_str();
    }

    IStringStream inStream(text);

    dynamicDict_.set(new dictionary(inStream));

    {
        fileName fName =
            obr_.time().path()/word(this->name()+".dictionaryDump");
        OFstream o(fName);
        o << dynamicDict_();
    }

    if (!dynamicDict_->found("functions"))
    {
        FatalErrorIn("dynamicFunctionObjectListProxy::initFunctions()")
            << "Dictionary for" << this->name()
            << " does not have an entry 'functions'" << nl << exit(FatalError);

    }

    functions_.set(new functionObjectList(time(), dynamicDict_()));

    if (writeDebug())
    {
        Info<< this->name() << " list initialized with "
            << functions_->size() << " FOs" << endl;
    }
}


autoPtr<dynamicFunctionObjectListProxy::dynamicDictionaryProvider>
dynamicFunctionObjectListProxy::dynamicDictionaryProvider::New(
    const word& type,
    const dictionary& dict,
    const dynamicFunctionObjectListProxy &owner
)
{
    const auto ctor =
        ctorTableLookup
        (
            "dynamicFunctionObjectListProxy::dynamicDictionaryProvider type",
            dictionaryConstructorTable_(),
            type
        );
    if (debug)
    {
        Pout<< "Creating dictionary provider of type " << type << endl;
    }

    return
        autoPtr<dynamicFunctionObjectListProxy::dynamicDictionaryProvider>
        (
            ctor(dict,owner)
        );
}

} // namespace Foam

// ************************************************************************* //
