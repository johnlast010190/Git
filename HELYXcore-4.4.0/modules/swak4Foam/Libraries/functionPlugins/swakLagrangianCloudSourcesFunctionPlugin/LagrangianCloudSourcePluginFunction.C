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
    (c) 1991-2008 OpenCFD Ltd.

Contributors/Copyright:
    2012-2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "LagrangianCloudSourcePluginFunction.H"
#include "FieldValueExpressionDriver.H"

#include "fields/cloud/cloud.H"
#include "db/IOstreams/IOstreams/IOmanip.H"

namespace Foam {

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

LagrangianCloudSourcePluginFunction::LagrangianCloudSourcePluginFunction(
    const FieldValueExpressionDriver &parentDriver,
    const word &name,
    const word &returnValueType,
    const string &additionalArgs
):
    FieldValuePluginFunction(
        parentDriver,
        name,
        returnValueType,
        string("cloudName primitive word")+ (
            additionalArgs==""
            ?
            additionalArgs
            :
            string(string(",")+additionalArgs)
        )
    )
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Ostream& LagrangianCloudSourcePluginFunction::listAvailableClouds(Ostream &o)
{
    o << nl << nl << "Available clouds in " << mesh().name() << endl;
    typedef HashTable<const cloud *> cloudTable;
    cloudTable clouds=mesh().lookupClass<cloud>();
    if (clouds.size()==0) {
        o << " No clouds available\n" << endl;
    }
    const label nameWidth=20;
    o << setw(nameWidth) << "Name" << " | " << "Type" << endl;
    o << "-------------------------------------------------------------" << endl;
    forAllConstIter(cloudTable,clouds,it) {
        o << setw(nameWidth) << it.key() << " | "
            << (*(*it)).type() << endl;
    }
    o << "-------------------------------------------------------------" << endl;
    if (mesh().foundObject<cloud>(cloudName())) {
        const cloud& obj = mesh().lookupObject<cloud>(cloudName());
        o << "Cloud " << cloudName() << " has type "
            << obj.type() << " typeid:"
            << typeid(obj).name()
            << endl;
    }
    return o;
}


void LagrangianCloudSourcePluginFunction::setArgument(
    label index,
    const word &name
)
{
    assert(index==0);
    cloudName_=name;
}


// * * * * * * * * * * * * * * * Concrete implementations * * * * * * * * * //


} // namespace

// ************************************************************************* //
