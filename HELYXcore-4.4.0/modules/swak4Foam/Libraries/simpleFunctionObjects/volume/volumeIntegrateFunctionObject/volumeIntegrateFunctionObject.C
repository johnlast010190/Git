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
    2008-2011, 2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "volumeIntegrateFunctionObject.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(volumeIntegrateFunctionObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        volumeIntegrateFunctionObject,
        dictionary
    );



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

volumeIntegrateFunctionObject::volumeIntegrateFunctionObject
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    volumeFieldFunctionObject(name,t,dict)
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

word volumeIntegrateFunctionObject::dirName()
{
    return typeName;
}

scalarField volumeIntegrateFunctionObject::process(const word& fieldName,scalar preset)
{
    return integrate(fieldName,preset);
}

Field<vector> volumeIntegrateFunctionObject::process(const word& fieldName,vector preset)
{
    return integrate(fieldName,preset);
}

Field<sphericalTensor> volumeIntegrateFunctionObject::process(const word& fieldName,sphericalTensor preset)
{
    return integrate(fieldName,preset);
}

Field<symmTensor> volumeIntegrateFunctionObject::process(const word& fieldName,symmTensor preset)
{
    return integrate(fieldName,preset);
}

Field<tensor> volumeIntegrateFunctionObject::process(const word& fieldName,tensor preset)
{
    return integrate(fieldName,preset);
}

stringList volumeIntegrateFunctionObject::columnNames()
{
    return stringList(1,"integrale");
}

} // namespace Foam

// ************************************************************************* //
