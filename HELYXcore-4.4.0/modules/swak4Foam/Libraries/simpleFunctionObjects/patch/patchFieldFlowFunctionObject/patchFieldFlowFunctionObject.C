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

#include "patchFieldFlowFunctionObject.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(patchFieldFlowFunctionObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        patchFieldFlowFunctionObject,
        dictionary
    );



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

patchFieldFlowFunctionObject::patchFieldFlowFunctionObject
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    patchFieldFunctionObject(name,t,dict),
    solver_(obr_,dict)
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

word patchFieldFlowFunctionObject::dirName()
{
    return typeName;
}

scalarField patchFieldFlowFunctionObject::process(const word& fieldName,scalar preset)
{
    return fieldFlow(fieldName,preset);
}

Field<vector> patchFieldFlowFunctionObject::process(const word& fieldName,vector preset)
{
    return fieldFlow(fieldName,preset);
}

Field<sphericalTensor> patchFieldFlowFunctionObject::process(const word& fieldName,sphericalTensor preset)
{
    return fieldFlow(fieldName,preset);
}

Field<symmTensor> patchFieldFlowFunctionObject::process(const word& fieldName,symmTensor preset)
{
    return fieldFlow(fieldName,preset);
}

Field<tensor> patchFieldFlowFunctionObject::process(const word& fieldName,tensor preset)
{
    return fieldFlow(fieldName,preset);
}


} // namespace Foam

// ************************************************************************* //
