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
    2010-2011, 2013-2014, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "expressionSource.H"
#include "meshes/polyMesh/polyMesh.H"
#include "sets/topoSets/cellSet.H"

#include "FieldValueExpressionDriver.H"

namespace Foam {

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
template<class T>
expressionSource<T>::expressionSource
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    FieldValueExpressionDriver(dict,mesh),
    expression_(
        dict.lookup("expression"),
        dict
    )
{
    createWriterAndRead(dict.name().name()+"_"+this->type()+"<"+pTraits<T>::typeName+">");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class T>
expressionSource<T>::~expressionSource()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
tmp<typename expressionSource<T>::resultField> expressionSource<T>::operator()()
{
    clearVariables();
    parse(expression_);
    if (!resultIsTyp<resultField>()) {
        FatalErrorIn("expressionSource<"+word(pTraits<T>::typeName)+">::operator()()")
            << "Result of " << expression_ << " is not a " << pTraits<T>::typeName
                << endl
                << exit(FatalError);
    }

    tmp<resultField> result(new resultField(getResult<resultField>()));

    return result;
}

template
class expressionSource<scalar>;

template
class expressionSource<vector>;

template
class expressionSource<tensor>;

template
class expressionSource<symmTensor>;

template
class expressionSource<sphericalTensor>;

} // end namespace

// ************************************************************************* //
