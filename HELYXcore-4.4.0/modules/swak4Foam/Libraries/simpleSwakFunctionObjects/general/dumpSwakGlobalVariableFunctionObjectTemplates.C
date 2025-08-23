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
    2011, 2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id:  $
\*---------------------------------------------------------------------------*/

#include "dumpSwakGlobalVariableFunctionObject.H"
#include "fields/volFields/volFields.H"
#include "db/IOstreams/IOstreams/IOmanip.H"
#include "fvMesh/fvMesh.H"
#include "repositories/GlobalVariablesRepository.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<>
void dumpSwakGlobalVariableFunctionObject::writeValue(
    Ostream &o,
    const scalar &val,
    unsigned int &w
)
{
    o << setw(w) << val;
}

template<class Type>
void dumpSwakGlobalVariableFunctionObject::writeValue(
    Ostream &o,
    const Type &val,
    unsigned int &w
)
{
    for (label j=0;j<Type::nComponents;j++) {
        o << setw(w) << val[j];
    }
}

template <class T>
void dumpSwakGlobalVariableFunctionObject::writeTheData(
    ExpressionResult &value
)
{
    Field<T> result(value.getResult<T>());

    if (Pstream::master()) {
        writeTime(name(),time().value());
        writeData(name(),result);
        endData(name());
    } else {
        Pout<< "My data is lost because for dumpSwakGlobalVariableFunctionObject"
            << " only the masters data gets written" << endl;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
