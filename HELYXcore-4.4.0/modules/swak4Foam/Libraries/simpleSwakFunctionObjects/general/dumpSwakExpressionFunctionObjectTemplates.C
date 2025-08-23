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
    2011, 2013, 2016 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "dumpSwakExpressionFunctionObject.H"
#include "fields/volFields/volFields.H"
#include "db/IOstreams/IOstreams/IOmanip.H"
#include "fvMesh/fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<>
void dumpSwakExpressionFunctionObject::writeValue(Ostream &o,const scalar &val,unsigned int &w)
{
    o << setw(w) << val;
}

template<class Type>
void dumpSwakExpressionFunctionObject::writeValue(Ostream &o,const Type &val,unsigned int &w)
{
    for (label j=0;j<Type::nComponents;j++) {
        o << setw(w) << val[j];
    }
}

template <class T>
void dumpSwakExpressionFunctionObject::writeTheData(CommonValueExpressionDriver &driver)
{
    List<Field<T>> results(Pstream::nProcs());
    results[Pstream::myProcNo()]=driver.getResult<T>();

    Pstream::gatherList(results);

    if (Pstream::master()) {
        writeTime(name(),time().value());
        forAll(results,procNo) {
            writeData(name(),results[procNo]);
        }
        endData(name());
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
