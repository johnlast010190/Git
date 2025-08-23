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
    (c) 2024 Engys Ltd.

Contributors/Copyright:
    2012-2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

namespace Foam {


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class T>
void StackExpressionResult::pushInternal(ExpressionResult &atEnd)
{
    Dbug<< "pushInternal(ExpressionResult &atEnd) "
            << pTraits<T>::typeName << endl;

    Field<T> oldValue(getResult<T>(false));
    autoPtr<Field<T>> pNewValue(new Field<T>(oldValue.size()+1));
    Field<T> &newValue=pNewValue();

    forAll(oldValue,i) {
        newValue[i]=oldValue[i];
    }
    newValue[oldValue.size()]=atEnd.getResult<T>(true)()[0];

    Dbug<< "New value:" << newValue << endl;

    this->setResult(
        pNewValue.ptr(),
        this->isPoint()
    );
}

template<class T>
ExpressionResult StackExpressionResult::popInternal()
{
    Field<T> oldValue(getResult<T>(false));
    autoPtr<Field<T>> pNewValue(new Field<T>(oldValue.size()-1));
    Field<T> &newValue=pNewValue();

    forAll(newValue,i) {
        newValue[i]=oldValue[i];
    }
    Field<T> result(1,oldValue[newValue.size()]);

    this->setResult(
        pNewValue.ptr(),
        this->isPoint()
    );

    return ExpressionResult(result);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}

// ************************************************************************* //
