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

#include "volumeFieldFunctionObject.H"
#include "fields/volFields/volFields.H"
#include "db/IOstreams/IOstreams/IOmanip.H"
#include "fvMesh/fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class T>
void volumeFieldFunctionObject::findFields(wordList& typeFieldNames, boolList& foundFields)
{
    typeFieldNames.setSize(fieldNames_.size());
    label typeFieldI = 0;

    forAll(fieldNames_, fieldI)
    {
        const word& fldName = fieldNames_[fieldI];

        if (obr_.foundObject<T>(fldName))
        {
            typeFieldNames[typeFieldI++] = fldName;
            foundFields[fieldI] = true;
        }
    }

    typeFieldNames.setSize(typeFieldI);
}


template <class T>
void volumeFieldFunctionObject::processAndWrite(const word& fieldName)
{
    const VolField<T>& fld =
        obr_.lookupObject<VolField<T>>
        (
            fieldName
        );

    // Make sure all processors call sample
    Field<T> vals(process(fieldName,-VGREAT*pTraits<T>::one));

    if (Pstream::master())
    {
        writeTime(fieldName,fld.time().value());
        writeData(fieldName,vals);
        endData(fieldName);
    }
}

template <class T>
void volumeFieldFunctionObject::processAndWrite(const wordList& typeFields)
{
    forAll(typeFields, i)
    {
        processAndWrite<T>(typeFields[i]);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
