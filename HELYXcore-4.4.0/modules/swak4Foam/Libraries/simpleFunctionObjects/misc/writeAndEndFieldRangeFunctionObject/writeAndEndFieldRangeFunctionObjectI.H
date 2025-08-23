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
    2012-2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "volume/volumeFieldFunctionObject/volumeFieldFunctionObject.H"
#include "fields/volFields/volFields.H"
#include "db/IOstreams/IOstreams/IOmanip.H"
#include "fvMesh/fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<>
scalar writeAndEndFieldRangeFunctionObject::getMin<volScalarField>() const
{
    const volScalarField &fld=obr_.lookupObject<volScalarField>(fieldName_);

    return min(fld).value();
}

template<>
scalar writeAndEndFieldRangeFunctionObject::getMax<volScalarField>() const
{
    const volScalarField &fld=obr_.lookupObject<volScalarField>(fieldName_);

    return max(fld).value();
}

template<class T>
scalar writeAndEndFieldRangeFunctionObject::getMin() const
{
    const T &fld=obr_.lookupObject<T>(fieldName_);

    return min(mag(fld)).value();
}

template<class T>
scalar writeAndEndFieldRangeFunctionObject::getMax() const
{
    const T &fld=obr_.lookupObject<T>(fieldName_);

    return max(mag(fld)).value();
}

template<class T>
bool writeAndEndFieldRangeFunctionObject::check() const
{
    if (obr_.foundObject<T>(fieldName_)) {
        scalar mini=getMin<T>();
        scalar maxi=getMax<T>();

        if (mini<minimum_ || maxi>maximum_) {
            // make sure this field gets written, even if its not set to AUTO_WRITE
            obr_.lookupObject<T>(fieldName_).write();

            bool result=const_cast<Time&>(time()).writeNow();

            // makes sure that all processes finished writing before dumping
            reduce(result,andOp<bool>());

            Info<< " Finishing run, because the field " << fieldName_
                << " exceeds the valid limits [ " << minimum_ << " , "
                << maximum_ << " ] with the current limits [ "
                << mini << " , " << maxi << " ]\n" << endl;

            return true;
        }
    }

    return false;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
