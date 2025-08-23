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
    2008-2011, 2013-2014, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id: writeIfFieldOutsideFunctionObjectTemplates.C,v 78b0d113b99b 2013-02-25 16:12:41Z bgschaid $
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
scalar writeIfFieldOutsideFunctionObject::getMin<volScalarField>() const
{
    const volScalarField &fld=obr_.lookupObject<volScalarField>(fieldName_);

    return min(fld).value();
}

template<>
scalar writeIfFieldOutsideFunctionObject::getMax<volScalarField>() const
{
    const volScalarField &fld=obr_.lookupObject<volScalarField>(fieldName_);

    return max(fld).value();
}

template<class T>
scalar writeIfFieldOutsideFunctionObject::getMin() const
{
    const T &fld=obr_.lookupObject<T>(fieldName_);

    return min(mag(fld)).value();
}

template<class T>
scalar writeIfFieldOutsideFunctionObject::getMax() const
{
    const T &fld=obr_.lookupObject<T>(fieldName_);

    return max(mag(fld)).value();
}

template<class T>
bool writeIfFieldOutsideFunctionObject::check() const
{
    if (obr_.foundObject<T>(fieldName_)) {
        scalar mini=getMin<T>();
        scalar maxi=getMax<T>();

        if (mini<minimum_ || maxi>maximum_) {
            return true;
        }
    }

    return false;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
