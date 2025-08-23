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
    2008-2011, 2013, 2015-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>
    2014 David Huckaby <e.david.huckaby@netl.doe.gov>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "volume/volumeFieldFunctionObject/volumeFieldFunctionObject.H"
#include "fields/volFields/volFields.H"
#include "db/IOstreams/IOstreams/IOmanip.H"
#include "fvMesh/fvMesh.H"

#include "clouds/derived/basicKinematicCloud/basicKinematicCloud.H"
#include "misc/objectRegistryUtility/objectRegistryUtility.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class T>
bool writeFieldsGeneralFunctionObject::writeField(const word &name) const
{
    if (obr_.foundObject<T>(name)) {
        obr_.lookupObject<T>(name).write();
        return true;
    } else {
        return false;
    }
}

template<class Type>
bool writeFieldsGeneralFunctionObject::writeCloud(const word &name) const
{
    if (obr_.foundObject<Type>(name)) {
        Info<< "\twriting Cloud: " << name << endl;
        lookupObject<Type>(obr_,name).write();
        return true;
    } else {
        return false;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
