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
    2008-2013, 2016-2017 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "fieldDistributionFunctionObject.H"
#include "fields/volFields/volFields.H"
#include "fields/GeometricFields/pointFields/pointFieldsFwd.H"
#include "db/IOstreams/IOstreams/IOmanip.H"
#include "fvMesh/fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <typename T>
void fieldDistributionFunctionObject::getDistributionInternal(
    autoPtr<SimpleDistribution<T>> &dist
) {
    const fvMesh &mesh=refCast<const fvMesh>(obr_);

    typedef VolField<T> volTField;
    typedef SurfaceField<T> surfaceTField;
    typedef PointField<T> pointTField;

    if (mesh.foundObject<volTField>(fieldName_)) {
        dist=setDataScalar(
            mesh.lookupObject<volTField>(
                fieldName_
            ).internalField(),
            mesh.V()
        );
        return;
    }
    if (mesh.foundObject<surfaceTField>(fieldName_)) {
        dist=setDataScalar(
            mesh.lookupObject<surfaceTField>(
                fieldName_
            ).internalField(),
            mesh.magSf()
        );
        return;
    }
    if (mesh.foundObject<pointTField>(fieldName_)) {
        dist=setDataScalar(
            mesh.lookupObject<pointTField>(
                fieldName_
            ).internalField(),
            scalarField(mesh.nPoints(),1)
        );
        return;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
