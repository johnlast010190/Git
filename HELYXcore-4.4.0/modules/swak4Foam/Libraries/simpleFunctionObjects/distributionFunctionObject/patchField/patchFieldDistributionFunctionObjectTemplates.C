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

#include "distributionFunctionObject/field/fieldDistributionFunctionObject.H"
#include "fields/volFields/volFields.H"
#include "db/IOstreams/IOstreams/IOmanip.H"
#include "fvMesh/fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <typename T>
void patchFieldDistributionFunctionObject::getDistributionInternal(
    autoPtr<SimpleDistribution<T>> &dist
) {
    const fvMesh &mesh=refCast<const fvMesh>(obr_);

    typedef VolField<T> volTField;
    typedef SurfaceField<T> surfaceTField;
    typedef PointField<T> pointTField;

    bool firstTime=true;

    forAll(patchIDs_,i) {
        label patchID=patchIDs_[i];
        autoPtr<SimpleDistribution<T>> partial;

        if (mesh.foundObject<volTField>(fieldName())) {
            partial=setDataScalar(
                mesh.lookupObject<volTField>(
                    fieldName()
                ).boundaryField()[patchID],
                mesh.boundary()[patchID].magSf()
            );
        }
        if (mesh.foundObject<surfaceTField>(fieldName())) {
            partial=setDataScalar(
                mesh.lookupObject<surfaceTField>(
                    fieldName()
                ).boundaryField()[patchID],
                mesh.boundary()[patchID].magSf()
            );
        }
        if (mesh.foundObject<pointTField>(fieldName())) {
            partial=setDataScalar(
                mesh.lookupObject<pointTField>(
                    fieldName()
                ).boundaryField()[patchID].patchInternalField()(),
                scalarField(
                    mesh.lookupObject<pointTField>(
                        fieldName()
                    ).boundaryField()[patchID].size(),
                    1
                )
            );
        }

        if (partial.valid()) {
            if (firstTime) {
                firstTime=false;
                dist=partial;
            } else {
                SimpleDistribution<T> &d=dist();
                d=d+partial();
            }
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
