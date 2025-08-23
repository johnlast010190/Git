/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  HELYX (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.4.0
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
    (c) 2011-2013 OpenFOAM Foundation
    (c) 2022 Engys Ltd.

\*----------------------------------------------------------------------------*/

#include "forceActuationDiskSource.H"
#include "fields/volFields/volFields.H"
#include "finiteVolume/fvc/fvcVolumeIntegrate.H"

// * * * * * * * * * * * * * * *  Member Functions * * * * * * * * * * * * * //

template<class RhoFieldType>
void Foam::fv::forceActuationDiskSource::addForceActuationDiskResistance
(
    vectorField& Usource,
    const labelList& cells,
    const scalarField& Vcells,
    const RhoFieldType& rho,
    const vectorField& U
) const
{
    // update force
    const vector force = force_->value(mesh_.time().value());

    // force field in global coordinates
    volVectorField forceSource
    (
        IOobject
        (
            "forceSource",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector(dimless, Zero)
    );

    const vectorField& C = mesh_.C();

    forAll(cells_, i)
    {
        const label celli = cells_[i];
        forceSource[celli] = csys().transform(C[celli], force);
    }

    //- special handling for force with swirl component
    if (
        csys().type() == "cylindrical"
     && radialDistribution_
     && mag(force.y()) > SMALL
    )
    {
        scalar totalWeighting = 0.0;
        // compute total weighting for correct scaling
        forAll(cells_, i)
        {
            const label celli = cells_[i];
            const scalar ri = csys().localPosition(C[celli]).x();
            totalWeighting += ri*Vcells[celli]/V_;
        }
        reduce(totalWeighting, sumOp<scalar>());

        // weight force theta component with radius
        forAll(cells_, i)
        {
            const label celli = cells_[i];
            const scalar ri = csys().localPosition(C[celli]).x();
            vector localForce = force;
            localForce.y() *= ri/totalWeighting;
            forceSource[celli] = csys().transform(C[celli], localForce);
        }
    }

    // apply momentum source
    forAll(cells_, i)
    {
        const label celli = cells_[i];
        Usource[celli] += (Vcells[celli]/VDash_)*forceSource[celli];
    }

    if (debug)
    {
        scalar totalSource(0.0);
        forAll(cells_, i)
        {
            const label celli = cells_[i];
            totalSource += mag((Vcells[celli]/VDash_)*forceSource[celli]);
        }
        reduce(totalSource, sumOp<scalar>());

        Info<< "  total zone volume     : " << V() << nl
             << "  total momentum source : " << totalSource << endl;

        if (mesh_.time().outputTime())
        {
            forceSource.primitiveFieldRef() *= (Vcells/VDash_);
            forceSource.write();
        }
    }
}


// ************************************************************************* //
