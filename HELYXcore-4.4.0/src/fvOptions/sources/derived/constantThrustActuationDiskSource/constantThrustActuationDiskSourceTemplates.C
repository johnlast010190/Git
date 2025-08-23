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
    (c) 2013 Kevin Maki
    (c) 2013 Engys Ltd.

\*----------------------------------------------------------------------------*/

#include "constantThrustActuationDiskSource.H"
#include "fields/volFields/volFields.H"
#include "finiteVolume/fvc/fvcVolumeIntegrate.H"

// * * * * * * * * * * * * * * *  Member Functions * * * * * * * * * * * * * //

template<class RhoFieldType>
void Foam::fv::constantThrustActuationDiskSource::addActuationDiskAxialInertialResistance
(
    vectorField& Usource,
    const labelList& cells,
    const scalarField& Vcells,
    const RhoFieldType& rho,
    const vectorField& U
) const
{
    // update thrust
    const_cast<constantThrustActuationDiskSource*>(this)
        ->updateThrust(TF1_->value(mesh_.time().value()));

    vector uniDiskDir = diskDir_/mag(diskDir_);

    volVectorField direction
    (
        IOobject
        (
            "direction",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector(dimless, Zero)
    );

    const vectorField& C = mesh_.C();

    forAll(cells, i)
    {
        const label celli = cells[i];

        // position in thruster co-ordinate system
        point x = csys().localPosition(C[celli]);

        // compute local transformation of uniDiskDir
        const scalar scaling = max(100*mag(x),1e-4);
        point x2 = x + uniDiskDir/scaling;
        point x2g = csys().globalPosition(x2);
        vector dir = x2g - C[celli];

        // scale direction vector back
        direction[celli] = dir*scaling;
    }

    // ensure global magnitude
    scalar dirScaling = (fvc::domainIntegrate(mag(direction))/V()).value();
    direction /= dirScaling;

    // apply momentum source
    forAll(cells, i)
    {
        const label celli = cells[i];
        Usource[celli] += (Vcells[celli]/VDash_)*T()*direction[celli];
    }

    if (debug)
    {
        scalar totalSource(0.0);
        forAll(cells, i)
        {
            const label celli = cells[i];
            totalSource += mag((Vcells[celli]/VDash_)*T()*direction[celli]);
        }
        reduce(totalSource, sumOp<scalar>());

        Info<< "  total thruster zone volume     : " << V() << nl
             << "  total thruster momentum source : " << totalSource << nl
             << "  mean direction scaling factor  : " << dirScaling << endl;

        if (mesh_.time().outputTime())
        {
            volVectorField thrusterSource
            (
                IOobject
                (
                    "thrusterSource",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedVector(dimless, Zero)
            );
            thrusterSource.primitiveFieldRef() = (Vcells/VDash_)*T()*direction;
            thrusterSource.write();
            direction.write();
        }
    }
}


// ************************************************************************* //
