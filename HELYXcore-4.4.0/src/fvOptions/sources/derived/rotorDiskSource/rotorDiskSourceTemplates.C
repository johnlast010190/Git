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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2024 Engys Ltd.
    (c) 2024 Alexandre Capitao Patrao (Chalmers) [3]

\*---------------------------------------------------------------------------*/

#include "rotorDiskSource.H"
#include "fields/volFields/volFields.H"
#include "global/unitConversion/unitConversion.H"

using namespace Foam::constant;

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class RhoFieldType>
void Foam::fv::rotorDiskSource::calculate
(
    const RhoFieldType& rho,
    const vectorField& U,
    const scalarField& thetag,
    vectorField& force,
    const bool divideVolume,
    const bool output
) const
{
    const scalarField& V = mesh_.V();

    // Logging info
    scalar dragEff = 0.0;
    scalar liftEff = 0.0;
    scalar thrustEff = 0.0;
    scalar torqueEff = 0.0;
    scalar powerEff = 0.0;
    scalar AOAmin = GREAT;
    scalar AOAmax = -GREAT;
    scalar epsMin = GREAT;
    scalar epsMax = -GREAT;

    // Cached position-dependent rotations available?
    const bool hasCache = Rcyl_.valid();

    forAll(cells_, i)
    {
        if (area_[i] > ROOTVSMALL)
        {
            const label celli = cells_[i];

            const scalar radius = x_[i].x();

            const tensor Rcyl =
            (
                hasCache
              ? (*Rcyl_)[i]
              : coordSys_.R(mesh_.C()[celli])
            );

            // Transform velocity into local cylindrical reference frame
            vector Uc = invTransform(Rcyl, U[celli]);

            // Transform velocity into local coning system
            Uc = transform(Rcone_[i], Uc);

            // Set radial component of velocity to zero
            Uc.x() = 0.0;

            // Set blade normal component of velocity
            Uc.y() = radius*omega_ - Uc.y();

            // Determine blade data for this radius
            // i2 = index of upper radius bound data point in blade list
            scalar twist = 0.0;
            scalar chord = 0.0;
            label i1 = -1;
            label i2 = -1;
            scalar invDr = 0.0;
            blade_.interpolate(radius, twist, chord, i1, i2, invDr);

            // Flip geometric angle if blade is spinning in reverse (clockwise)
            scalar alphaGeom = thetag[i] + twist;
            if (omega_ < 0)
            {
                alphaGeom = mathematical::pi - alphaGeom;
            }

            // Effective angle of attack
            scalar alphaEff = alphaGeom - atan2(-Uc.z(), Uc.y());
            if (alphaEff > mathematical::pi)
            {
                alphaEff -= mathematical::twoPi;
            }
            if (alphaEff < -mathematical::pi)
            {
                alphaEff += mathematical::twoPi;
            }

            AOAmin = min(AOAmin, alphaEff);
            AOAmax = max(AOAmax, alphaEff);

            // Determine profile data for this radius and angle of attack
            const label profile1 = blade_.profileID()[i1];
            const label profile2 = blade_.profileID()[i2];

            scalar Cd1 = 0.0;
            scalar Cl1 = 0.0;
            profiles_[profile1].Cdl(alphaEff, Cd1, Cl1);

            scalar Cd2 = 0.0;
            scalar Cl2 = 0.0;
            profiles_[profile2].Cdl(alphaEff, Cd2, Cl2);

            scalar Cd = invDr*(Cd2 - Cd1) + Cd1;
            scalar Cl = invDr*(Cl2 - Cl1) + Cl1;

            // Apply tip effect for blade lift
            scalar tipFactor = 0.0;
            if (!prandtlDrella_)
            {
                tipFactor = neg(radius/rMax_ - tipEffect_);
            }

            // Calculate forces perpendicular to blade
            scalar pDyn = 0.5*rho[celli]*magSqr(Uc);

            scalar f = pDyn*chord*nBlades_*area_[i]/radius/mathematical::twoPi;
            vector localForce(Zero);
            if (!prandtlDrella_)
            {
                localForce = vector(0.0, -f*Cd, tipFactor*f*Cl);
            }
            else
            {
                // Flow angle
                scalar eps = atan2(-Uc.z(), Uc.y());
                if (eps < -mathematical::pi)
                {
                    eps = (2.0*mathematical::pi + eps);
                }
                if (eps > mathematical::pi)
                {
                    eps = (eps - 2.0*mathematical::pi);
                }
                epsMin = min(epsMin, eps);
                epsMax = max(epsMax, eps);

                //Drela tip factor:
                //scalar lambdaVal = radius/(0.5*diameterRef_)*tan(eps);
                //
                // The original implementation does not have mag. But if eps
                // becomes negative the acos(exp(-1/tan(x))) has discontinuity
                // at zero causing it to crash.

                scalar nRadius = radius/(0.5*diameterRef_);

                //- additional check if a cell radius is larger than half of
                //  the specified diameter
                nRadius = min(nRadius, 1);

                scalar lambdaVal = nRadius*mag(tan(eps));

                scalar tipFactor_f = 0;

                if (mag(lambdaVal)<SMALL)
                {
                    lambdaVal = 0.0;
                    tipFactor_f = 0.0;
                    tipFactor = 0.0;
                }
                else
                {
                    tipFactor_f = (nBlades_/2)*(1-nRadius)*(1/lambdaVal);
                    tipFactor = 2/mathematical::pi*acos(exp(-tipFactor_f));
                }

                // Tangential and axial forces
                scalar fTang = (f*Cd*cos(eps) + tipFactor*f*Cl*sin(eps));
                scalar fAxial = (-f*Cd*sin(eps) + tipFactor*f*Cl*cos(eps));

                thrustEff += rhoRef_*fAxial;
                torqueEff += rhoRef_*fTang*radius;
                powerEff += rhoRef_*fTang*radius*omega_;

                localForce = vector(0.0,-fTang,fAxial);

            }
            // Accumulate forces
            dragEff += rhoRef_*localForce.y();
            liftEff += rhoRef_*localForce.z();

            // Transform force from local coning system into rotor
            // cylindrical
            localForce = invTransform(Rcone_[i], localForce);

            // Transform force into global Cartesian coordinate system
            force[celli] = transform(Rcyl, localForce);

            if (divideVolume)
            {
                force[celli] /= V[celli];
            }
        }
    }

    if (output)
    {
        reduce
        (
            std::tie(AOAmin, AOAmax, dragEff, liftEff),
            ParallelOp
            <
                minOp<scalar>, maxOp<scalar>,
                sumOp<scalar>, sumOp<scalar>
            >{}
        );

        if (prandtlDrella_)
        {
            reduce
            (
                std::tie(epsMin, epsMax, torqueEff, powerEff, thrustEff),
                ParallelOp
                <
                    minOp<scalar>, maxOp<scalar>,
                    sumOp<scalar>, sumOp<scalar>, sumOp<scalar>
                >{}
            );

            scalar etaProp = thrustEff*refVelEta_/powerEff;

            Info<< type() << " output:" << nl
                 << tab << "min/max(AOA) = " << radToDeg(AOAmin) << ", "
                 << radToDeg(AOAmax) << nl
                 << tab << "min/max(eps) = " << radToDeg(epsMin) << ", "
                 << radToDeg(epsMax) << nl
                 << tab << "Rotor thrust = " << thrustEff << nl
                 << tab << "Rotor torque = " << torqueEff << nl
                 << tab << "Rotor power = " << powerEff << nl
                 << tab << "Rotor propeller efficiency = " << etaProp << nl
                 << tab << "Effective drag = " << dragEff << nl
                 << tab << "Effective lift = " << liftEff << endl;
        }
        else
        {
            Info<< type() << " output:" << nl
                 << tab << "min/max(AOA)   = " << radToDeg(AOAmin) << ", "
                 << radToDeg(AOAmax) << nl
                 << tab << "Effective drag = " << dragEff << nl
                 << tab << "Effective lift = " << liftEff
                 << endl;
        }
    }
}


template<class Type>
void Foam::fv::rotorDiskSource::writeField
(
    const word& name,
    const List<Type>& values,
    const bool writeNow
) const
{
    typedef VolField<Type> fieldType;

    if (mesh_.time().writeTime() || writeNow)
    {
        tmp<fieldType> tfield
        (
            new fieldType
            (
                IOobject
                (
                    name,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensioned<Type>("zero", dimless, Zero)
            )
        );

        Field<Type>& field = tfield.ref().primitiveFieldRef();

        if (cells_.size() != values.size())
        {
            FatalErrorInFunction
                << abort(FatalError);
        }

        forAll(cells_, i)
        {
            const label celli = cells_[i];
            field[celli] = values[i];
        }

        tfield().write();
    }
}


// ************************************************************************* //
