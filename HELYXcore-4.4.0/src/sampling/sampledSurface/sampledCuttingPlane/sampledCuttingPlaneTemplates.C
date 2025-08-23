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

\*---------------------------------------------------------------------------*/

#include "sampledSurface/sampledCuttingPlane/sampledCuttingPlane.H"
#include "fields/volFields/volFieldsFwd.H"
#include "fields/GeometricFields/pointFields/pointFields.H"
#include "interpolation/volPointInterpolation/volPointInterpolation.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::sampledCuttingPlane::sampleField
(
    const VolField<Type>& vField
) const
{
    return tmp<Field<Type>>(new Field<Type>(vField, surface().meshCells()));
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::sampledCuttingPlane::interpolateField
(
    const interpolation<Type>& interpolator
) const
{
    // Get fields to sample. Assume volPointInterpolation!
    const VolField<Type>& volFld =
        interpolator.psi();

    if (subMeshPtr_.valid())
    {
        tmp<VolField<Type>> tvolSubFld =
            subMeshPtr_().interpolate(volFld);

        const VolField<Type>& volSubFld =
            tvolSubFld();

        tmp<PointField<Type>> tpointSubFld =
            volPointInterpolation::New(volSubFld.mesh()).interpolate(volSubFld);

        // Sample.
        return surface().interpolate
        (
            (
                average_
              ? pointAverage(tpointSubFld())()
              : volSubFld
            ),
            tpointSubFld()
        );
    }
    else
    {
        tmp<PointField<Type>> tpointFld
        (
            volPointInterpolation::New(volFld.mesh()).interpolate(volFld)
        );

        // Sample.
        return surface().interpolate
        (
            (
                average_
              ? pointAverage(tpointFld())()
              : volFld
            ),
            tpointFld()
        );
    }
}


// ************************************************************************* //
