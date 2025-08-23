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

#include "finiteVolume/fvc/fvcVolumeIntegrate.H"
#include "fvMesh/fvMesh.H"
#include "fields/Fields/Field/Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<Field<Type>>
volumeIntegrate
(
    const VolField<Type>& vf
)
{
    return vf.mesh().V()*vf.primitiveField();
}


template<class Type>
tmp<Field<Type>>
volumeIntegrate
(
    const tmp<VolField<Type>>& tvf
)
{
    tmp<Field<Type>> tvivf = tvf().mesh().V()*tvf().primitiveField();
    tvf.clear();
    return tvivf;
}


template<class Type>
tmp<Field<Type>> volumeIntegrate(const DimensionedField<Type, volMesh>& df)
{
    return df.mesh().V()*df.field();
}


template<class Type>
tmp<Field<Type>>
volumeIntegrate(const tmp<DimensionedField<Type, volMesh>>& tdf)
{
    tmp<Field<Type>> tdidf = tdf().mesh().V()*tdf().field();
    tdf.clear();
    return tdidf;
}


template<class Type>
dimensioned<Type>
domainIntegrate
(
    const VolField<Type>& vf
)
{
    return dimensioned<Type>
    (
        "domainIntegrate(" + vf.name() + ')',
        dimVol*vf.dimensions(),
        gSum(fvc::volumeIntegrate(vf))
    );
}


template<class Type>
dimensioned<Type> domainIntegrate
(
    const tmp<VolField<Type>>& tvf
)
{
    dimensioned<Type> integral = domainIntegrate(tvf());
    tvf.clear();
    return integral;
}


template<class Type>
dimensioned<Type> domainIntegrate
(
    const DimensionedField<Type, volMesh>& df
)
{
    return dimensioned<Type>
    (
        "domainIntegrate(" + df.name() + ')',
        dimVol*df.dimensions(),
        gSum(fvc::volumeIntegrate(df))
    );
}


template<class Type>
dimensioned<Type> domainIntegrate
(
    const tmp<DimensionedField<Type, volMesh>>& tdf
)
{
    dimensioned<Type> integral = domainIntegrate(tdf());
    tdf.clear();
    return integral;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
