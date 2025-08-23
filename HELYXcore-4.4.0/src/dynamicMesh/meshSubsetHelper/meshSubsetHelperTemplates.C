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
    (c) 2016 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "meshSubsetHelper/meshSubsetHelper.H"
#include "fvMesh/fvMesh.H"
#include "fields/volFields/volFields.H"
#include "meshes/polyMesh/globalMeshData/globalIndex.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchField.H"


// * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::VolField<Type>>
Foam::meshSubsetHelper::zeroGradientField
(
    const typename GeometricField
    <
        Type,
        fvPatchField,
        volMesh
    >::Internal& df
)
{
    IOobject io(df);
    io.readOpt()  = IOobject::NO_READ;
    io.writeOpt() = IOobject::NO_WRITE;
    io.registerObject() = false;

    tmp<VolField<Type>> tvf
    (
        new VolField<Type>
        (
            io,
            df.mesh(),
            dimensioned<Type>("0", df.dimensions(), Zero),
            zeroGradientFvPatchField<Type>::typeName
        )
    );
    tvf.ref().primitiveFieldRef() = df;
    tvf.ref().correctBoundaryConditions();

    return tvf;
}


template<class Type>
Foam::tmp<Foam::VolField<Type>>
Foam::meshSubsetHelper::interpolate
(
    const VolField<Type>& vf
) const
{
    if (subsetter_.hasSubMesh())
    {
        tmp<VolField<Type>> tfld
        (
            subsetter_.interpolate(vf, true)
        );
        tfld.ref().checkOut();
        tfld.ref().rename(vf.name());
        return tfld;
    }
    else
    {
        return vf;
    }
}


template<class Type>
Foam::tmp<Foam::VolField<Type>>
Foam::meshSubsetHelper::interpolate
(
    const typename GeometricField
    <
        Type,
        fvPatchField,
        volMesh
    >::Internal& df
) const
{
    tmp<VolField<Type>> tvf =
        zeroGradientField<Type>(df);

    if (subsetter_.hasSubMesh())
    {
        return interpolate<Type>(tvf());
    }
    else
    {
        return tvf;
    }
}


template<class GeoField>
Foam::tmp<GeoField>
Foam::meshSubsetHelper::interpolate
(
    const GeoField& fld
) const
{
    if (subsetter_.hasSubMesh())
    {
        tmp<GeoField> subFld = subsetter_.interpolate(fld, true);
        subFld.ref().checkOut();
        subFld.ref().rename(fld.name());
        return subFld;
    }
    else
    {
        return fld;
    }
}


// ************************************************************************* //
