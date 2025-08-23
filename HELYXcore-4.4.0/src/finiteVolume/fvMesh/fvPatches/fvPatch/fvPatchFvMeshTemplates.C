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

#include "fvMesh/fvPatches/fvPatch/fvPatch.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class GeometricField, class Type>
const typename GeometricField::Patch& Foam::fvPatch::lookupPatchField
(
    const word& name,
    const GeometricField*,
    const Type*
) const
{
    return patchField<GeometricField, Type>
    (
        boundaryMesh().mesh().objectRegistry::template
            lookupObject<GeometricField>(name)
    );
}

template<class GeometricField, class Type>
const typename GeometricField::Patch& Foam::fvPatch::lookupPatchFieldInDb
(
const objectRegistry& obr,
const word& name,
const GeometricField*,
const Type*
) const
{
    return patchField<GeometricField, Type>
    (
        obr.objectRegistry::template
            lookupObject<GeometricField>(name)
    );
}
// ************************************************************************* //
