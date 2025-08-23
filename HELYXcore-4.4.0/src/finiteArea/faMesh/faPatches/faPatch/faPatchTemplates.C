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
    (c) 2018 Engys Ltd.
    (c) 2016-2017 Wikki Ltd

\*----------------------------------------------------------------------------*/

#include "faMesh/faPatches/faPatch/faPatch.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::faPatch::patchInternalField
(
    const UList<Type>& f
) const
{
    tmp<Field<Type>> tpif(new Field<Type>(size()));
    Field<Type>& pif = tpif.ref();

    const labelUList& edgeFaces = this->edgeFaces();

    forAll(pif, facei)
    {
        pif[facei] = f[edgeFaces[facei]];
    }

    return tpif;
}

template<class Type>
void Foam::faPatch::patchInternalField
(
    const UList<Type>& f,
    Field<Type>& pif
) const
{
    pif.setSize(size());

    const labelUList& edgeFaces = this->edgeFaces();

    forAll(pif, facei)
    {
        pif[facei] = f[edgeFaces[facei]];
    }
}


template<class GeometricField, class Type>
const typename GeometricField::Patch& Foam::faPatch::patchField
(
    const GeometricField& gf
) const
{
    return gf.boundaryField()[index()];
}


// ************************************************************************* //
