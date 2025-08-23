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
    (c) 1991-2005 OpenCFD Ltd.
    (c) 2016 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/faPatchFields/basic/coupled/coupledFaPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
coupledFaePatchField<Type>::coupledFaePatchField
(
    const faPatch& p,
    const DimensionedField<Type, faEdgeMesh>& iF
)
:
    faePatchField<Type>(p, iF)
{}


template<class Type>
coupledFaePatchField<Type>::coupledFaePatchField
(
    const faPatch& p,
    const DimensionedField<Type, faEdgeMesh>& iF,
    const Field<Type>& f
)
:
    faePatchField<Type>(p, iF, f)
{}


template<class Type>
coupledFaePatchField<Type>::coupledFaePatchField
(
    const coupledFaePatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, faEdgeMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    faePatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
coupledFaePatchField<Type>::coupledFaePatchField
(
    const faPatch& p,
    const DimensionedField<Type, faEdgeMesh>& iF,
    const dictionary& dict
)
:
    faePatchField<Type>(p, iF, dict)
{}


template<class Type>
coupledFaePatchField<Type>::coupledFaePatchField
(
    const coupledFaePatchField<Type>& ptf
)
:
    faePatchField<Type>(ptf)
{}


template<class Type>
coupledFaePatchField<Type>::coupledFaePatchField
(
    const coupledFaePatchField<Type>& ptf,
    const DimensionedField<Type, faEdgeMesh>& iF
)
:
    faePatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
