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
    (c) 2021-2022 OpenFOAM Foundation
    (c) 2022 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/constraint/nonConformalCyclic/nonConformalCyclicFvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::nonConformalCyclicFvPatchField<Type>::nonConformalCyclicFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    cyclicFvPatchField<Type>(p, iF),
    nonConformalCyclicFvPatch_(refCast<const nonConformalCyclicFvPatch>(p))
{}


template<class Type>
Foam::nonConformalCyclicFvPatchField<Type>::nonConformalCyclicFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    cyclicFvPatchField<Type>(p, iF, dict),
    nonConformalCyclicFvPatch_(refCast<const nonConformalCyclicFvPatch>(p))
{}


template<class Type>
Foam::nonConformalCyclicFvPatchField<Type>::nonConformalCyclicFvPatchField
(
    const nonConformalCyclicFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    cyclicFvPatchField<Type>(ptf, p, iF, mapper),
    nonConformalCyclicFvPatch_(refCast<const nonConformalCyclicFvPatch>(p))
{}


template<class Type>
Foam::nonConformalCyclicFvPatchField<Type>::nonConformalCyclicFvPatchField
(
    const nonConformalCyclicFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    cyclicFvPatchField<Type>(ptf, iF),
    nonConformalCyclicFvPatch_(ptf.nonConformalCyclicFvPatch_)
{}


// ************************************************************************* //
