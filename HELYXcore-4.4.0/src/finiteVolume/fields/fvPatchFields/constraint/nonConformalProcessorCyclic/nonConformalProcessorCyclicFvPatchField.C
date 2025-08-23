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

#include "fields/fvPatchFields/constraint/nonConformalProcessorCyclic/nonConformalProcessorCyclicFvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

template<class Type>
Foam::nonConformalProcessorCyclicFvPatchField<Type>::
nonConformalProcessorCyclicFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    processorCyclicFvPatchField<Type>(p, iF),
    procPatch_(refCast<const nonConformalProcessorCyclicFvPatch>(p))
{}


template<class Type>
Foam::nonConformalProcessorCyclicFvPatchField<Type>::
nonConformalProcessorCyclicFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    processorCyclicFvPatchField<Type>(p, iF, dict),
    procPatch_(refCast<const nonConformalProcessorCyclicFvPatch>(p))
{}


template<class Type>
Foam::nonConformalProcessorCyclicFvPatchField<Type>::
nonConformalProcessorCyclicFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const Field<Type>& f
)
:
    processorCyclicFvPatchField<Type>(p, iF, f),
    procPatch_(refCast<const nonConformalProcessorCyclicFvPatch>(p))
{}


template<class Type>
Foam::nonConformalProcessorCyclicFvPatchField<Type>::
nonConformalProcessorCyclicFvPatchField
(
    const nonConformalProcessorCyclicFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    processorCyclicFvPatchField<Type>(ptf, p, iF, mapper),
    procPatch_(refCast<const nonConformalProcessorCyclicFvPatch>(p))
{}


template<class Type>
Foam::nonConformalProcessorCyclicFvPatchField<Type>::
nonConformalProcessorCyclicFvPatchField
(
    const nonConformalProcessorCyclicFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    processorCyclicFvPatchField<Type>(ptf, iF),
    procPatch_(refCast<const nonConformalProcessorCyclicFvPatch>(ptf.patch()))
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class Type>
Foam::nonConformalProcessorCyclicFvPatchField<Type>::
~nonConformalProcessorCyclicFvPatchField()
{}


// ************************************************************************* //
