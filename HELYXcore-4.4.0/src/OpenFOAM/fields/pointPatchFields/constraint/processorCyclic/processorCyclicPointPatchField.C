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
    (c) 2011-2022 OpenFOAM Foundation
    (c) 2022-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/pointPatchFields/constraint/processorCyclic/processorCyclicPointPatchField.H"
#include "meshes/polyMesh/polyPatches/constraint/processor/processorPolyPatch.H"

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

template<class Type>
Foam::processorCyclicPointPatchField<Type>::processorCyclicPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
:
    coupledPointPatchField<Type>(p, iF),
    procPatch_(refCast<const processorCyclicPointPatch>(p)),
    receiveBuf_(0)
{}


template<class Type>
Foam::processorCyclicPointPatchField<Type>::processorCyclicPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const dictionary& dict
)
:
    coupledPointPatchField<Type>(p, iF, dict),
    procPatch_(refCast<const processorCyclicPointPatch>(p)),
    receiveBuf_(0)
{}


template<class Type>
Foam::processorCyclicPointPatchField<Type>::processorCyclicPointPatchField
(
    const processorCyclicPointPatchField<Type>& ptf,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    coupledPointPatchField<Type>(ptf, p, iF, mapper),
    procPatch_(refCast<const processorCyclicPointPatch>(ptf.patch())),
    receiveBuf_(0)
{}


template<class Type>
Foam::processorCyclicPointPatchField<Type>::processorCyclicPointPatchField
(
    const processorCyclicPointPatchField<Type>& ptf,
    const DimensionedField<Type, pointMesh>& iF
)
:
    coupledPointPatchField<Type>(ptf, iF),
    procPatch_(refCast<const processorCyclicPointPatch>(ptf.patch())),
    receiveBuf_(0)
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class Type>
Foam::processorCyclicPointPatchField<Type>::~processorCyclicPointPatchField()
{}


// ************************************************************************* //
