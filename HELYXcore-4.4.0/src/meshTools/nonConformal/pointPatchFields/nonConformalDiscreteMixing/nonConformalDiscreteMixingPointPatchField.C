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
    (c) 2024-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "nonConformal/pointPatchFields/nonConformalDiscreteMixing/nonConformalDiscreteMixingPointPatchField.H"
#include "fields/GeometricFields/pointFields/pointFields.H"
#include "fields/GeometricFields/pointFields/pointFieldsFwd.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::nonConformalDiscreteMixingPointPatchField<Type>::
nonConformalDiscreteMixingPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
:
    coupledPointPatchField<Type>(p, iF),
    nonConformalPatch_(refCast<const nonConformalDiscreteMixingPointPatch>(p)),
    ppiPtr_(nullptr),
    nbrPpiPtr_(nullptr)
{}


template<class Type>
Foam::nonConformalDiscreteMixingPointPatchField<Type>::
nonConformalDiscreteMixingPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const dictionary& dict
)
:
    coupledPointPatchField<Type>(p, iF, dict),
    nonConformalPatch_(refCast<const nonConformalDiscreteMixingPointPatch>(p)),
    ppiPtr_(nullptr),
    nbrPpiPtr_(nullptr)
{
    if (!isType<nonConformalDiscreteMixingPointPatch>(p))
    {
        FatalIOErrorInFunction(dict)
            << "patch " << this->patch().index() << " not "
            << nonConformalDiscreteMixingPointPatch::typeName << " type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }
}


template<class Type>
Foam::nonConformalDiscreteMixingPointPatchField<Type>::
nonConformalDiscreteMixingPointPatchField
(
    const nonConformalDiscreteMixingPointPatchField<Type>& ptf,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    coupledPointPatchField<Type>(ptf, p, iF, mapper),
    nonConformalPatch_(refCast<const nonConformalDiscreteMixingPointPatch>(p)),
    ppiPtr_(nullptr),
    nbrPpiPtr_(nullptr)
{
    if (!isType<nonConformalDiscreteMixingPointPatch>(this->patch()))
    {
        FatalErrorInFunction
            << "Field type does not correspond to patch type for patch "
            << this->patch().index() << "." << endl
            << "Field type: " << typeName << endl
            << "Patch type: " << this->patch().type()
            << exit(FatalError);
    }
}


template<class Type>
Foam::nonConformalDiscreteMixingPointPatchField<Type>::
nonConformalDiscreteMixingPointPatchField
(
    const nonConformalDiscreteMixingPointPatchField<Type>& ptf,
    const DimensionedField<Type, pointMesh>& iF
)
:
    coupledPointPatchField<Type>(ptf, iF),
    nonConformalPatch_(ptf.nonConformalPatch_),
    ppiPtr_(nullptr),
    nbrPpiPtr_(nullptr)
{}


// ************************************************************************* //
