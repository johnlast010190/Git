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

#include "AMIInterpolation/patches/cyclicAMI/cyclicAMIPointPatchField/cyclicAMIPointPatchField.H"
#include "fields/Fields/transformField/transformField.H"
#include "fields/GeometricFields/pointFields/pointFields.H"
#include "fields/GeometricFields/pointFields/pointFieldsFwd.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::cyclicAMIPointPatchField<Type>::cyclicAMIPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
:
    coupledPointPatchField<Type>(p, iF),
    cyclicAMIPatch_(refCast<const cyclicAMIPointPatch>(p)),
    ppiPtr_(nullptr),
    nbrPpiPtr_(nullptr)
{}


template<class Type>
Foam::cyclicAMIPointPatchField<Type>::cyclicAMIPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const dictionary& dict
)
:
    coupledPointPatchField<Type>(p, iF, dict),
    cyclicAMIPatch_(refCast<const cyclicAMIPointPatch>(p)),
    ppiPtr_(nullptr),
    nbrPpiPtr_(nullptr)
{
    if (!isType<cyclicAMIPointPatch>(p))
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "patch " << this->patch().index() << " not cyclicAMI type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }
}


template<class Type>
Foam::cyclicAMIPointPatchField<Type>::cyclicAMIPointPatchField
(
    const cyclicAMIPointPatchField<Type>& ptf,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    coupledPointPatchField<Type>(ptf, p, iF, mapper),
    cyclicAMIPatch_(refCast<const cyclicAMIPointPatch>(p)),
    ppiPtr_(nullptr),
    nbrPpiPtr_(nullptr)
{
    if (!isType<cyclicAMIPointPatch>(this->patch()))
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
Foam::cyclicAMIPointPatchField<Type>::cyclicAMIPointPatchField
(
    const cyclicAMIPointPatchField<Type>& ptf,
    const DimensionedField<Type, pointMesh>& iF
)
:
    coupledPointPatchField<Type>(ptf, iF),
    cyclicAMIPatch_(ptf.cyclicAMIPatch_),
    ppiPtr_(nullptr),
    nbrPpiPtr_(nullptr)
{}


// ************************************************************************* //
