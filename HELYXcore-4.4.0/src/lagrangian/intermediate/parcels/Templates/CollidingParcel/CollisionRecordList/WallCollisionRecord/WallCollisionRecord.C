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

#include "parcels/Templates/CollidingParcel/CollisionRecordList/WallCollisionRecord/WallCollisionRecord.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class Type>
const Foam::scalar Foam::WallCollisionRecord<Type>::errorCosAngle(1.0 + 1e-6);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::WallCollisionRecord<Type>::WallCollisionRecord()
:
    accessed_(false),
    pRel_(),
    data_(Zero)
{}


template<class Type>
Foam::WallCollisionRecord<Type>::WallCollisionRecord
(
    bool accessed,
    const vector& pRel,
    const Type& data
)
:
    accessed_(accessed),
    pRel_(pRel),
    data_(data)
{}


template<class Type>
Foam::WallCollisionRecord<Type>::WallCollisionRecord
(
    const WallCollisionRecord<Type>& wCR
)
:
    accessed_(wCR.accessed_),
    pRel_(wCR.pRel_),
    data_(wCR.data_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::WallCollisionRecord<Type>::~WallCollisionRecord()
{}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

template<class Type>
void Foam::WallCollisionRecord<Type>::operator=
(
    const WallCollisionRecord<Type>& rhs
)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorInFunction
            << "Attempted assignment to self"
            << abort(FatalError);
    }

    accessed_ = rhs.accessed_;
    pRel_ = rhs.pRel_;
    data_ = rhs.data_;
}


// * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * * * //

#include "parcels/Templates/CollidingParcel/CollisionRecordList/WallCollisionRecord/WallCollisionRecordIO.C"


// ************************************************************************* //
