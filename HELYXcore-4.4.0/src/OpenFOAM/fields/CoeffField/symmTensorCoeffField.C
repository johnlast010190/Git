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
    (c) 2004-6 H. Jasak All rights reserved

\*---------------------------------------------------------------------------*/

#include "fields/CoeffField/symmTensorCoeffField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CoeffField<Foam::symmTensor>::CoeffField(const label size)
:
    DecoupledCoeffField<symmTensor>(size)
{}


Foam::CoeffField<Foam::symmTensor>::CoeffField
(
    const CoeffField<symmTensor>& f
)
:
    DecoupledCoeffField<symmTensor>(f)
{}


Foam::CoeffField<Foam::symmTensor>::CoeffField
(
    const DecoupledCoeffField<symmTensor>& f
)
:
    DecoupledCoeffField<symmTensor>(f)
{}


Foam::CoeffField<Foam::symmTensor>::CoeffField
(
    const tmp<DecoupledCoeffField<symmTensor>>& tf
)
:
    DecoupledCoeffField<symmTensor>(tf())
{}


Foam::CoeffField<Foam::symmTensor>::CoeffField(Istream& is)
:
    DecoupledCoeffField<symmTensor>(is)
{}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::CoeffField<Foam::symmTensor>::operator=
(
    const CoeffField<symmTensor>& f
)
{
    DecoupledCoeffField<symmTensor>::operator=(f);
}


void Foam::CoeffField<Foam::symmTensor>::operator=
(
    const tmp<CoeffField<symmTensor>>& f
)
{
    DecoupledCoeffField<symmTensor>::operator=(f);
}


void Foam::CoeffField<Foam::symmTensor>::operator=
(
    const CoeffField<symmTensor>::scalarTypeField& f
)
{
    DecoupledCoeffField<symmTensor>::operator=(f);
}


void Foam::CoeffField<Foam::symmTensor>::operator=
(
    const tmp<CoeffField<symmTensor>::scalarTypeField>& f
)
{
    DecoupledCoeffField<symmTensor>::operator=(f);
}


void Foam::CoeffField<Foam::symmTensor>::operator=
(
    const CoeffField<symmTensor>::linearTypeField& f
)
{
    DecoupledCoeffField<symmTensor>::operator=(f);
}


void Foam::CoeffField<Foam::symmTensor>::operator=
(
    const tmp<CoeffField<symmTensor>::linearTypeField>& f
)
{
    DecoupledCoeffField<symmTensor>::operator=(f);
}


/* * * * * * * * * * * * * * * * Global functions  * * * * * * * * * * * * * */

Foam::tmp<Foam::CoeffField<Foam::symmTensor>> Foam::inv
(
    const CoeffField<symmTensor>& f
)
{
    const DecoupledCoeffField<symmTensor>& df = f;

    return tmp<CoeffField<symmTensor>>
    (
        new CoeffField<symmTensor>(inv(df)())
    );
}


// ************************************************************************* //
