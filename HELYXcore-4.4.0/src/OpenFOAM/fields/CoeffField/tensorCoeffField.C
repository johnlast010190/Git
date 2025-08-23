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

#include "fields/CoeffField/tensorCoeffField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CoeffField<Foam::tensor>::CoeffField(const label size)
:
    DecoupledCoeffField<tensor>(size)
{}


Foam::CoeffField<Foam::tensor>::CoeffField(const CoeffField<tensor>& f)
:
    DecoupledCoeffField<tensor>(f)
{}


Foam::CoeffField<Foam::tensor>::CoeffField(const DecoupledCoeffField<tensor>& f)
:
    DecoupledCoeffField<tensor>(f)
{}


Foam::CoeffField<Foam::tensor>::CoeffField
(
    const tmp<DecoupledCoeffField<tensor>>& tf
)
:
    DecoupledCoeffField<tensor>(tf())
{}


Foam::CoeffField<Foam::tensor>::CoeffField(Istream& is)
:
    DecoupledCoeffField<tensor>(is)
{}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::CoeffField<Foam::tensor>::operator=(const CoeffField<tensor>& f)
{
    DecoupledCoeffField<tensor>::operator=(f);
}


void Foam::CoeffField<Foam::tensor>::operator=
(
    const tmp<CoeffField<tensor>>& f
)
{
    DecoupledCoeffField<tensor>::operator=(f);
}


void Foam::CoeffField<Foam::tensor>::operator=
(
    const CoeffField<tensor>::scalarTypeField& f
)
{
    DecoupledCoeffField<tensor>::operator=(f);
}


void Foam::CoeffField<Foam::tensor>::operator=
(
    const tmp<CoeffField<tensor>::scalarTypeField>& f
)
{
    DecoupledCoeffField<tensor>::operator=(f);
}


void Foam::CoeffField<Foam::tensor>::operator=
(
    const CoeffField<tensor>::linearTypeField& f
)
{
    DecoupledCoeffField<tensor>::operator=(f);
}


void Foam::CoeffField<Foam::tensor>::operator=
(
    const tmp<CoeffField<tensor>::linearTypeField>& f
)
{
    DecoupledCoeffField<tensor>::operator=(f);
}


/* * * * * * * * * * * * * * * * Global functions  * * * * * * * * * * * * * */

Foam::tmp<Foam::CoeffField<Foam::tensor>> Foam::inv
(
    const CoeffField<tensor>& f
)
{
    const DecoupledCoeffField<tensor>& df = f;

    return tmp<CoeffField<tensor>>(new CoeffField<tensor>(inv(df)()));
}


// ************************************************************************* //
