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
    (c) 2011 OpenFOAM Foundation
    (c) 2010-2019 Engys Ltd.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::mappedFvPatch::interpolate
(
    const Field<Type>& fld,
    const UList<Type>& defaultValues
) const
{
    tmp<Foam::Field<Type>> tfld(new Field<Type>(fld));
    // If direct mapped, map from this side regardless of whether we are owner
    // or neighbour, otherwise gaps are possible
    if (!usingAMI())
    {
        mappedPolyPatch_.distribute(tfld.ref(), defaultValues);
    }
    else
    {
        if (owner())
        {
            mappedPolyPatch_.distribute(tfld.ref(), defaultValues);
        }
        else
        {
            nbrPatch().mappedPolyPatch_.reverseDistribute
            (
                tfld.ref(), defaultValues
            );
        }
    }
    return tfld;
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::mappedFvPatch::interpolate
(
    tmp<Field<Type>> tfld,
    const UList<Type>& defaultValues
) const
{
    return interpolate(tfld(), defaultValues);
}

// ************************************************************************* //
