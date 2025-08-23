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
    (c) 2019 OpenFOAM Foundation
    (c) 2022 Engys Ltd.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::distributedWeightedFvPatchFieldMapper::map
(
    Field<Type>& f,
    const Field<Type>& mapF,
    const bool applyFlip
) const
{
    if (distributed())
    {
        // Fetch remote parts of mapF
        const distributionMapBase& distMap = *distMapPtr_;
        Field<Type> newMapF(mapF);

        if (applyFlip)
        {
            distMap.distribute(newMapF);
        }
        else
        {
            distMap.distribute(newMapF, noOp());
        }

        f.map(newMapF, addressing(), weights());
    }
    else
    {
        f.map(mapF, addressing(), weights());
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::distributedWeightedFvPatchFieldMapper::map
(
    const Field<Type>& mapF,
    const bool applyFlip
) const
{
    tmp<Field<Type>> tf(new Field<Type>(size()));
    map(tf.ref(), mapF, applyFlip);

    return tf;
}


// ************************************************************************* //
