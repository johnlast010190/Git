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
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::nonConformalDiscreteMixingFvPatch::interpolate
(
    const Field<Type>& fld
) const
{
    if (!returnReduce(size(), sumOp<label>()))
    {
        return tmp<Field<Type>>(new Field<Type>());
    }

    tmp<Field<Type>> tresult(new Field<Type>(size(), Zero));

    if (owner())
    {
        if (fld.size() != nbrPatch().size())
        {
            FatalErrorInFunction
                << "Supplied field size is not equal to neighbour patch size:"
                << nl << "    owner patch     = " << size()
                << nl << "    neighbour patch = " << nbrPatch().size()
                << nl << "    supplied field  = " << fld.size()
                << abort(FatalError);
        }

        tresult.ref() = intersection().interpolateToOwner(fld);
    }
    else
    {
        if (fld.size() != nbrPatch().size())
        {
            FatalErrorInFunction
                << "Supplied field size is not equal to owner patch size:"
                << nl << "    owner patch     = " << nbrPatch().size()
                << nl << "    neighbour patch = " << size()
                << nl << "    supplied field  = " << fld.size()
                << abort(FatalError);
        }

        tresult.ref() = nbrPatch().intersection().interpolateToNeighbour(fld);
    }

    return tresult;
}


template<class Type, class CombineOp>
void Foam::nonConformalDiscreteMixingFvPatch::interpolate
(
    const UList<Type>& fld,
    List<Type>& result,
    const CombineOp& cop
) const
{
    if (owner())
    {
        intersection().interpolateToOwner(fld, result, cop);
    }
    else
    {
        nbrPatch().intersection().interpolateToNeighbour(fld, result, cop);
    }
}


// ************************************************************************* //
