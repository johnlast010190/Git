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
    (c) 2019-2025 Engys Ltd.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline void Foam::fieldSample::getSolidField
(
    label i,
    label j,
    label k,
    const std::vector<label>& nb,
    const std::vector<scalar>& dist,
    bool setLevel
)
{
    if (!setLevel)
    {
        gridField_(i, j, k) = 0;
    }
    else
    {
        gridField_(i, j, k) = -dist[0];
    }
}


inline void Foam::fieldSample::getGridField
(
    label i,
    label j,
    label k,
    const std::vector<label>& nb,
    const std::vector<scalar>& dist
)
{
    label ic = nb[0];

    scalar node = 0;

    if (!interp_)
    {
        node = (*field_)[ic];
    }
    else
    {
        scalar sum = 0;
        scalar sum1 = 1.0e-30;

        label sz = nb.size();

        for (label ik = 0; ik < sz; ik++)
        {
            label ick = nb[ik];

            scalar tmp = 1.0/(dist[ik] + 1.0e-30);

            sum += ((*field_)[ick]*tmp);
            sum1 += tmp;
        }

        node = sum/sum1;
    }

    gridField_(i, j, k) = node;
}


// ************************************************************************* //
