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

inline Foam::scalar Foam::helyxSample::interp
(
    const word& name,
    const std::vector<label>& nb,
    const std::vector<scalar>& dist
)
{
    scalar sum1 = 0;
    scalar sum2 = 0;
    scalar eps = 1.0e-25;

    std::vector<scalar> prop(nb.size());

    for (unsigned i = 0; i < nb.size(); i++)
    {
        label ic = nb[i];
        if (!Pstream::parRun())
        {
            prop[i] = source().scalarFields_[name][ic];
        }
        else
        {
            prop[i] = sourceScalarFields_[name][ic];
        }
    }

    for (unsigned i = 0; i < nb.size(); i++)
    {
        scalar disti = dist[i] + eps;

        sum1 += prop[i]/disti;
        sum2 += 1.0/disti;
    }

    return sum1/sum2;
}


inline Foam::point Foam::helyxSample::vInterp
(
    const word& name,
    const std::vector<label>& nb,
    const std::vector<scalar>& dist
)
{
    point sum1(0, 0, 0);
    scalar sum2 = 0;
    scalar eps = 1.0e-25;

    std::vector<point> prop(nb.size());

    for (unsigned i = 0; i < nb.size(); i++)
    {
        label ic = nb[i];
        if (!Pstream::parRun())
        {
            prop[i] = source().vectorFields_[name][ic];
        }
        else
        {
            prop[i] = sourceVectorFields_[name][ic];
        }
    }

    for (unsigned i = 0; i < nb.size(); i++)
    {
        scalar disti = dist[i] + eps;

        sum1 += prop[i]/disti;
        sum2 += 1.0/disti;
    }

    return sum1/sum2;
}


inline void Foam::helyxSample::outDomainGrid
(
    gridField& grid,
    label i,
    label j,
    label k
)
{
    label nsc = grid.scalarSamples_.size();

    if (nsc > 0)
    {
        gridField::nodeField node(nsc, 1.0e20);
        grid.gridFields_(i, j, k) = node;
    }

    label nvc = grid.vectorSamples_.size();

    if (nvc > 0)
    {
        point bigVector(1.0e20, 1.0e20, 1.0e20);

        gridField::nodeVectors vnode(nvc, bigVector);
        grid.gridVectors_(i, j, k) = vnode;
    }
}


inline void Foam::helyxSample::getGridField
(
    gridField& grid,
    label i,
    label j,
    label k,
    const std::vector<label>& nb,
    const std::vector<scalar>& dist,
    bool insolid
)
{
    label nsc = grid.scalarSamples_.size();

    grid.dists_(i, j, k) = dist[0];

    if (nsc > 0)
    {
        gridField::nodeField node(nsc, 0.0);

        if (insolid)
        {
            grid.gridFields_(i, j, k) = node;
        }
        else
        {
            for (label isc = 0; isc < nsc; isc++)
            {
                word name = grid.scalarSamples_[isc];
                label ic = nb[0];

                if (!grid.interp_)
                {
                    scalar sval;

                    if (!Pstream::parRun())
                    {
                        sval = source().scalarFields_[name][ic];
                    }
                    else
                    {
                        sval = sourceScalarFields_[name][ic];
                    }

                    node[isc] = sval;
                }
                else
                {
                    scalar sval = interp(name, nb, dist);
                    node[isc] = sval;
                }
            }

            grid.gridFields_(i, j, k) = node;
        }
    }
}


inline void Foam::helyxSample::getGridVectorField
(
    gridField& grid,
    label i,
    label j,
    label k,
    const std::vector<label>& nb,
    const std::vector<scalar>& dist,
    bool insolid
)
{
    label nvc = grid.vectorSamples_.size();

    if (nvc > 0)
    {
        point zeroVector(0, 0, 0);
        gridField::nodeVectors node(nvc, zeroVector);

        if (insolid)
        {
            grid.gridVectors_(i, j, k) = node;
        }
        else
        {
            for (label ivc = 0; ivc < nvc; ivc++)
            {
                word name = grid.vectorSamples_[ivc];
                label ic = nb[0];

                point svec;
                if (!grid.interp_)
                {
                    if (!Pstream::parRun())
                    {
                        svec = source().vectorFields_[name][ic];
                    }
                    else
                    {
                        svec = sourceVectorFields_[name][ic];
                    }

                    node[ivc] = svec;
                }
                else
                {
                    point svec = vInterp(name, nb, dist);
                    node[ivc] = svec;
                }
            }

            grid.gridVectors_(i, j, k) = node;
        }
    }
}


// ************************************************************************* //
