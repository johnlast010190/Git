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
    (c) 2014 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "hessianMeshOptimization/smoother/generalMeshSmoother.H"


namespace Foam
{
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
generalMeshSmoother::generalMeshSmoother(const meshMetric& metric)
:
    mesh_(metric.getMesh()),
    metric_(metric)
{
}

generalMeshSmoother::~generalMeshSmoother(){}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
vector generalMeshSmoother::displacement(const label& pI) const
{
    const scalarField& anisotropy = metric_.metric().aspectRatio();

    const labelList& pCells = mesh_.pointCells()[pI];

    vector centroid = vector::zero;

    forAll(pCells, cI)
    {
        if (anisotropy[pCells[cI]]>1)
        {
            return vector::zero;
        }
    }

    forAll(pCells, cI)
    {
        const label& cL = pCells[cI];
        centroid += metric_.cellCtrs()[cL];
    }

    centroid /= pCells.size();

    vector disp = 0.15*(metric_.getPoints()[pI]-centroid);

    return disp;
}

// ************************************************************************* //

} /* namespace Foam */
