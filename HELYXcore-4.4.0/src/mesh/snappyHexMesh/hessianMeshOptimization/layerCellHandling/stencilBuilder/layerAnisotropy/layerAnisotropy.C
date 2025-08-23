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
    (c) 2017 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "hessianMeshOptimization/layerCellHandling/stencilBuilder/layerAnisotropy/layerAnisotropy.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(layerAnisotropy, 0);
    addToRunTimeSelectionTable
    (
        stencilBuilder,
        layerAnisotropy,
        dictionary
    );
}


void Foam::layerAnisotropy::calculateAnisotropy()
{
    forAll(mesh_.cells(), cI)
    {
        if (!layers_.isALayerCell(cI))
        {
            anisotropy_[cI] = 1;
        }
    }
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::layerAnisotropy::layerAnisotropy
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    stencilBuilder(mesh),
    mesh_(mesh),
    layers_(mesh, 1.6),
    anisotropy_(layers_.aspectRatio())
{
    calculateAnisotropy();
}

Foam::layerAnisotropy::~layerAnisotropy()
{}

const Foam::scalarField& Foam::layerAnisotropy::anisotropy() const
{
    return anisotropy_;
}

const Foam::vectorField& Foam::layerAnisotropy::anisotropicDirection() const
{
    return layers_.layerDirection();
}
