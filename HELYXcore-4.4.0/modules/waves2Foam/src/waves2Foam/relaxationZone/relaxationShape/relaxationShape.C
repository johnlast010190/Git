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
    (c) held by original author

\*---------------------------------------------------------------------------*/

#include "relaxationShape.H"

#include "fields/UniformDimensionedFields/uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relaxationShapes
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(relaxationShape, 0);
defineRunTimeSelectionTable(relaxationShape, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

relaxationShape::relaxationShape
(
    const word& subDictName,
    const fvMesh& mesh
)
:
    IOdictionary
    (
        mesh.thisDb().lookupObject<IOobject>("waveProperties")
    ),

    mesh_(mesh),
    coeffDict_(subDict(subDictName + "Coeffs").subDict("relaxationZone")),
    PI_( M_PI ),

    refreshIndexCells_(-1),

    refreshIndexSigma_(-1)
{
    vector g( uniformDimensionedVectorField
        (
            mesh_.thisDb()
            .lookupObject<uniformDimensionedVectorField>("g")).value()
        );

    direction_ = g/mag(g);
}


relaxationShape::~relaxationShape()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void relaxationShape::refreshCells()
{
    if (refreshIndexCells_ != mesh_.time().timeIndex() && mesh_.changing())
    {
        findComputationalCells();
    }

    refreshIndexCells_ = mesh_.time().timeIndex();
}


void relaxationShape::refreshSigma()
{
    refreshCells();

    if (refreshIndexSigma_ != mesh_.time().timeIndex())
    {
        computeSigmaCoordinate();
    }

    refreshIndexSigma_ = mesh_.time().timeIndex();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace relaxationShapes
} // End namespace Foam

// ************************************************************************* //
