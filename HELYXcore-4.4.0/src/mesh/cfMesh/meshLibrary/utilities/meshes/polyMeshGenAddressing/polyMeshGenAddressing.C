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
    (c) Creative Fields, Ltd.

Authors
    Franjo Juretic (franjo.juretic@c-fields.com)

\*---------------------------------------------------------------------------*/

#include "utilities/meshes/polyMeshGenAddressing/polyMeshGenAddressing.H"
#include "include/demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(polyMeshGenAddressing, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

polyMeshGenAddressing::polyMeshGenAddressing(const polyMeshGenCells& mesh)
:
    mesh_(mesh),
    edgesPtr_(nullptr),
    ccPtr_(nullptr),
    ecPtr_(nullptr),
    pcPtr_(nullptr),
    efPtr_(nullptr),
    pfPtr_(nullptr),

    cePtr_(nullptr),
    fePtr_(nullptr),
    pePtr_(nullptr),
    ppPtr_(nullptr),
    cpPtr_(nullptr),
    cellCentresPtr_(nullptr),
    faceCentresPtr_(nullptr),
    cellVolumesPtr_(nullptr),
    faceAreasPtr_(nullptr),
    globalPointLabelPtr_(nullptr),
    globalFaceLabelPtr_(nullptr),
    globalCellLabelPtr_(nullptr),
    globalEdgeLabelPtr_(nullptr),
    pProcsPtr_(nullptr),
    globalToLocalPointAddressingPtr_(nullptr),
    pointNeiProcsPtr_(nullptr),
    eProcsPtr_(nullptr),
    globalToLocalEdgeAddressingPtr_(nullptr),
    edgeNeiProcsPtr_(nullptr)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyMeshGenAddressing::~polyMeshGenAddressing()
{
    clearAll();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
