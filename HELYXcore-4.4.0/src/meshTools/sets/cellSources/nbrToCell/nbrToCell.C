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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2024 Engys Ltd.

\*---------------------------------------------------------------------------*/

#include "sets/cellSources/nbrToCell/nbrToCell.H"
#include "meshes/polyMesh/polyMesh.H"

#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(nbrToCell, 0);

addToRunTimeSelectionTable(topoSetSource, nbrToCell, word);

addToRunTimeSelectionTable(topoSetSource, nbrToCell, istream);

}


Foam::topoSetSource::addToUsageTable Foam::nbrToCell::usage_
(
    nbrToCell::typeName,
    "\n    Usage: nbrToCell <nNeighbours>\n\n"
    "    Select all cells with <= nNeighbours neighbouring cells\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::nbrToCell::combine(topoSet& set, const bool add) const
{
    const cellList& cells = mesh().cells();
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    boolList isCoupled(mesh_.nFaces()-mesh_.nInternalFaces(), false);

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        if (pp.coupled())
        {
            label facei = pp.start();
            forAll(pp, i)
            {
                isCoupled[facei-mesh_.nInternalFaces()] = true;
                facei++;
            }
        }
    }

    forAll(cells, celli)
    {
        const cell& cFaces = cells[celli];

        label nNbrCells = 0;

        forAll(cFaces, i)
        {
            label facei = cFaces[i];

            if (mesh_.isInternalFace(facei))
            {
                nNbrCells++;
            }
            else if (isCoupled[facei-mesh_.nInternalFaces()])
            {
                nNbrCells++;
            }
        }

        if (nNbrCells <= minNbrs_)
        {
            addOrDelete(set, celli, add);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::nbrToCell::nbrToCell
(
    const polyMesh& mesh,
    const label minNbrs
)
:
    topoSetSource(mesh),
    minNbrs_(minNbrs)
{}


// Construct from dictionary
Foam::nbrToCell::nbrToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    minNbrs_(dict.lookup<label>("neighbours"))
{}


// Construct from Istream
Foam::nbrToCell::nbrToCell
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetSource(mesh),
    minNbrs_(readLabel(checkIs(is)))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nbrToCell::~nbrToCell()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::nbrToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
    {
        Info<< "    Adding cells with only " << minNbrs_ << " or less"
                " neighbouring cells" << " ..." << endl;

        combine(set, true);
    }
    else if (action == topoSetSource::DELETE)
    {
        Info<< "    Removing cells with only " << minNbrs_ << " or less"
                " neighbouring cells" << " ..." << endl;

        combine(set, false);
    }
}


// ************************************************************************* //
