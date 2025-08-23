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

#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "autoOptimize/helyxOptimize/helyxOptimize.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(helyxOptimize, 0);
    addToRunTimeSelectionTable
    (
        autoOptimize,
        helyxOptimize,
        dictionary
    );
}


void Foam::helyxOptimize::movePoints(Foam::pointField& newPoints)
{
    optimPtr_.reset
    (
        new hessianMeshOptimization
        (mesh_,newPoints, coeffsDict_, true)
    );
    newPoints = optimPtr_->newPoints();
}


void Foam::helyxOptimize::optimize()
{
    optimPtr_.reset
    (
        new hessianMeshOptimization
        (mesh_, coeffsDict_, true)
    );
    tmp<pointField> newPoints = optimPtr_->newPoints();

    mesh_.movePoints(newPoints);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::helyxOptimize::helyxOptimize
(
    fvMesh& mesh,
    const dictionary& dict
)
:
    autoOptimize(mesh, dict),
    mesh_(mesh),
    coeffsDict_(dict.subDict("helyxOptimizeCoeffs")),
    optimPtr_(nullptr)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::helyxOptimize::~helyxOptimize()
{}

// ************************************************************************* //
