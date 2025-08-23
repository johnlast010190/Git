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
    (c) 2011-2012 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "motionDiffusivity/manipulators/quadratic/quadraticDiffusivity.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(quadraticDiffusivity, 0);

    addToRunTimeSelectionTable
    (
        motionDiffusivity,
        quadraticDiffusivity,
        Istream
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::quadraticDiffusivity::quadraticDiffusivity
(
    const fvMesh& mesh,
    Istream& mdData
)
:
    motionDiffusivity(mesh),
    basicDiffusivityPtr_(motionDiffusivity::New(mesh, mdData))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::quadraticDiffusivity::~quadraticDiffusivity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField>
Foam::quadraticDiffusivity::operator()() const
{
    return sqr(basicDiffusivityPtr_->operator()());
}


void Foam::quadraticDiffusivity::correct()
{
    basicDiffusivityPtr_->correct();
}


// ************************************************************************* //
