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
    (c) 2014-2019 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "noBlending.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace eulerianBlendingMethods
{
    defineTypeNameAndDebug(noBlending, 0);

    addToRunTimeSelectionTable
    (
        eulerianBlendingMethod,
        noBlending,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eulerianBlendingMethods::noBlending::noBlending
(
    const dictionary& dict,
    const wordList& phaseNames
)
:
    eulerianBlendingMethod(dict),
    continuousPhase_(dict.lookup("continuousPhase"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::eulerianBlendingMethods::noBlending::~noBlending()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::eulerianBlendingMethods::noBlending::f1
(
    const eulerianPhaseModel& phase1,
    const eulerianPhaseModel& phase2
) const
{
    const fvMesh& mesh = phase1.mesh();
    return volScalarField::New
    (
        "f",
        mesh,
        dimensionedScalar(dimless, phase2.name() == continuousPhase_)
    );
}


Foam::tmp<Foam::volScalarField> Foam::eulerianBlendingMethods::noBlending::f2
(
    const eulerianPhaseModel& phase1,
    const eulerianPhaseModel& phase2
) const
{
    const fvMesh& mesh = phase1.mesh();
    return volScalarField::New
    (
        "f",
        mesh,
        dimensionedScalar(dimless, phase1.name() == continuousPhase_)
    );
}


// ************************************************************************* //
