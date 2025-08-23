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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2024-2025 Engys Ltd

\*---------------------------------------------------------------------------*/

#include "FSD/reactionRateFlameAreaModels/reactionRateFlameArea/reactionRateFlameArea.H"
#include "fvSolutionRegistry/fvSolutionRegistry.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(reactionRateFlameArea, 0);
    defineRunTimeSelectionTable(reactionRateFlameArea, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reactionRateFlameArea::reactionRateFlameArea
(
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr,
    const combustionModel& combModel
)
:
    coeffDict_(dict.optionalSubDict(modelType + "Coeffs")),
    obr_(obr),
    mesh_(fvSolutionRegistry::getMesh(obr)),
    combModel_(combModel),
    fuel_(dict.lookup("fuel")),
    omega_
    (
        IOobject
        (
            "FSDomega",
            mesh_.time().timeName(),
            obr_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::reactionRateFlameArea::~reactionRateFlameArea()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::reactionRateFlameArea::read(const dictionary& dict)
{
    fuel_ = dict.lookup<word>("fuel");

    return true;
}


// ************************************************************************* //
