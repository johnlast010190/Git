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
    (c) 2011-2022 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "basicChemistryModel/basicChemistryModel.H"
#include "fvMesh/fvMesh.H"
#include "db/Time/Time.H"
#include "fvSolutionRegistry/fvSolutionRegistry.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
    template<>
    const char* NamedEnum<basicChemistryModel::jacobianType, 2>::names[] =
    {
        "fast",
        "exact"
    };
}


const Foam::NamedEnum
<
    Foam::basicChemistryModel::jacobianType,
    2
> Foam::basicChemistryModel::jacobianTypeNames_;


namespace Foam
{
    defineTypeNameAndDebug(basicChemistryModel, 0);
    defineRunTimeSelectionTable(basicChemistryModel, thermo);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basicChemistryModel::basicChemistryModel
(
    const fluidMulticomponentThermo& thermo
)
:
    IOdictionary
    (
        IOobject
        (
            thermo.phasePropertyName("chemistryProperties"),
            thermo.db().time().constant(),
            thermo.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    mesh_(fvSolutionRegistry::getMesh(thermo.db())),
    obr_(thermo.db()),
    thermo_(&thermo),
    solidThermo_(nullptr),
    chemistry_(lookup("chemistry")),
    deltaTChemIni_(lookup<scalar>("initialChemicalTimeStep")),
    deltaTChemMax_(lookupOrDefault("maxChemicalTimeStep", GREAT)),
    deltaTChem_
    (
        IOobject
        (
            thermo.phasePropertyName("deltaTChem"),
            mesh().time().constant(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimTime, deltaTChemIni_)
    )
{}


Foam::basicChemistryModel::basicChemistryModel
(
    const solidMulticomponentThermo& thermo
)
:
    IOdictionary
    (
        IOobject
        (
            thermo.phasePropertyName("chemistryProperties"),
            thermo.db().time().constant(),
            thermo.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    mesh_(fvSolutionRegistry::getMesh(thermo.db())),
    obr_(thermo.db()),
    thermo_(nullptr),
    solidThermo_(&thermo),
    chemistry_(lookup("chemistry")),
    deltaTChemIni_(lookup<scalar>("initialChemicalTimeStep")),
    deltaTChemMax_(lookupOrDefault("maxChemicalTimeStep", GREAT)),
    deltaTChem_
    (
        IOobject
        (
            thermo.phasePropertyName("deltaTChem"),
            mesh().time().constant(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimTime, deltaTChemIni_)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::basicChemistryModel::~basicChemistryModel()
{}


// ************************************************************************* //
