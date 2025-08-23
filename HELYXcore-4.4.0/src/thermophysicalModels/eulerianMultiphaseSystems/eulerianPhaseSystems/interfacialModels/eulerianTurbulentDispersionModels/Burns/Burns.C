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
    (c) 2014-2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "Burns.H"
#include "../../../eulerianPhasePair/eulerianPhasePair/eulerianPhasePair.H"
#include "../../../eulerianPhaseModel/MovingPhaseModel/phaseCompressibleTurbulenceModel.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

#include "../../eulerianDragModels/eulerianDragModel/eulerianDragModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace eulerianTurbulentDispersionModels
{
    defineTypeNameAndDebug(Burns, 0);
    addToRunTimeSelectionTable
    (
        eulerianTurbulentDispersionModel,
        Burns,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eulerianTurbulentDispersionModels::Burns::Burns
(
    const dictionary& dict,
    const eulerianPhasePair& pair
)
:
    eulerianTurbulentDispersionModel(dict, pair),
    sigma_("sigma", dimless, dict),
    residualAlpha_
    (
        "residualAlpha",
        dimless,
        dict.lookupOrDefault<scalar>
        (
            "residualAlpha",
            pair_.dispersed().residualAlpha().value()
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::eulerianTurbulentDispersionModels::Burns::~Burns()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::eulerianTurbulentDispersionModels::Burns::D() const
{
    const fvMesh& mesh(pair_.phase1().mesh());
    const eulerianDragModel&
        drag
        (
            mesh.lookupObject<eulerianDragModel>
            (
                IOobject::groupName(eulerianDragModel::typeName, pair_.name())
            )
        );

    return
        0.75
       *drag.CdRe()
       *pair_.continuous().nu()
       *pair_.continuous().turbulence().nut()
       /(
            sigma_
           *sqr(pair_.dispersed().d())
        )
       *pair_.continuous().rho()
       *(
            1.0 + pair_.dispersed().volFrac()
           /max(pair_.continuous().volFrac(), residualAlpha_)
        );
}


// ************************************************************************* //
