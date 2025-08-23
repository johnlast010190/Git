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

#include "Gosman.H"
#include "../../../eulerianPhasePair/eulerianPhasePair/eulerianPhasePair.H"
#include "../../../eulerianPhaseModel/MovingPhaseModel/phaseCompressibleTurbulenceModel.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

#include "../../eulerianDragModels/eulerianDragModel/eulerianDragModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace eulerianTurbulentDispersionModels
{
    defineTypeNameAndDebug(Gosman, 0);
    addToRunTimeSelectionTable
    (
        eulerianTurbulentDispersionModel,
        Gosman,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eulerianTurbulentDispersionModels::Gosman::Gosman
(
    const dictionary& dict,
    const eulerianPhasePair& pair
)
:
    eulerianTurbulentDispersionModel(dict, pair),
    sigma_("sigma", dimless, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::eulerianTurbulentDispersionModels::Gosman::~Gosman()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::eulerianTurbulentDispersionModels::Gosman::D() const
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
       *pair_.dispersed().volFrac()
       *pair_.continuous().nu()
       *pair_.continuous().turbulence().nut()
       /(
            sigma_
           *sqr(pair_.dispersed().d())
        )
       *pair_.continuous().rho();
}


// ************************************************************************* //
