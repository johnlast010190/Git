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

#include "constantTurbulentDispersionCoefficient.H"
#include "../../../eulerianPhasePair/eulerianPhasePair/eulerianPhasePair.H"
#include "../../../eulerianPhaseModel/MovingPhaseModel/phaseCompressibleTurbulenceModel.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace eulerianTurbulentDispersionModels
{
    defineTypeNameAndDebug(constantTurbulentDispersionCoefficient, 0);
    addToRunTimeSelectionTable
    (
        eulerianTurbulentDispersionModel,
        constantTurbulentDispersionCoefficient,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eulerianTurbulentDispersionModels::constantTurbulentDispersionCoefficient::
constantTurbulentDispersionCoefficient
(
    const dictionary& dict,
    const eulerianPhasePair& pair
)
:
    eulerianTurbulentDispersionModel(dict, pair),
    Ctd_("Ctd", dimless, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::eulerianTurbulentDispersionModels::constantTurbulentDispersionCoefficient::
~constantTurbulentDispersionCoefficient()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::eulerianTurbulentDispersionModels::constantTurbulentDispersionCoefficient::
D() const
{
    return
        Ctd_
       *pair_.dispersed().volFrac()
       *pair_.continuous().rho()
       *pair_.continuous().turbulence().k();
}


// ************************************************************************* //
