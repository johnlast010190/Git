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
    (c) 2011-2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "isothermalDiameter.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace eulerianDiameterModels
{
    defineTypeNameAndDebug(isothermal, 0);

    addToRunTimeSelectionTable
    (
        eulerianDiameterModel,
        isothermal,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eulerianDiameterModels::isothermal::isothermal
(
    const dictionary& diameterProperties,
    const eulerianPhaseModel& phase
)
:
    eulerianDiameterModel(diameterProperties, phase),
    d0_("d0", dimLength, diameterProperties_),
    p0_("p0", dimPressure, diameterProperties_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::eulerianDiameterModels::isothermal::~isothermal()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::eulerianDiameterModels::isothermal::d() const
{
    const volScalarField& p =
        phase_.volFrac().db().lookupObject<volScalarField>("p");
    dimensionedScalar pRef = basicThermo::lookupOrCreate(p.db()).pRef();

    return d0_*pow((p0_+pRef)/(p+pRef), 1.0/3.0);
}


bool Foam::eulerianDiameterModels::isothermal::read(const dictionary& phaseProperties)
{
    eulerianDiameterModel::read(phaseProperties);

    diameterProperties_.lookup("d0") >> d0_;
    diameterProperties_.lookup("p0") >> p0_;

    return true;
}


// ************************************************************************* //
