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
    (c) held by original author

\*---------------------------------------------------------------------------*/

#include "porosityCoefficient.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(porosityCoefficient, 0);
defineRunTimeSelectionTable(porosityCoefficient, porosityCoefficient);

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //


// This method has been added, such that it is possible to have porosity and
// moving meshes at the same time. If resistancePorosity is applied in the
// input dictionary, it is simply possible to apply "porosity 1.0" for the
// hydrodynamic calculations. This option is required for boundedness issues
// in MULES (tested in 1.6-ext and foam-extend-3.1)
// NGJ 16/03/2015
scalar porosityCoefficient::readResistancePorosity
(
	const dictionary& dict
) const
{
    if (dict.found("resistancePorosity"))
    {
    	Info<< "Resistance coefficients are based on the porosity given by the"
    	     << " keyword 'resistancePorosity'.\n" << endl;
        return dict.lookup<scalar>("resistancePorosity");
    }
    else
    {
    	Info<< "Resistance coefficients are based on the porosity given by the"
    	     << " keyword 'porosity'.\n" << endl;
    	return dict.lookup<scalar>("porosity");
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

porosityCoefficient::porosityCoefficient
(
    const dictionary & poroProp
)
:
    poroProperties_(poroProp),

    linearCoefficient_(dimensionedVector(dimless, vector::zero)),

    quadraticCoefficient_(dimensionedVector(dimless, vector::zero)),

    KCQuadraticCoefficient_(dimensionedVector(dimless, vector::zero)),

    scaledKC_(dimensionedScalar(dimless, 0))
{

}


porosityCoefficient::~porosityCoefficient()
{}

autoPtr<porosityCoefficient> porosityCoefficient::New
(
    const dictionary & poroProp
)
{
    word poroForm = poroProp.lookup("resistanceFormulation");

    const auto ctor =
        ctorTableLookup
        (
            "resistance formulation",
            porosityCoefficientConstructorTable_(),
            poroForm
        );
    return autoPtr<porosityCoefficient>(ctor( poroProp ));
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
