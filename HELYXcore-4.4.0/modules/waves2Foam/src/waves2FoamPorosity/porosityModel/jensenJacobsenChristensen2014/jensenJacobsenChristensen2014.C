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

#include "jensenJacobsenChristensen2014.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(jensenJacobsenChristensen2014, 0);
addToRunTimeSelectionTable
(
    wavesPorosityModel,
    jensenJacobsenChristensen2014,
    wavesPorosityModel
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


jensenJacobsenChristensen2014::jensenJacobsenChristensen2014
(
    const fvMesh& mesh
)
:
    wavesPorosityModel(mesh),

    pZones_(mesh)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<fvMatrix<vector>> jensenJacobsenChristensen2014::ddt
(
    VolField<vector>& U
)
{
    return pZones_.ddt(U);
}


tmp<fvMatrix<vector>> jensenJacobsenChristensen2014::ddt
(
    const geometricOneField& rho,
    VolField<vector>& U
)
{
    return pZones_.ddt(rho, U);
}


tmp<fvMatrix<vector>> jensenJacobsenChristensen2014::ddt
(
    const dimensionedScalar& rho,
    VolField<vector>& U
)
{
	return pZones_.ddt(rho, U);
}


tmp<fvMatrix<vector>> jensenJacobsenChristensen2014::ddt
(
    const volScalarField& rho,
    VolField<vector>& U
)
{
    return pZones_.ddt(rho, U);
}

void jensenJacobsenChristensen2014::updatePorosity()
{
	// Store the porosity from the last time step. Needed for moving porosity
	// fields
	porosity_.storeOldTimes();

	// Obtain the new porosity field from the pZones as a tmp<volScalarField>
	tmp<volScalarField> tporosity = pZones_.porosity();
	const volScalarField& poro = tporosity();

    // Set the internal field values
#if EXTBRANCH==1
	porosity_.internalField() = poro.internalField();
#elif OFPLUSBRANCH==1
    #if OFVERSION<1706
        porosity_.internalField() = poro.internalField();
    #else
        porosity_.ref() = poro.internalField();
    #endif
#else
    #if OFVERSION<400
	porosity_.internalField() = poro.internalField();
    #else
	porosity_.ref() = poro.internalField();
    #endif
#endif

	// Update boundary conditions
	porosity_.correctBoundaryConditions();
}


const volScalarField& jensenJacobsenChristensen2014::porosity() const
{
//	// Store the porosity from the last time step. Needed for moving porosity
//	// fields
////	porosity_.storeOldTimes();
//
//	tmp<volScalarField> tporosity = pZones_.porosity();
//
//
//
////    return pZones_.porosity();
////    return tporosity;
    return porosity_;
}


void jensenJacobsenChristensen2014::addResistance(fvVectorMatrix& UEqn) const
{
    pZones_.addResistance(UEqn);
}


void jensenJacobsenChristensen2014::addResistance
(
    const fvVectorMatrix& UEqn,
    volTensorField& AU
) const
{
    pZones_.addResistance(UEqn, AU);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
